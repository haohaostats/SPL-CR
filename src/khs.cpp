// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ============================================================================
//                               Kernel family
// ============================================================================

inline double evaluate_kernel(double u, int kernel_id) {
  switch (kernel_id) {
  case 1: // Uniform
    return (std::abs(u) <= 1.0) ? 0.5 : 0.0;
  case 2: // Triangular
    return (std::abs(u) <= 1.0) ? (1.0 - std::abs(u)) : 0.0;
  case 3: // Epanechnikov
    return (std::abs(u) <= 1.0) ? (0.75 * (1.0 - u*u)) : 0.0;
  case 4: { // Quartic (Biweight)
      if (std::abs(u) > 1.0) return 0.0;
      double t = 1.0 - u*u;
      return 0.9375 * t * t;  // 15/16
    }
  default:
    return 0.0;
  }
}

inline double K_primitive(double z, int kernel_id) {
  if (z <= -1.0) z = -1.0;
  if (z >=  1.0) z =  1.0;

  switch (kernel_id) {
  case 1:
    return 0.5 * z;

  case 2:
    if (z >= 0) {
      return z - 0.5 * z * z;
    } else {
      return z + 0.5 * z * z;
    }

  case 3:
    return 0.75 * (z - (z*z*z)/3.0);

  case 4:
    return (15.0/16.0) * ( z - (2.0/3.0)*z*z*z + (1.0/5.0)*std::pow(z,5) );

  default:
    return 0.0;
  }
}

inline double kernel_overlap_integral(int kernel_id, double L, double U) {
  if (L > U) std::swap(L, U);
  double low  = std::max(-1.0, L);
  double high = std::min( 1.0, U);
  if (high <= low) return 0.0;
  return K_primitive(high, kernel_id) - K_primitive(low, kernel_id);
}

// ============================================================================
//                 kNN kernel smoothing for baseline hazard jumps
// ============================================================================
inline double smooth_hazard_knn_bandwidth(arma::mat hazard_jumps,
                                          int n_neighbors,
                                          double time_point,
                                          int kernel_id) {
  if (time_point < 0) return 0.0;
  int num_jumps = hazard_jumps.n_rows;
  if (num_jumps == 0) return 0.0;

  double tau_max = hazard_jumps(num_jumps - 1, 0);

  if (time_point >= hazard_jumps(num_jumps - 1, 0)) {
    arma::rowvec new_row(2); new_row(0) = time_point; new_row(1) = 0.0;
    hazard_jumps.insert_rows(num_jumps, new_row);
    num_jumps++;
  }

  int nearest_left_idx = -1;
  for (int i = 0; i < num_jumps; ++i) {
    if (hazard_jumps(i,0) < time_point) nearest_left_idx = i; else break;
  }
  int left_k_idx = -1, found = 0;
  for (int i = nearest_left_idx; i >= 0; --i) {
    if (hazard_jumps(i,1) > 0) { found++; if (found >= n_neighbors) { left_k_idx = i; break; } }
  }

  int nearest_right_idx = num_jumps - 1;
  for (int i = 0; i < num_jumps; ++i) {
    if (hazard_jumps(i,0) > time_point) { nearest_right_idx = i; break; }
  }
  int right_k_idx = -1; found = 0;
  for (int i = nearest_right_idx; i < num_jumps; ++i) {
    if (hazard_jumps(i,1) > 0) { found++; if (found >= n_neighbors) { right_k_idx = i; break; } }
  }

  double tL = (left_k_idx  >= 0) ? hazard_jumps(left_k_idx,  0) : 0.0;
  double tR = (right_k_idx >= 0) ? hazard_jumps(right_k_idx, 0) : tau_max;

  double left_dist  = std::max(0.0, time_point - tL);
  double right_dist = std::max(0.0, tR - time_point);
  double bandwidth  = std::max(left_dist, right_dist);
  if (bandwidth <= 0) return 0.0;

  double winL = std::max(0.0, time_point - bandwidth);
  double winR = std::min(tau_max, time_point + bandwidth);

  int iL = 0; while (iL < num_jumps && hazard_jumps(iL,0) < winL) ++iL;
  int iR = iL; while (iR < num_jumps && hazard_jumps(iR,0) <= winR) ++iR; --iR;
  if (iR < iL) return 0.0;

  double numer = 0.0;
  for (int i = iL; i <= iR; ++i) {
    double u = (time_point - hazard_jumps(i,0)) / bandwidth;
    double w = evaluate_kernel(u, kernel_id);
    if (w != 0.0) numer += w * hazard_jumps(i,1);
  }

  double L = (0.0     - time_point) / bandwidth;
  double U = (tau_max - time_point) / bandwidth;
  double c = kernel_overlap_integral(kernel_id, L, U);
  const double EPS = 1e-12;
  if (c < EPS) return 0.0;

  double lambda_s = (numer / bandwidth) / c;
  return (lambda_s > 0.0) ? lambda_s : 0.0;
}

// C++ implementation used by the .Call wrapper
static double calc_smoothed_hazard_cpp(Rcpp::NumericMatrix Jumps_H0,
                                       double bw, int nn, double t, int method_id) {
  arma::mat hazard_jumps(Jumps_H0.begin(), Jumps_H0.nrow(), Jumps_H0.ncol(), false);
  if (bw > 0) {
    Rcpp::stop("Fixed bandwidth method not implemented here. Use k-NN (bw <= 0).");
    return NA_REAL;
  }
  return smooth_hazard_knn_bandwidth(hazard_jumps, nn, t, method_id);
}

// ============================================================================
//        kNN kernel smoothing for baseline step function H(t) (e.g., F0)
// ============================================================================
inline std::pair<double,double>
smooth_level_and_density_knn(const arma::mat& jumps,
                             int n_neighbors,
                             double time_point,
                             int kernel_id) {
  const double EPS = 1e-12;

  int J = jumps.n_rows;
  if (J == 0 || time_point < 0.0) return {0.0, 0.0};

  double tau_max = jumps(J - 1, 0);

  int left_k_idx = -1, found = 0;
  for (int i = J - 1; i >= 0; --i) {
    if (jumps(i,0) < time_point && jumps(i,1) > 0) { found++; if (found >= n_neighbors) { left_k_idx = i; break; } }
  }
  int right_k_idx = -1; found = 0;
  for (int i = 0; i < J; ++i) {
    if (jumps(i,0) > time_point && jumps(i,1) > 0) { found++; if (found >= n_neighbors) { right_k_idx = i; break; } }
  }

  double tL = (left_k_idx  >= 0) ? jumps(left_k_idx,  0) : 0.0;
  double tR = (right_k_idx >= 0) ? jumps(right_k_idx, 0) : tau_max;

  double left_dist  = std::max(0.0, time_point - tL);
  double right_dist = std::max(0.0, tR - time_point);
  double bw = std::max(left_dist, right_dist);
  if (bw <= 0.0) return {0.0, 0.0};

  double L = (0.0     - time_point) / bw;
  double U = (tau_max - time_point) / bw;
  double Z = kernel_overlap_integral(kernel_id, L, U);
  if (Z < EPS) return {0.0, 0.0};

  double winL = std::max(0.0, time_point - bw);
  double winR = std::min(tau_max, time_point + bw);

  int iL = 0; while (iL < J && jumps(iL,0) < winL) ++iL;
  int iR = iL; while (iR < J && jumps(iR,0) <= winR) ++iR; --iR;
  if (iR < iL) return {0.0, 0.0};

  double numer_level = 0.0;
  double numer_dens  = 0.0;

  for (int i = iL; i <= iR; ++i) {
    double u = (time_point - jumps(i,0)) / bw;
    double dH = jumps(i,1);
    double Kpdf = evaluate_kernel(u, kernel_id);
    double Kcdf = 0.5 + K_primitive(u, kernel_id);
    numer_dens  += Kpdf * dH;
    numer_level += Kcdf * dH;
  }

  double Hs = numer_level / Z;
  double hs = (numer_dens / bw) / Z;

  if (Hs < 0.0) Hs = 0.0;
  if (Hs > 1.0) Hs = 1.0;
  if (hs < 0.0) hs = 0.0;

  return {Hs, hs};
}

static Rcpp::NumericVector smooth_level_density_cpp(Rcpp::NumericMatrix Jumps_H,
                                                    int nn, double t, int kernel_id) {
  arma::mat J(Jumps_H.begin(), Jumps_H.nrow(), Jumps_H.ncol(), false);
  auto pr = smooth_level_and_density_knn(J, nn, t, kernel_id);
  Rcpp::NumericVector out(2);
  out[0] = pr.first;
  out[1] = pr.second;
  return out;
}

// ---------------------------------------------------------------------------
// .Call wrappers with registration
// ---------------------------------------------------------------------------
extern "C" SEXP _splcr_calculate_smoothed_hazard(SEXP Jumps_H0SEXP,
                                                 SEXP bwSEXP,
                                                 SEXP nnSEXP,
                                                 SEXP tSEXP,
                                                 SEXP method_idSEXP) {
  BEGIN_RCPP;
  NumericMatrix Jumps_H0(Jumps_H0SEXP);
  double bw = as<double>(bwSEXP);
  int nn = as<int>(nnSEXP);
  double t = as<double>(tSEXP);
  int method_id = as<int>(method_idSEXP);
  double res = calc_smoothed_hazard_cpp(Jumps_H0, bw, nn, t, method_id);
  return wrap(res);
  END_RCPP;
}

extern "C" SEXP _splcr_smooth_level_and_density(SEXP Jumps_HSEXP,
                                                SEXP nnSEXP,
                                                SEXP tSEXP,
                                                SEXP kernel_idSEXP) {
  BEGIN_RCPP;
  NumericMatrix Jumps_H(Jumps_HSEXP);
  int nn = as<int>(nnSEXP);
  double t = as<double>(tSEXP);
  int kernel_id = as<int>(kernel_idSEXP);
  NumericVector out = smooth_level_density_cpp(Jumps_H, nn, t, kernel_id);
  return out;
  END_RCPP;
}

static const R_CallMethodDef CallEntries[] = {
  {"_splcr_calculate_smoothed_hazard", (DL_FUNC) &_splcr_calculate_smoothed_hazard, 5},
  {"_splcr_smooth_level_and_density", (DL_FUNC) &_splcr_smooth_level_and_density, 4},
  {NULL, NULL, 0}
};

extern "C" void R_init_splcr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
