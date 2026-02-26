# Internal wrappers for registered C++ routines

.calculate_smoothed_hazard_cpp <- function(Jumps_H0, bw, nn, t, method_id) {
  .Call(`_splcr_calculate_smoothed_hazard`, Jumps_H0, bw, as.integer(nn), t, as.integer(method_id))
}

.smooth_level_and_density_cpp <- function(Jumps_H, nn, t, kernel_id) {
  .Call(`_splcr_smooth_level_and_density`, Jumps_H, as.integer(nn), t, as.integer(kernel_id))
}
