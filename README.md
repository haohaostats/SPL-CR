# SPL-CR <img src="https://img.shields.io/badge/R-package-276DC3?logo=r&logoColor=white" alt="R package" align="right"/>

<p align="center">
  <b>Smoothed Predictive Likelihood for Competing Risks</b><br>
  <i>Likelihood-based validation for cause-specific and Fine–Gray prediction models</i>
</p>

<p align="center">
  <a href="https://github.com/haohaostats/SPL-CR"><img src="https://img.shields.io/badge/GitHub-haohaostats%2FSPL--CR-black?logo=github" alt="GitHub repo"></a>
  <img src="https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-4caf50" alt="platform">
  <img src="https://img.shields.io/badge/status-active%20development-0a7ea4" alt="status">
  <img src="https://img.shields.io/badge/methods-SPL--CR%20%7C%20SPL--FG-7b1fa2" alt="methods">
</p>

---

## ✨ Overview

`SPL-CR` provides **likelihood-based, estimand-aware validation tools for competing-risks prediction models**.

It implements two cross-validated smoothed predictive likelihood scores:

- **SPL-CR**: for **cause-specific hazards (CSH)** models (Cox-type per cause)
- **SPL-FG**: for **Fine–Gray (FG)** subdistribution hazard models

The package addresses a key technical issue in survival/competing-risks validation:

> Standard Cox / Fine–Gray baseline estimators are step functions.  
> For out-of-sample events occurring between training jump times, naive predictive densities can be zero, causing \(\log(0)\) failures in cross-validation.

`SPL-CR` resolves this by applying a **post-fit k-nearest-neighbour kernel smoother** to training-fold baseline components only, yielding finite and stable out-of-sample log-scores.

---

## 🔬 Methods implemented

### 1) SPL-CR (Cause-Specific Hazards path)
- Fit one Cox model per cause on each training fold
- Smooth baseline cumulative hazard jumps \( \Delta \widehat{\Lambda}_{0c}(t) \)
- Evaluate out-of-fold predictive log-score under the cause-specific hazard representation

### 2) SPL-FG (Fine–Gray path)
- Fit one Fine–Gray model per cause on each training fold
- Extract baseline CIF jumps from `predictRisk()` on a reference covariate profile
- Smooth baseline CIF and induced density
- Evaluate out-of-fold predictive log-score under the Fine–Gray representation

### Shared smoothing engine
- **k-NN adaptive bandwidth**
- **Boundary renormalization**
- Kernels:
  - Uniform
  - Triangular
  - Epanechnikov
  - Quartic (Biweight)

---

## 📦 Installation

### From GitHub (recommended)

```r
# install.packages("remotes")
remotes::install_github("haohaostats/SPL-CR")
```

Then load:

```r
library(splcr)
```

---

## 🧰 Dependencies

Core dependencies include:

- `survival`
- `riskRegression`
- `prodlim`
- `Rcpp`
- `RcppArmadillo`

---

## 🚀 Quick start (PBC example)

```r
library(splcr)
library(survival)

data(pbc, package = "survival")

# PBC preprocessing (competing risks: 1=death, 2=transplant, 0=censoring)
pbc_data <- subset(pbc, !is.na(trt))
pbc_data$status <- ifelse(pbc_data$status == 2, 1L,
                          ifelse(pbc_data$status == 1, 2L, 0L))
pbc_data$time <- as.numeric(pbc_data$time)  # use original unit (days)
pbc_data$logbili <- log(pbc_data$bili)
pbc_data$sex <- factor(pbc_data$sex)
pbc_data$stage <- factor(pbc_data$stage)
pbc_data$ascites <- factor(pbc_data$ascites)

pbc_data <- pbc_data[complete.cases(
  pbc_data[, c("time","status","logbili","albumin","protime","ascites","age","sex","stage")]
), ]

# SPL-CR
res_cr <- spl_cr(
  ~ logbili + albumin + protime + ascites,
  data   = pbc_data,
  folds  = 10,
  k_nn   = 5,
  kernel = "epanechnikov",
  seed   = 20250101
)
print(res_cr)

# SPL-FG
res_fg <- spl_fg(
  ~ logbili + albumin + protime + ascites,
  data   = pbc_data,
  folds  = 10,
  k_nn   = 5,
  kernel = "epanechnikov",
  seed   = 20250101
)
print(res_fg)
```

---

## 📊 Compare multiple model specifications

Use `spl_compare()` to compare candidate formulas across kernels under a shared fold assignment.

```r
models <- list(
  "M1: Baseline Clinical"    = ~ age + sex + stage,
  "M2: Core Biochemical"     = ~ logbili + albumin,
  "M3: Extended Biochemical" = ~ logbili + albumin + protime + ascites
)

# Compare under SPL-CR
tab_cr <- spl_compare(
  models  = models,
  data    = pbc_data,
  method  = "cr",
  kernels = c("uniform", "triangular", "epanechnikov", "quartic"),
  folds   = 10,
  k_nn    = 5,
  seed    = 20250101
)

# Compare under SPL-FG
tab_fg <- spl_compare(
  models  = models,
  data    = pbc_data,
  method  = "fg",
  kernels = c("uniform", "triangular", "epanechnikov", "quartic"),
  folds   = 10,
  k_nn    = 5,
  seed    = 20250101
)

# Convert to manuscript-style wide tables (kernels as columns)
tab_cr_wide <- spl_compare_wide(tab_cr)
tab_fg_wide <- spl_compare_wide(tab_fg)

# Pretty print
print(tab_cr_wide)
print(tab_fg_wide)
```

---

## ⚠️ Important notes (time unit & comparability)

### Time unit is **user-defined** (days / years / months all allowed)
Both `spl_cr()` and `spl_fg()` are **time-unit agnostic**.

However:

- Absolute log-score values **depend on the time unit**
- For fair model comparison, **all models must use the same time unit**
- The manuscript real-data examples (e.g., PBC / melanoma) use **time in days**

### Fine–Gray implementation note (manuscript-consistent behavior)
The package implementation of `spl_fg()` uses a **common all-cause event-time grid within each training fold** when extracting baseline CIF jumps, matching the manuscript reference scripts.

---

## 📁 Repository structure

```text
.
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── spl_cr.R
│   ├── spl_fg.R
│   ├── spl_compare.R
│   ├── utils.R
│   └── ...
├── src/
│   └── khs.cpp
├── man/               
├── README.md
└── ...
```

---

## 📬 Contact

- **Hao Chen** — University of Sydney  
- GitHub: <https://github.com/haohaostats>

---

<p align="center">
  <i>If this package is useful in your work, a ⭐ on the repository helps visibility and reproducibility.</i>
</p>
