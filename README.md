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

- **SPL-CSH**: for **cause-specific hazards (CSH)** models (Cox-type per cause)
- **SPL-FG**: for **Fine–Gray (FG)** subdistribution hazard models

The package addresses a key technical issue in survival/competing-risks validation:

> Standard Cox / Fine–Gray baseline estimators are step functions.  
> For out-of-sample events occurring between training jump times, naive predictive densities can be zero, causing `log(0)` failures in cross-validation.

`SPL-CR` resolves this by applying a **post-fit k-nearest-neighbour kernel smoother** to training-fold baseline components only, yielding finite and stable out-of-sample log-scores.

---

## 🔬 Methods implemented

### 1) SPL-CSH (Cause-Specific Hazards path)
- Fit one Cox model per cause on each training fold
- Smooth baseline cumulative hazard jumps \( \Delta \widehat{\Lambda}_{0c}(t) \)
- Evaluate out-of-fold predictive log-score under the cause-specific hazards representation

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

# SPL-CSH
res_csh <- spl_csh(
  ~ logbili + albumin + protime + ascites,
  data   = pbc_data,
  folds  = 10,
  k_nn   = 5,
  kernel = "epanechnikov",
  seed   = 20250101
)
print(res_csh)

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

**Example output (PBC cohort):**

```text
Smoothed Predictive Likelihood (SPL-CR)
Kernel: epanechnikov  kNN: 5  folds: 10
CV mean log-score (Lhat): -3.9810
Jackknife SE: 0.0358
95% CI: [-4.0619, -3.9000]

Smoothed Predictive Likelihood (SPL-FG)
Kernel: epanechnikov  kNN: 5  folds: 10
CV mean log-score (Lhat): -3.7391
Jackknife SE: 0.0389
95% CI: [-3.8270, -3.6512]
```

---

## 📊 Compare multiple model specifications

Use `spl_compare()` to compare candidate formulas across kernels under a shared fold assignment.

```r
models <- list(
  "M1: Baseline Clinical"    = ~ age + sex + stage,
  "M2: Core Biochemical"     = ~ logbili + albumin,
  "M3: Extended Biochemical" = ~ logbili + albumin + protime + ascites,
  "M4: Clinical + Core Bio"  = ~ age + sex + stage + logbili + albumin,
  "M5: Full"                 = ~ age + sex + stage + logbili + albumin + protime + ascites
)

# Compare under SPL-CR
tab_csh <- spl_compare(
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
tab_csh_wide <- spl_compare_wide(tab_csh)
tab_fg_wide <- spl_compare_wide(tab_fg)

# Pretty print
print(tab_csh_wide)
print(tab_fg_wide)
```

**Example output (wide tables, `metric = Lhat`):**

```text
📊 SPL model comparison  •  metric = Lhat
────────────────────────────────────────────────────────────────────────────
Model                        Uniform   Triangular   Epanechnikov     Quartic
────────────────────────────────────────────────────────────────────────────
M1: Baseline Clinical        -4.6304      -4.6792        -4.6735     -4.6807
M2: Core Biochemical       ★ -2.6720    ★ -2.7210      ★ -2.7155   ★ -2.7222
M3: Extended Biochemical     -3.9374      -3.9865        -3.9810     -3.9876
M4: Clinical + Core Bio      -3.1065      -3.1565        -3.1511     -3.1579
M5: Full                     -4.0118      -4.0620        -4.0568     -4.0634
────────────────────────────────────────────────────────────────────────────
Note: ★ marks the highest value in each kernel column (larger is better).

📊 SPL model comparison  •  metric = Lhat
────────────────────────────────────────────────────────────────────────────
Model                        Uniform   Triangular   Epanechnikov     Quartic
────────────────────────────────────────────────────────────────────────────
M1: Baseline Clinical        -4.1256      -4.1721        -4.1677     -4.1732
M2: Core Biochemical         -4.3483      -4.4351        -4.4223     -4.4452
M3: Extended Biochemical   ★ -3.6968    ★ -3.7459      ★ -3.7391   ★ -3.7545
M4: Clinical + Core Bio      -3.8935      -3.9405        -3.9358     -3.9421
M5: Full                     -4.1578      -4.2067        -4.1996     -4.2121
────────────────────────────────────────────────────────────────────────────
Note: ★ marks the highest value in each kernel column (larger is better).
```

---

## ⚠️ Important notes (time unit & comparability)

### Time unit is **user-defined** (days / years / months all allowed)
Both `spl_csh()` and `spl_fg()` are **time-unit agnostic**.

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
│   ├── spl_csh.R
│   ├── spl_fg.R
│   ├── spl_compare.R
│   ├── spl_compare_wide.R
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
