# SPL-CR <img src="https://img.shields.io/badge/R-package-276DC3?logo=r&logoColor=white" alt="R package" align="right"/>

<p align="center">
  <img src="https://img.shields.io/badge/Competing%20Risks-Likelihood%20Validation-0f4c5c?style=for-the-badge" alt="Competing Risks Validation">
  <img src="https://img.shields.io/badge/SPL--CSH%20%7C%20SPL--FG-Methodology-c6903d?style=for-the-badge" alt="Methods">
</p>

<p align="center">
  <b>Smoothed Predictive Likelihood for Competing Risks</b><br>
  <i>Estimand-aware cross-validated validation for cause-specific and Fine-Gray prediction models</i>
</p>

<p align="center">
  <a href="https://github.com/haohaostats/SPL-CR"><img src="https://img.shields.io/badge/GitHub-Repository-black?logo=github" alt="GitHub repo"></a>
  <a href="https://haostats.shinyapps.io/splcr-web/"><img src="https://img.shields.io/badge/Web-App%20Live-0f766e?logo=shiny&logoColor=white" alt="Web app"></a>
  <img src="https://img.shields.io/badge/platform-Windows%20%7C%20macOS%20%7C%20Linux-4caf50" alt="platform">
  <img src="https://img.shields.io/badge/status-active%20development-0a7ea4" alt="status">
</p>

---

## Why SPL-CR

A core challenge in competing-risks validation is that baseline estimators from Cox and Fine-Gray models are step functions.
Out-of-sample event times falling between training jump points can yield zero predictive densities and unstable `log(0)` terms.

`SPL-CR` resolves this through post-fit k-nearest-neighbour smoothing of training-fold baseline components only, yielding finite and stable predictive log-scores while preserving estimand alignment.

---

## Methods Implemented

### SPL-CSH

- One Cox model per cause in each training fold
- Smoothing of baseline cumulative hazard jumps \(\Delta \widehat{\Lambda}_{0c}(t)\)
- Out-of-fold predictive log-score under the cause-specific hazards representation

### SPL-FG

- One Fine-Gray model per cause in each training fold
- Baseline CIF jump extraction via `predictRisk()` on a reference profile
- Smoothing of baseline CIF and induced event density
- Out-of-fold predictive log-score under the Fine-Gray representation

### Shared smoothing engine

- k-NN adaptive bandwidth
- boundary renormalization
- supported kernels: `uniform`, `triangular`, `epanechnikov`, `quartic`

---

## Clinical-Research Oriented Access

`SPL-CR` now supports two parallel workflows:

- **Interactive Web App (recommended for applied users):** run analyses in a browser with guided inputs, status-label mapping, and formatted result tables.
- **R Package (recommended for methodological users):** fully scriptable workflows for reproducible analysis, simulation, and manuscript pipelines.

**Web app:** [https://haostats.shinyapps.io/splcr-web/](https://haostats.shinyapps.io/splcr-web/)

![SPL-CR Web App Screenshot](https://image.thum.io/get/width/1400/https://haostats.shinyapps.io/splcr-web/)

---

## Installation (R Package)

```r
# install.packages("remotes")
remotes::install_github("haohaostats/SPL-CR")
library(splcr)
```

---

## Quick Start (PBC)

```r
library(splcr)
library(survival)

data(pbc, package = "survival")

pbc_data <- subset(pbc, !is.na(trt))
pbc_data$status <- ifelse(pbc_data$status == 2, 1L,
                          ifelse(pbc_data$status == 1, 2L, 0L))
pbc_data$time <- as.numeric(pbc_data$time)
pbc_data$logbili <- log(pbc_data$bili)
pbc_data$sex <- factor(pbc_data$sex)
pbc_data$stage <- factor(pbc_data$stage)
pbc_data$ascites <- factor(pbc_data$ascites)

pbc_data <- pbc_data[complete.cases(
  pbc_data[, c("time","status","logbili","albumin","protime","ascites","age","sex","stage")]
), ]

res_csh <- spl_csh(
  ~ logbili + albumin + protime + ascites,
  data   = pbc_data,
  folds  = 10,
  k_nn   = 5,
  kernel = "epanechnikov",
  seed   = 20250101
)

res_fg <- spl_fg(
  ~ logbili + albumin + protime + ascites,
  data   = pbc_data,
  folds  = 10,
  k_nn   = 5,
  kernel = "epanechnikov",
  seed   = 20250101
)

print(res_csh)
print(res_fg)
```

---

## Results Snapshot

### Single-model validation output

```text
Smoothed Predictive Likelihood (SPL-CSH)
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

### Multi-model comparison (5 models x 4 kernels)

```r
models <- list(
  "M1: Baseline Clinical"    = ~ age + sex + stage,
  "M2: Core Biochemical"     = ~ logbili + albumin,
  "M3: Extended Biochemical" = ~ logbili + albumin + protime + ascites,
  "M4: Clinical + Core Bio"  = ~ age + sex + stage + logbili + albumin,
  "M5: Full"                 = ~ age + sex + stage + logbili + albumin + protime + ascites
)

tab_csh <- spl_compare(
  models  = models,
  data    = pbc_data,
  method  = "csh",
  kernels = c("uniform", "triangular", "epanechnikov", "quartic"),
  folds   = 10,
  k_nn    = 5,
  seed    = 20250101
)

tab_fg <- spl_compare(
  models  = models,
  data    = pbc_data,
  method  = "fg",
  kernels = c("uniform", "triangular", "epanechnikov", "quartic"),
  folds   = 10,
  k_nn    = 5,
  seed    = 20250101
)

print(spl_compare_wide(tab_csh))
print(spl_compare_wide(tab_fg))
```

```text
SPL-CSH model comparison (metric = Lhat)
Model                        Uniform   Triangular   Epanechnikov     Quartic
M1: Baseline Clinical        -4.6304      -4.6792        -4.6735     -4.6807
M2: Core Biochemical      *  -2.6720   *  -2.7210     *  -2.7155  *  -2.7222
M3: Extended Biochemical     -3.9374      -3.9865        -3.9810     -3.9876
M4: Clinical + Core Bio      -3.1065      -3.1565        -3.1511     -3.1579
M5: Full                     -4.0118      -4.0620        -4.0568     -4.0634

SPL-FG model comparison (metric = Lhat)
Model                        Uniform   Triangular   Epanechnikov     Quartic
M1: Baseline Clinical        -4.1256      -4.1721        -4.1677     -4.1732
M2: Core Biochemical         -4.3483      -4.4351        -4.4223     -4.4452
M3: Extended Biochemical  *  -3.6968   *  -3.7459     *  -3.7391  *  -3.7545
M4: Clinical + Core Bio      -3.8935      -3.9405        -3.9358     -3.9421
M5: Full                     -4.1578      -4.2067        -4.1996     -4.2121

* marks the highest value in each kernel column (larger is better).
```

---

## Interpretation Notes

### Time unit comparability

Both `spl_csh()` and `spl_fg()` are time-unit agnostic.
You may use days, months, or years.

For valid model comparison:

- keep a consistent time unit across all models;
- note that absolute log-score values depend on the chosen unit;
- manuscript real-data examples use **time in days**.

### Fine-Gray implementation detail

`spl_fg()` uses a common all-cause event-time grid within each training fold when extracting baseline CIF jumps, matching the manuscript reference implementation.

---

## Repository Structure

```text
.
|-- DESCRIPTION
|-- NAMESPACE
|-- R/
|   |-- spl_csh.R
|   |-- spl_fg.R
|   |-- spl_compare.R
|   |-- spl_compare_wide.R
|   `-- ...
|-- src/
|   `-- khs.cpp
|-- man/
|-- vignettes/
`-- README.md
```

---

## Dependencies

Core computational dependencies:

- `survival`
- `riskRegression`
- `prodlim`
- `Rcpp`
- `RcppArmadillo`

Interface and deployment dependencies:

- `shiny`
- `rsconnect`

---

## Contact

- **Hao Chen** - University of Sydney
- GitHub: [https://github.com/haohaostats](https://github.com/haohaostats)

---

<p align="center">
  <i>If this repository contributes to your work, a star helps visibility and reproducibility.</i>
</p>



