# Replication code for "Slab and Shrinkage Linear Regression Estimation"

This repository contains the R code and the stored results for the simulation
and real-data tables and figures of the paper and its Supplementary Material.
Each numbered folder computes one group of outputs.

| Folder | Output | Scripts |
|---|---|---|
| `01_Tables_1_2` | Tables 1 and 2 (latent-factor design) | `run_simulation.R`, `make_tables.R` |
| `03_Tables_B1_B3` | Tables B1, B2, and B3 (Gaussian Toeplitz covariates) | `run_simulation.R`, `make_tables.R` |
| `04_Table_B4` | Table B4 (binomial covariates) | `run_simulation.R`, `make_tables.R` |
| `05_Table_B5` | Table B5 (computing time) | `run_simulation.R`, `make_tables.R` |
| `06_Table_B6` | Table B6 (oracle RR and SR) | `run_simulation.R`, `make_tables.R` |
| `07_Tables_D7_D9` | Tables D7, D8, and D9 and Figures D1 and D2 (cyber-sickness data) | `run_real_data.R`, `make_tables.R`, `make_figures.R` |
| `08_Figures_SM2` | Figures SM.2.1 to SM.2.9 | `run_simulation.R`, `make_figures.R` |

Table 3 lists the estimators and has no numerical content.


## Requirements

The results were computed with R 4.4.2 and R 4.5.2 and the packages glmnet,
future, future.apply, and parallelly. The figures also need ggplot2 and
corrplot, and Figures SM.2.1 and SM.2.2 use MASS, which comes with R.

Table B5 times the savvySh package from CRAN. The timings in the paper were
taken with version 0.1.3, and `05_Table_B5/run_simulation.R` installs savvySh
from CRAN when it is missing or older than that.


## Running the code

Every script is run from its own folder, for example

```
cd 03_Tables_B1_B3
Rscript run_simulation.R
Rscript make_tables.R
```

`run_simulation.R` (`run_real_data.R` in `07_Tables_D7_D9`) writes one CSV
file per setting to `results/`, with one row per replication.
`make_tables.R` reads `results/` and writes each table to `tables/` as a
LaTeX file, together with a CSV file of the numbers it shows.
`make_figures.R` writes every figure to `figures/` as EPS and PNG, under the
file names used in the paper.

The `results/` folders hold the results reported in the paper, so
`make_tables.R` and `make_figures.R` can be run straight away and give the
tables and figures of the paper. `run_simulation.R` skips every setting that
is already complete in `results/`. To compute a setting again, delete its CSV
file (or the whole `results/` folder) first. A run that stops can be started
again with the same command, and only the missing replications are computed.

The replications run in parallel on at most 10 workers. This limit is
`max_workers` at the top of `R/simulation_tools.R`. Each worker is a separate
R process, so with a multithreaded BLAS (OpenBLAS or MKL) set
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `OMP_NUM_THREADS=1`.

Approximate computing cost on one core:

| Folder | Cost |
|---|---|
| `01_Tables_1_2` | 5 core-hours |
| `03_Tables_B1_B3` | 250 core-hours, of which about 210 are for n = 2,500 with p/n of 50% or more |
| `04_Table_B4` | 1 core-hour |
| `05_Table_B5` | 10 hours on a single core |
| `06_Table_B6` | 1 core-hour |
| `07_Tables_D7_D9` | 12 core-hours |
| `08_Figures_SM2` | up to an hour and 9 GB of memory per replication for n = 10,000 |

In `03_Tables_B1_B3`, `Rscript run_simulation.R B2` runs one table and
`Rscript run_simulation.R B2 2500` runs it for n = 2,500 only. A replication
with n = 2,500 and p = 2,375 needs about 0.8 GB of memory. Set `memory_gb` at
the top of the script to the memory the run may use, and each setting then
runs on as many workers as fit in it. `08_Figures_SM2/run_simulation.R` has
four parts (`SM2_1_2`, `SM2_3_4`, `SM2_5_8`, and `SM2_9`), and
`Rscript run_simulation.R SM2_5_8` runs one of them. Its `memory_gb` works in
the same way. Results computed on several machines can be combined by copying
the CSV files into one `results/` folder.

Table B5 runs its replications one after another on a single core, since
they are timings, and should be run on an otherwise idle machine. The time
ratios of a new run depend on the machine.


## Shared code

The folder `R/` holds the code used by all tables.

- `All_eight_estimators.R` has OLS, St, DSh, Sh, SR, GSR, LSh, and SRR, with
  Sigma = X'X as in the paper. The comments give the theorem and equation
  numbers of the paper.
- `benchmarks.R` has RR, EN, and LASSO through `cv.glmnet` with 10-fold
  cross-validation, and POET (Fan, Liao, and Mincheva, 2013) with K chosen to
  explain 99% of the variance or given.
- `data_generation.R` has the covariate designs, the regression coefficients,
  the responses, and the seeds.
- `simulation_tools.R` runs one setting in parallel batches and saves the
  results after every batch.
- `latex_tables.R` has the helpers used by every `make_tables.R`.

`07_Tables_D7_D9/glm_irls.R` fits the logistic and Poisson GLMs by IRLS, with
the weighted least squares step done by OLS, St, DSh, Sh, SR, or GSR.


## Random numbers

Each replication sets its own seed, so its results do not depend on the
number of workers or on restarts.

| Output | Generator | Seed of replication k |
|---|---|---|
| Tables 1 and 2 | Mersenne-Twister | 123 + k |
| Tables B1-B4 and B6 | L'Ecuyer-CMRG | 123 + k |
| Table B5 | Mersenne-Twister | 123 + k for the data, 100123 + k for the folds |
| Tables D7-D9 | Mersenne-Twister | 123 + k for the 70/30 split, 10000 + k for the folds |
| Figures SM.2.1 and SM.2.2 | Mersenne-Twister | one stream per experiment, started with 123 |
| Figures SM.2.3 and SM.2.4 | Mersenne-Twister | s + k, with s = 123, 1123, 2123, and 3123 for (beta, lambda_max) = (alternating, 1), (alternating, 10), (ones, 1), and (ones, 10) |
| Figures SM.2.5-SM.2.9 | Mersenne-Twister | 123 + k |

Gaussian covariates with a non-diagonal covariance matrix (Tables B1-B4 and
B6 with rho different from 0, and Figures SM.2.1 and SM.2.2) are drawn as in
`MASS::mvrnorm`, through the eigendecomposition of the covariance matrix. The
signs of the eigenvectors depend on the LAPACK library that R uses. On another
installation these replications can differ from the stored ones while
following the same distribution, and a full rerun then agrees with the paper
up to Monte Carlo error. All other stored results are reproduced to rounding
error, apart from the timings of Table B5.


## Cyber-sickness data

The data in `07_Tables_D7_D9/data/` were built from the physiological
recordings of the CyberSicknessClassification repository
(https://github.com/shovonis/CyberSicknessClassification), as described in
Appendix D. `data_for_regression1.csv` (25,870 rows) has the 13 features at
lag t-1, and `data_for_regression10.csv` (25,663 rows) has the same features
at lags t-10 to t-1 (130 features). In both files the first column is the
intercept and the last column is the FMS score at time t. The 67 scores of
2.5 and the 74 scores of 3.5 are rounded down to 2 and 3 when the data are
read, as in Appendix D.

With the raw features and M = 10, the weighted least squares step of DSh
fails in three splits of the FMS 0 vs 2 and 0 vs 3 comparisons, because
almost every fitted probability reaches 0 or 1. IRLS then stops at the last
coefficients with a warning, and the split is kept.


## Figures

| Figure | Files in `figures/` |
|---|---|
| D1 | `initialclass_distribution`, `4class_distribution` |
| D2 | `heatmap13features`, `heatmap39of130features` |
| SM.2.1 | `plot_1_p=5_rho=*`, `plot_1_p=10_rho=*` |
| SM.2.2 | `plot_1_p=*_Independent`, `plot_2_p=*_Dependent` |
| SM.2.3 | `Boxplot_Ratio_Beta_theta_alt_*`, `Boxplot_Ratio_Lambda_theta_alt_*` |
| SM.2.4 | `Boxplot_Ratio_Beta_theta_all_*`, `Boxplot_Ratio_Lambda_theta_all_*` |
| SM.2.5-SM.2.8 | `Faceted_OLS_Boxplot_Ratio_Beta_n1000_n10000` and the same for St, DSh, and GSR |
| SM.2.9 | `Faceted_Boxplot_Ratio_St_vs_OLS` and the same for DSh and GSR |

`08_Figures_SM2/make_figures.R` also writes the median, mean, and standard
deviation of every box in Figures SM.2.5-SM.2.9 to
`summaries/SM2_5_9_ratio_summaries.csv`.
