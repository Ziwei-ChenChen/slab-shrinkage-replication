# Tables B1, B2 and B3: Gaussian Toeplitz covariates, N = 250 replications.
#   B1: Gaussian errors, sigma = 1
#   B2: Gaussian errors, sigma = 5
#   B3: t errors with v = 50/24 degrees of freedom
# n = 500, 1000 and 2500 with p/n = 5%, 10%, 25%, 50%, 75% and 95%.
#
# Run from this folder:  Rscript run_simulation.R            (all three tables)
#                        Rscript run_simulation.R B2         (one table)
#                        Rscript run_simulation.R B2 2500    (one table, n = 2500 only)
#
# A replication with n = 2500 and p = 2375 takes about a minute on one core and
# 0.8 GB of memory, so set memory_gb below to the memory the run may use. The
# number of workers for each setting follows from it.

source("../R/All_eight_estimators.R")
source("../R/benchmarks.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")

designs <- list(
  B1 = list(sigma = 1, df = NULL, tag = "normal_sigma_1"),
  B2 = list(sigma = 5, df = NULL, tag = "normal_sigma_5"),
  B3 = list(sigma = 1, df = 50 / 24, tag = "t_df_50over24")
)
N <- 250
n_vals <- c(500, 1000, 2500)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
memory_gb <- 10

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% names(designs))) designs <- designs[intersect(args, names(designs))]
if (any(args %in% n_vals)) n_vals <- n_vals[n_vals %in% args]
grid <- expand.grid(ratio = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95), n = n_vals)

dir.create("results", showWarnings = FALSE)
start_workers()

one_rep <- function(k, n, p, A, beta, design) {
  seed_rep(123 + k)
  X <- draw_gaussian_x(n, A)
  Xt <- cbind(1, X)
  y <- draw_y(Xt, beta, design$sigma, design$df)
  fit <- ols_parts(Xt, y)
  foldid <- make_folds(n)
  pen <- penalized_fits(X, y, foldid)
  est <- list(
    OLS = fit$beta,
    St = St_est(fit),
    DSh = DSh_est(fit),
    Sh = Sh_est(fit),
    SR = SR_est(fit),
    GSR = GSR_est(fit),
    RR = pen$RR$beta,
    EN = pen$EN$beta,
    LASSO = pen$LASSO$beta
  )
  L2 <- sapply(est, l2, b = beta)
  data.frame(
    replication = k, seed = 123 + k, n = n, p = p,
    t(setNames(L2, paste0("L2_", names(L2)))),
    lambda_RR = pen$RR$lambda, lambda_EN = pen$EN$lambda, lambda_LASSO = pen$LASSO$lambda,
    winner = names(which.min(L2))
  )
}

for (design in designs) {
  for (i in seq_len(nrow(grid))) {
    n <- grid$n[i]
    p <- round(n * grid$ratio[i])
    beta <- beta_alt(p + 1)
    # at its peak a replication holds about nine n x p and p x p matrices
    fit_workers(0.25 + 9 * 8 * (n * p + p^2) / 1e9, memory_gb)
    for (rho in rho_vals) {
      file <- sprintf("results/%s_rho_%g_n_%d_p_%d.csv", design$tag, rho, n, p)
      # the covariance root is slow for large p, so it is skipped for finished settings
      A <- if (case_done(file, N)) NULL else mvrnorm_root(toeplitz_cor(rho, p))
      run_case(file, N, function(k) one_rep(k, n, p, A, beta, design))
    }
  }
}

plan(sequential)
