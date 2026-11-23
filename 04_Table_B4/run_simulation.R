# Table B4: multivariate binomial covariates through a Gaussian copula,
# X_j ~ Binomial(2, q0), with q0 going from 0.01 to 0.25 over the replications
# and beta = (0, U(0.01, 0.30), ...). Gaussian errors with sigma = 1 and 5,
# n = 1000, 2500, 5000 and p = 1, 2, 5; N = 250 replications.
# glmnet needs at least two covariates, so RR, EN and LASSO are NA for p = 1.
#
# Run from this folder:  Rscript run_simulation.R

source("../R/All_eight_estimators.R")
source("../R/benchmarks.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")

N <- 250
sigmas <- c(1, 5)
n_vals <- c(1000, 2500, 5000)
p_vals <- c(1, 2, 5)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
q_range <- c(0.01, 0.25)

dir.create("results", showWarnings = FALSE)
start_workers()

one_rep <- function(k, n, p, A, sigma) {
  seed_rep(123 + k)
  q0 <- q_range[1] + (k - 1) * diff(q_range) / (N - 1)
  X <- draw_binomial_x(n, A, q0)
  beta <- c(0, runif(p, 0.01, 0.30))
  Xt <- cbind(1, X)
  y <- draw_y(Xt, beta, sigma)
  fit <- ols_parts(Xt, y)
  est <- list(
    OLS = fit$beta,
    St = St_est(fit),
    DSh = DSh_est(fit),
    Sh = Sh_est(fit),
    SR = SR_est(fit),
    GSR = GSR_est(fit),
    RR = NA, EN = NA, LASSO = NA
  )
  if (p >= 2) {
    pen <- penalized_fits(X, y, make_folds(n))
    est[c("RR", "EN", "LASSO")] <- lapply(pen, `[[`, "beta")
  }
  L2 <- sapply(est, function(b) if (anyNA(b)) NA else l2(b, beta))
  data.frame(
    replication = k, seed = 123 + k, n = n, p = p, sigma = sigma, q0 = q0,
    t(setNames(L2, paste0("L2_", names(L2)))),
    winner = names(which.min(replace(L2, is.na(L2), Inf)))
  )
}

for (sigma in sigmas) {
  for (n in n_vals) {
    for (p in p_vals) {
      # rho plays no role when p = 1
      for (rho in if (p == 1) NA else rho_vals) {
        A <- mvrnorm_root(toeplitz_cor(if (p == 1) 0 else rho, p))
        file <- sprintf("results/binomial_sigma_%g_rho_%s_n_%d_p_%d.csv", sigma, rho, n, p)
        run_case(file, N, function(k) one_rep(k, n, p, A, sigma))
      }
    }
  }
}

plan(sequential)
