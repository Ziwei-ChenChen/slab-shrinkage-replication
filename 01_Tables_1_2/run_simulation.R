# Tables 1 and 2: SRR against OLS, RR, EN, LASSO, POET_99% and POET_f under the
# latent space features design (n = 1000, sigma = 1 and 5, N = 250).
#
# Run from this folder:  Rscript run_simulation.R
# One CSV per setting is written to results/ (one row per replication), and
# make_tables.R turns them into Tables 1 and 2.

source("../R/All_eight_estimators.R")
source("../R/benchmarks.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")

N <- 250
n <- 1000
sigmas <- c(1, 5)
p_n_ratios <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95)
f_p_ratios <- c(0.25, 0.50, 0.75)

dir.create("results", showWarnings = FALSE)
start_workers()

one_rep <- function(k, p, f, sigma) {
  seed_rep(123 + k, kind = "Mersenne-Twister")
  X <- draw_latent_x(n, f, p)
  Xt <- cbind(1, X)
  beta <- beta_alt(p + 1)
  y <- draw_y(Xt, beta, sigma)
  foldid <- make_folds(n)

  pen <- penalized_fits(X, y, foldid)
  # POET of Fan, Liao and Mincheva (2013) with soft thresholding, see
  # poet_sigma() in R/benchmarks.R
  poet99 <- poet_est(X, y, share = 0.99)
  poetf <- poet_est(X, y, K = f)
  est <- list(
    OLS = lm.fit(Xt, y)$coefficients,
    RR = pen$RR$beta,
    EN = pen$EN$beta,
    LASSO = pen$LASSO$beta,
    SRR = SRR_est(gram_parts(Xt, y)),
    POET_99 = poet99$beta,
    POET_f = poetf$beta
  )
  L2 <- sapply(est, l2, b = beta)
  data.frame(
    replication = k, seed = 123 + k, n = n, p = p, f = f, sigma = sigma,
    t(setNames(L2, paste0("L2_", names(L2)))),
    K_POET_99 = poet99$K, C_POET_99 = poet99$C, C_POET_f = poetf$C
  )
}

for (sigma in sigmas) {
  for (ratio in p_n_ratios) {
    p <- max(1, round(n * ratio))
    for (f in ceiling(p * f_p_ratios)) {
      file <- sprintf("results/latent_sigma_%g_n_%d_p_%d_f_%d.csv", sigma, n, p, f)
      run_case(file, N, function(k) one_rep(k, p, f, sigma))
    }
  }
}

plan(sequential)
