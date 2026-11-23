# Table B5: computational time relative to OLS (n = 500, Gaussian Toeplitz
# covariates, sigma = 1 and 5, N = 250 replications).
#
# The table shows the computing time of the savvySh package, so St, DSh, Sh,
# SR, GSR and SRR are timed as calls to savvySh (version 0.1.3 or later, from
# CRAN).
#
# The timings are taken one replication after the other on a single core, so
# run this on an otherwise idle machine:  Rscript run_simulation.R
#
# What is timed, as in the notes of Table B5:
#   OLS          lm.fit
#   RR/EN/LASSO  one cv.glmnet call each (10 folds, shared fold labels)
#   St/DSh/Sh    savvySh(X, y, model_class = "Multiplicative"), which returns
#                St, DSh and Sh (include_Sh = TRUE by default)
#   SR/GSR       savvySh(X, y, model_class = "Slab")
#   SRR          savvySh(X, y, model_class = "ShrinkageRR")
#   POET         POET_99% with soft thresholding

# savvySh is installed from CRAN when it is missing or older than 0.1.3, before
# any timing starts
if (!nzchar(system.file(package = "savvySh")) || packageVersion("savvySh") < "0.1.3") {
  install.packages("savvySh", repos = "https://cloud.r-project.org")
}
library(savvySh)
if (package_version(getNamespaceVersion("savvySh")) < "0.1.3") stop("Table B5 needs savvySh 0.1.3 or later")
cat("savvySh", getNamespaceVersion("savvySh"), "\n")
source("../R/benchmarks.R")
source("../R/data_generation.R")

N <- 250
n <- 500
sigmas <- c(1, 5)
p_n_ratios <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)

dir.create("results", showWarnings = FALSE)

# The clock is read with Sys.time(). On Windows proc.time() moves in steps of
# 10-15 ms, longer than an OLS fit with n = 500.
time_call <- function(expr) {
  gc(FALSE)
  start <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), start, units = "secs"))
}

glmnet_fit <- function(X, y, alpha, foldid) {
  cv <- cv.glmnet(X, y, alpha = alpha, foldid = foldid, type.measure = "mse", intercept = TRUE)
  as.numeric(coef(cv, s = "lambda.min"))
}

one_rep <- function(k, p, A, beta, sigma) {
  # separate seeds for the data and the folds
  seed_rep(123 + k, kind = "Mersenne-Twister")
  X <- draw_gaussian_x(n, A)
  y <- draw_y(cbind(1, X), beta, sigma)
  seed_rep(123 + 100000 + k, kind = "Mersenne-Twister")
  foldid <- make_folds(n)

  T_OLS <- time_call(lm.fit(cbind(1, X), y)$coefficients)
  T_RR <- time_call(glmnet_fit(X, y, 0, foldid))
  T_EN <- time_call(glmnet_fit(X, y, 0.5, foldid))
  T_LASSO <- time_call(glmnet_fit(X, y, 1, foldid))
  T_St_DSh_Sh <- time_call(savvySh(X, y, model_class = "Multiplicative"))
  T_SR_GSR <- time_call(savvySh(X, y, model_class = "Slab"))
  T_SRR <- time_call(savvySh(X, y, model_class = "ShrinkageRR"))
  T_POET <- time_call(poet_est(X, y, share = 0.99))

  data.frame(
    replication = k, n = n, p = p, sigma = sigma,
    T_OLS, T_RR, T_EN, T_LASSO, T_POET, T_St_DSh_Sh, T_SR_GSR, T_SRR
  )
}

# The first call of each procedure in a new R session carries loading costs,
# so the session first runs one untimed replication (k = 0, not saved).
# one_rep() sets its own seed, so this does not change any saved replication.
warmed_up <- FALSE

for (sigma in sigmas) {
  for (ratio in p_n_ratios) {
    p <- round(n * ratio)
    beta <- beta_alt(p + 1)
    for (rho in rho_vals) {
      file <- sprintf("results/runtime_sigma_%g_rho_%g_n_%d_p_%d.csv", sigma, rho, n, p)
      done <- if (file.exists(file)) read.csv(file) else NULL
      todo <- setdiff(1:N, done$replication)
      if (length(todo) == 0) next
      cat("running", basename(file), "\n")
      A <- mvrnorm_root(toeplitz_cor(rho, p))
      if (!warmed_up) {
        invisible(one_rep(0, p, A, beta, sigma))
        warmed_up <- TRUE
      }
      for (k in todo) {
        done <- rbind(done, one_rep(k, p, A, beta, sigma))
        if (k %% 25 == 0 || k == N) write.csv(done, file, row.names = FALSE)
      }
    }
  }
}
