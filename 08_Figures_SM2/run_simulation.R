# Figures SM.2.1-SM.2.9 of the Supplementary Material. They show the
# eigenvalues of the sample covariance matrix, and the L1 errors of OLS, St,
# DSh and GSR when the covariates are independent or dependent with the same
# eigenvalues.
#
#   SM2_1_2  SM.2.1  sample eigenvalues, Toeplitz covariance, p = 5, 10
#            SM.2.2  sample eigenvalues, covariance Q diag(lambda) Q' with
#                    Q = I or a random orthogonal matrix
#   SM2_3_4  SM.2.3-SM.2.4  OLS L1 error ratio (independent / dependent
#            covariates), n = 250
#   SM2_5_8  SM.2.5-SM.2.8  the same ratio for OLS, St, DSh and GSR,
#            n = 1,000 to 10,000
#   SM2_9    SM.2.9  L1 error of St, DSh and GSR over that of OLS
#
# Every setting has N = 250 replications, and sigma = 5 in SM.2.3-SM.2.9.
#
# Run from this folder:  Rscript run_simulation.R           (all four parts)
#                        Rscript run_simulation.R SM2_5_8   (one part)
#
# The settings with n = 10,000 are large. One replication with p = 9,500
# needs about 9 GB of memory and an hour on one core, so set memory_gb below
# to the memory of the machine. The number of workers follows from it.

source("../R/All_eight_estimators.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")

N <- 250
sigma <- 5
n_small <- 250
n_large <- c(1000, 2500, 5000, 10000)
p_n_large <- c(0.25, 0.5, 0.75, 0.9, 0.95)
methods <- c("OLS", "St", "DSh", "GSR")
memory_gb <- 64

args <- commandArgs(trailingOnly = TRUE)
parts <- if (length(args) > 0) args else c("SM2_1_2", "SM2_3_4", "SM2_5_8", "SM2_9")
for (part in parts) dir.create(file.path("results", part), recursive = TRUE, showWarnings = FALSE)


# Shared pieces ---------------------------------------------------------------

l1 <- function(a, b) sum(abs(a - b))

sample_eigen <- function(X) eigen(cov(X), symmetric = TRUE, only.values = TRUE)$values

# Q factor of B B' for a random symmetric B = (A + A') / 2, A ~ N(0, 1)
random_rotation <- function(p) {
  A <- matrix(rnorm(p * p), p)
  A <- 0.5 * (A + t(A))
  qr.Q(qr(tcrossprod(A)))
}

estimate <- function(method, fit) {
  switch(method, OLS = OLS_est(fit), St = St_est(fit), DSh = DSh_est(fit), GSR = GSR_est(fit))
}

l1_errors <- function(X, y, beta, methods) {
  fit <- ols_parts(X, y, eigen = "GSR" %in% methods)
  sapply(methods, function(m) l1(estimate(m, fit), beta))
}

# One replication of SM.2.3-SM.2.8. The population eigenvalues lambda are
# U(0, lambda_max) and shared by both designs: X1 = Z1 diag(sqrt(lambda))
# (independent) and X2 = Z2 diag(sqrt(lambda)) Q', so Cov(X2) = Q diag(lambda) Q'.
# Replication k uses the seed seed + k.
indep_dep_rep <- function(k, n, p, lambda_max, beta, methods, seed = 123) {
  seed_rep(seed + k, kind = "Mersenne-Twister")
  lambda <- runif(p, 0, lambda_max)
  Q <- random_rotation(p)
  errors <- list()
  eigen_l1 <- c()
  for (design in c("independent", "dependent")) {
    X <- matrix(rnorm(n * p), n) * rep(sqrt(lambda), each = n)
    if (design == "dependent") X <- tcrossprod(X, Q)
    eigen_l1[design] <- l1(sample_eigen(X), sort(lambda, decreasing = TRUE))
    X <- cbind(1, X)
    errors[[design]] <- l1_errors(X, draw_y(X, beta, sigma), beta, methods)
  }
  ratio <- errors$independent / errors$dependent
  data.frame(replication = k, ratio_lambda = eigen_l1[[1]] / eigen_l1[[2]],
             t(setNames(ratio, paste0("ratio_beta_", methods))))
}

# One replication of SM.2.9: independent N(0, 1) covariates and beta = 1
vs_ols_rep <- function(k, n, p) {
  seed_rep(123 + k, kind = "Mersenne-Twister")
  X <- cbind(1, matrix(rnorm(n * p), n))
  beta <- rep(1, p + 1)
  e <- l1_errors(X, draw_y(X, beta, sigma), beta, methods)
  data.frame(replication = k, t(setNames(e[-1] / e[1], paste0("ratio_", methods[-1], "_OLS"))),
             t(setNames(e, paste0("L1_", methods))))
}

# At its peak a replication holds about six n x p and p x p matrices, so the
# large settings run on fewer workers
use_workers <- function(n, p) {
  fit_workers(0.2 + 6 * 8 * (n * p + p^2) / 1e9, memory_gb)
}


# SM.2.1 and SM.2.2 -----------------------------------------------------------
# These are small, and each experiment runs as one stream of random numbers
# started with set.seed(123).

half_eigen <- function(p) rep(c(1, 0.5), c(ceiling(p / 2), p - ceiling(p / 2)))

eigen_experiment <- function(file, p, make_cov) {
  E <- t(sapply(1:N, function(k) {
    Sigma <- make_cov(p)
    sample_eigen(MASS::mvrnorm(n_small, rep(0, p), Sigma))
  }))
  colnames(E) <- paste0("E", 1:p)
  write.csv(data.frame(replication = 1:N, E), file, row.names = FALSE)
}

# Each experiment is one stream, so the part is run whole or skipped whole
sm2_1_2_files <- c(
  sprintf("results/SM2_1_2/toeplitz_p%d_rho%g.csv", rep(c(5, 10), each = 3), c(-0.5, 0, 0.5)),
  sprintf("results/SM2_1_2/%s_p%d.csv", rep(c("independent", "dependent"), each = 2), c(5, 10))
)
if ("SM2_1_2" %in% parts && all(sapply(sm2_1_2_files, case_done, N = N))) {
  cat("done already: SM2_1_2\n")
} else if ("SM2_1_2" %in% parts) {
  seed_rep(123, kind = "Mersenne-Twister")
  for (p in c(5, 10)) for (rho in c(-0.5, 0, 0.5)) {
    eigen_experiment(sprintf("results/SM2_1_2/toeplitz_p%d_rho%g.csv", p, rho), p,
                     function(p) toeplitz_cor(rho, p))
  }
  seed_rep(123, kind = "Mersenne-Twister")
  for (p in c(5, 10)) {
    eigen_experiment(sprintf("results/SM2_1_2/independent_p%d.csv", p), p,
                     function(p) diag(half_eigen(p)))
  }
  seed_rep(123, kind = "Mersenne-Twister")
  for (p in c(5, 10)) {
    eigen_experiment(sprintf("results/SM2_1_2/dependent_p%d.csv", p), p, function(p) {
      Q <- random_rotation(p)
      Q %*% diag(half_eigen(p)) %*% t(Q)
    })
  }
}


# SM.2.3 to SM.2.9 ------------------------------------------------------------

# Each (beta, lambda_max) setting has its own seeds, 123, 1123, 2123 and 3123
# plus k.
if ("SM2_3_4" %in% parts) {
  start_workers()
  seed <- 123
  for (beta_design in c("alt", "ones")) for (lambda_max in c(1, 10)) {
    for (p in ceiling(n_small * c(0.05, 0.1, 0.25, 0.5, 0.75))) {
      beta <- if (beta_design == "alt") rep(c(1, -1), length.out = p + 1) else rep(1, p + 1)
      file <- sprintf("results/SM2_3_4/beta_%s_lambdamax%g_p%d.csv", beta_design, lambda_max, p)
      run_case(file, N, function(k) indep_dep_rep(k, n_small, p, lambda_max, beta, "OLS", seed))
    }
    seed <- seed + 1000
  }
}

if ("SM2_5_8" %in% parts) {
  for (n in n_large) for (p in ceiling(n * p_n_large)) {
    w <- use_workers(n, p)
    file <- sprintf("results/SM2_5_8/n%d_p%d.csv", n, p)
    run_case(file, N, function(k) indep_dep_rep(k, n, p, 1, rep(1, p + 1), methods), batch_size = max(w, 10))
  }
}

if ("SM2_9" %in% parts) {
  for (n in n_large) for (p in ceiling(n * p_n_large)) {
    w <- use_workers(n, p)
    file <- sprintf("results/SM2_9/n%d_p%d.csv", n, p)
    run_case(file, N, function(k) vs_ols_rep(k, n, p), batch_size = max(w, 10))
  }
}

plan(sequential)
