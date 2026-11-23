# Table B6: oracle comparison of RR and SR (Appendix B.4). The true beta and
# sigma are used to tune each estimator for the realized design.
#   ORR      ridge with the lambda that minimizes its MSE
#   OSR(1)   SR with u = 1 and mu**(1) = sigma^2 / (1'beta)^2 of Theorem 2
#   OSR(U)   SR with u chosen among the eigenvectors of X'X by smallest MSE
# n = 500, Gaussian Toeplitz covariates, sigma = 1 and 5, N = 500 replications.
# When p is odd, beta = (1, -1, 2, -2, ...) has 1'beta = 0, so mu**(1) = Inf
# and OSR(1) restricts OLS to the hyperplane 1'beta = 0.
#
# Run from this folder:  Rscript run_simulation.R

source("../R/All_eight_estimators.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")

N <- 500
n <- 500
sigmas <- c(1, 5)
p_n_ratios <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)

dir.create("results", showWarnings = FALSE)
start_workers()

# MSE of ridge with penalty lambda, written in the eigenbasis of X'X
ridge_mse <- function(lambda, d, beta_eig, sigma2) {
  if (is.infinite(lambda)) return(sum(beta_eig^2))
  sum((lambda^2 * beta_eig^2 + sigma2 * d) / (d + lambda)^2)
}

# lambda = median(d) t / (1 - t) maps t in [0, 1) onto [0, Inf); the end
# points lambda = 0 (OLS) and lambda = Inf (beta = 0) are checked as well
oracle_lambda <- function(d, beta_eig, sigma2) {
  s <- median(d)
  opt <- optimize(function(t) ridge_mse(s * t / (1 - t), d, beta_eig, sigma2), c(0, 1 - 1e-10), tol = 1e-10)
  cand <- c(0, s * opt$minimum / (1 - opt$minimum), Inf)
  mse <- sapply(cand, ridge_mse, d = d, beta_eig = beta_eig, sigma2 = sigma2)
  list(lambda = cand[which.min(mse)], mse = min(mse))
}

one_rep <- function(k, p, A, beta, sigma) {
  seed_rep(123 + k)
  Xt <- cbind(1, draw_gaussian_x(n, A))
  y <- draw_y(Xt, beta, sigma)
  sigma2 <- sigma^2
  fit <- gram_parts(Xt, y)
  d <- fit$lambda
  U <- fit$U
  if (min(d) <= 1e-12 * max(d)) stop("X'X is numerically singular")
  beta_eig <- drop(crossprod(U, beta))
  w <- drop(crossprod(U, fit$Xty))
  ols <- drop(U %*% (w / d))
  ols_eig <- w / d
  trace_inv <- sum(1 / d)

  rr <- oracle_lambda(d, beta_eig, sigma2)
  orr <- if (is.infinite(rr$lambda)) 0 * beta else drop(U %*% (w / (d + rr$lambda)))

  # OSR(1)
  u <- rep(1, length(beta))
  z <- drop(crossprod(U, u))
  a <- SR_a(z, d)
  mu1 <- SR_mu(sigma2, sum(u * beta))
  c1 <- SR_c(mu1, a$a1)
  osr1 <- ols - c1 * sum(u * ols) * drop(U %*% (z / d))
  mse1 <- SR_mse(c1, a, trace_inv, sigma2, sum(u * beta))

  # OSR(U): for u = u_j, a_l = 1 / d_j^l and u'beta = beta_eig[j]
  aU <- list(a1 = 1 / d, a2 = 1 / d^2)
  muU <- SR_mu(sigma2, beta_eig)
  cU <- SR_c(muU, aU$a1)
  mseU <- SR_mse(cU, aU, trace_inv, sigma2, beta_eig)
  j <- which.min(mseU)
  osrU <- ols - cU[j] * ols_eig[j] * U[, j] / d[j]

  data.frame(
    replication = k, seed = 123 + k, n = n, p = p, sigma = sigma,
    SE2_ORR = sum((orr - beta)^2), SE2_OSR1 = sum((osr1 - beta)^2), SE2_OSRU = sum((osrU - beta)^2),
    MSE_ORR = rr$mse, MSE_OSR1 = mse1, MSE_OSRU = mseU[j],
    lambda_ORR = rr$lambda, mu_OSR1 = mu1, mu_OSRU = muU[j], direction_OSRU = j
  )
}

for (sigma in sigmas) {
  for (ratio in p_n_ratios) {
    p <- round(n * ratio)
    beta <- beta_alt(p + 1)
    for (rho in rho_vals) {
      A <- mvrnorm_root(toeplitz_cor(rho, p))
      file <- sprintf("results/oracle_sigma_%g_rho_%g_n_%d_p_%d.csv", sigma, rho, n, p)
      run_case(file, N, function(k) one_rep(k, p, A, beta, sigma))
    }
  }
}

plan(sequential)
