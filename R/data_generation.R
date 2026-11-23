# data_generation.R
#
# Covariate designs, regression parameters and responses of Section 2.2 and
# Appendix B.

# beta = (1, -1, 2, -2, 3, -3, ...) of length m, eq. (B.7)
beta_alt <- function(m) {
  rep(c(1, -1), length.out = m) * rep(1:ceiling(m / 2), each = 2)[1:m]
}

# Toeplitz correlation matrix rho^|i - j|, eq. (B.1)
toeplitz_cor <- function(rho, p) {
  outer(1:p, 1:p, function(i, j) rho^abs(i - j))
}

# MASS::mvrnorm(n, 0, Sigma) computes t(A %*% t(Z)) with A = V diag(sqrt(ev))
# from eigen(Sigma). Computing A once per setting gives the same draws as
# mvrnorm, without repeating the eigen decomposition in every replication.
mvrnorm_root <- function(Sigma) {
  eS <- eigen(Sigma, symmetric = TRUE)
  eS$vectors %*% diag(sqrt(pmax(eS$values, 0)), ncol(Sigma))
}

draw_gaussian_x <- function(n, A) {
  t(A %*% t(matrix(rnorm(n * ncol(A)), n)))
}

# Binomial(2, q) covariates through a Gaussian copula, eq. (B.2)
draw_binomial_x <- function(n, A, q) {
  Z <- draw_gaussian_x(n, A)
  matrix(qbinom(pnorm(Z), size = 2, prob = q), nrow = n)
}

# Latent space features (Tables 1 and 2): X = A Z + noise, A is n x f, Z is f x p
draw_latent_x <- function(n, f, p) {
  A <- matrix(rnorm(n * f), nrow = n)
  Z <- matrix(rnorm(f * p), nrow = f)
  A %*% Z + matrix(rnorm(n * p, sd = sqrt(1e-6)), nrow = n)
}

# Response with Gaussian errors, or t errors when df is given, eqs. (B.5)-(B.6)
draw_y <- function(Xtilde, beta, sigma = 1, df = NULL) {
  mean_y <- drop(Xtilde %*% beta)
  if (is.null(df)) {
    rnorm(length(mean_y), mean = mean_y, sd = sigma)
  } else {
    mean_y + rt(length(mean_y), df = df)
  }
}

# Replication k uses set.seed(123 + k), with L'Ecuyer-CMRG for Tables B1-B4
# and B6 and Mersenne-Twister for the other tables and figures. The generator
# is named explicitly so that the results do not depend on how the work is
# split over cores.
seed_rep <- function(seed, kind = "L'Ecuyer-CMRG") {
  set.seed(seed, kind = kind, normal.kind = "Inversion", sample.kind = "Rejection")
}

l2 <- function(a, b) sqrt(sum((a - b)^2))
