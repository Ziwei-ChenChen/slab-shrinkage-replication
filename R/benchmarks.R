# benchmarks.R
#
# Benchmark estimators: RR, EN and LASSO from glmnet with 10-fold
# cross-validation, and the POET plug-in estimator.

library(glmnet)

# RR (alpha = 0), EN (alpha = 0.5) and LASSO (alpha = 1) on the same folds.
# X has no intercept column, and glmnet fits the intercept itself.
penalized_fits <- function(X, y, foldid, alphas = c(RR = 0, EN = 0.5, LASSO = 1)) {
  lapply(alphas, function(alpha) {
    cv <- cv.glmnet(X, y, alpha = alpha, foldid = foldid, type.measure = "mse", intercept = TRUE)
    list(beta = as.numeric(coef(cv, s = "lambda.min")), lambda = cv$lambda.min)
  })
}

# Fold labels, drawn the same way as inside cv.glmnet
make_folds <- function(n, nfolds = 10) sample(rep(1:nfolds, length.out = n))


# POET -----------------------------------------------------------------------

# POET covariance estimate of Fan, Liao and Mincheva (2013, JRSS-B), their
# eq. (2.7), for the p x n data matrix Y with K factors (n is T in their paper).
# The covariance sigma_ij of the residuals u_it left after the first K
# principal components keeps its diagonal, and each off-diagonal entry is
# soft-thresholded, sign(sigma_ij) (|sigma_ij| - tau_ij)_+, with the adaptive
# threshold of their eq. (3.2),
#   tau_ij = C omega_n sqrt(theta_ij),  omega_n = 1/sqrt(p) + sqrt(log(p)/n),
#   theta_ij = (1/n) sum_t (u_it u_jt - sigma_ij)^2.
# C = 0.5 as in their simulations. If this does not give a positive definite
# matrix, the smallest such constant C_min is found by bisection and
# C_min + 0.1 is used (their Section 4.1). The C used is returned.
#
# K = NULL picks the smallest K whose eigenvalues explain at least `share` of
# the total variance (POET_99%).
poet_sigma <- function(Y, K = NULL, share = 0.99, C = 0.5) {
  p <- nrow(Y)
  n <- ncol(Y)
  Y <- Y - rowMeans(Y)
  if (p <= n) {
    eig <- eigen(tcrossprod(Y), symmetric = TRUE)
    ev <- pmax(eig$values, 0)
  } else {
    eig <- eigen(crossprod(Y), symmetric = TRUE)
    ev <- pmax(eig$values, 0)
  }
  if (is.null(K)) K <- which(cumsum(ev) / sum(ev) >= share)[1]
  K <- max(1, min(K, p))

  if (p <= n) {
    G <- eig$vectors[, 1:K, drop = FALSE]
    uhat <- Y - G %*% crossprod(G, Y)
    Lowrank <- G %*% (ev[1:K] / n * t(G))
  } else {
    V <- eig$vectors[, 1:K, drop = FALSE]
    YV <- Y %*% V
    uhat <- Y - tcrossprod(YV, V)
    Lowrank <- tcrossprod(YV) / n
  }

  Su <- tcrossprod(uhat) / n

  # threshold for C = 1
  omega <- 1 / sqrt(p) + sqrt(log(p) / n)
  tau1 <- omega * sqrt(pmax(tcrossprod(uhat^2) / n - Su^2, 0))
  soft <- function(C) {
    S <- sign(Su) * pmax(abs(Su) - C * tau1, 0)
    diag(S) <- diag(Su)
    S
  }
  is_pd <- function(S) !inherits(try(chol(S), silent = TRUE), "try-error")

  SuT <- soft(C)
  if (!is_pd(SuT)) {
    # Above C_max every off-diagonal entry is set to zero.
    off <- row(Su) != col(Su) & tau1 > 0
    C_max <- max(C, abs(Su[off]) / tau1[off])
    lo <- C
    hi <- C_max
    if (is_pd(soft(hi))) {
      while (hi - lo > 0.01) {
        mid <- (lo + hi) / 2
        if (is_pd(soft(mid))) hi <- mid else lo <- mid
      }
      C <- if (is_pd(soft(hi + 0.1))) hi + 0.1 else hi
    } else {
      C <- C_max
    }
    SuT <- soft(C)
  }
  list(SigmaY = Lowrank + SuT, K = K, C = C)
}

# POET plug-in regression estimate. X has no intercept column, and the output
# is (intercept, slopes) with the K and C that POET used.
poet_est <- function(X, y, K = NULL, share = 0.99) {
  n <- nrow(X)
  xbar <- colMeans(X)
  Xc <- sweep(X, 2, xbar)
  yc <- y - mean(y)
  S <- poet_sigma(t(Xc), K = K, share = share)
  slope <- solve(S$SigmaY, drop(crossprod(Xc, yc)) / n)
  list(beta = c(mean(y) - sum(xbar * slope), slope), K = S$K, C = S$C)
}
