# All_eight_estimators.R
#
# The eight estimators of the paper for the linear model y = X beta + e:
# OLS, St, DSh, Sh, SR, GSR, LSh and SRR.
#
# X is the n x (p+1) design matrix and its first column is the intercept.
# Sigma = X'X as in the paper, and theorem and equation numbers refer to the
# main text.
#
# Typical use:
#   fit <- ols_parts(X, y)
#   St_est(fit); DSh_est(fit); Sh_est(fit); SR_est(fit); GSR_est(fit); SRR_est(fit)
#
# ols_parts() runs the least squares fit and the eigen decomposition of X'X
# once, so the estimators that need them share the work.


# Shared quantities -----------------------------------------------------------

# X'X, X'y, y'y and, if requested, the eigen decomposition X'X = U diag(lambda) U'
gram_parts <- function(X, y, eigen = TRUE) {
  fit <- list(
    XtX = crossprod(X),
    Xty = drop(crossprod(X, y)),
    yty = drop(crossprod(y)),
    n = nrow(X)
  )
  if (eigen) fit <- add_eigen(fit)
  fit
}

add_eigen <- function(fit) {
  if (is.null(fit$U)) {
    eig <- eigen(fit$XtX, symmetric = TRUE)
    fit$lambda <- eig$values
    fit$U <- eig$vectors
  }
  fit
}

# gram_parts() plus the OLS estimate and sigma^2 = RSS / (n - p - 1)
ols_parts <- function(X, y, eigen = TRUE) {
  ls_fit <- lm.fit(X, y)
  fit <- gram_parts(X, y, eigen = eigen)
  fit$beta <- unname(ls_fit$coefficients)
  fit$sigma2 <- sum(ls_fit$residuals^2) / (nrow(X) - ncol(X))
  fit
}

# Diagonal of (X'X)^{-1}
inv_diag <- function(fit) {
  if (is.null(fit$U)) {
    diag(chol2inv(chol(fit$XtX)))
  } else {
    drop(fit$U^2 %*% (1 / fit$lambda))
  }
}


# OLS -------------------------------------------------------------------------

OLS_est <- function(fit) fit$beta


# Multiplicative shrinkage: St, DSh and Sh ------------------------------------

# Stein estimator, a* = beta'beta / (beta'beta + sigma^2 tr(Sigma^{-1}))
St_est <- function(fit) {
  b2 <- sum(fit$beta^2)
  a <- b2 / (b2 + fit$sigma2 * sum(inv_diag(fit)))
  a * fit$beta
}

# Diagonal shrinkage, b_k = beta_k^2 / (beta_k^2 + sigma^2 (Sigma^{-1})_kk)
DSh_est <- function(fit) {
  b2 <- fit$beta^2
  b2 / (b2 + fit$sigma2 * inv_diag(fit)) * fit$beta
}

# Sh, Theorem 7 iii). C* = beta beta' Sigma / (sigma^2 + beta' Sigma beta), so
# the plug-in estimate C* beta is beta q / (sigma^2 + q) with q = beta' Sigma beta
Sh_est <- function(fit) {
  b <- fit$beta
  q <- sum(b * drop(fit$XtX %*% b))
  q / (fit$sigma2 + q) * b
}


# Slab regression: SR and GSR -------------------------------------------------

# a_l(u) = u' Sigma^{-l} u for l = 1, 2, written with z = U'u
SR_a <- function(z, lambda) {
  list(a1 = sum(z^2 / lambda), a2 = sum(z^2 / lambda^2))
}

# mu**(u) of Theorem 2, eq. (6), with beta_u = u'beta. Delta(u) = a2 (u'beta)^2,
# so mu** = sigma^2 / (u'beta)^2, and mu** = Inf when u'beta = 0.
# The arguments may be vectors (one entry per direction), as in Table B6.
SR_mu <- function(sigma2, beta_u) {
  ifelse(beta_u != 0, sigma2 / beta_u^2, Inf)
}

# SR is beta - c Sigma^{-1} u u'beta with c = mu / (1 + mu a1)
SR_c <- function(mu, a1) {
  ifelse(is.finite(mu), mu / (1 + mu * a1), 1 / a1)
}

# MSE of SR in eq. (4) for a given c, with trace_inv = tr(Sigma^{-1}).
# With A = Sigma^{-1} u u', tr(A Sigma^{-1}) = a2, tr(A Sigma^{-1} A') = a1 a2
# and beta'A'A beta = a2 (u'beta)^2.
SR_mse <- function(c, a, trace_inv, sigma2, beta_u) {
  sigma2 * (trace_inv - 2 * c * a$a2 + c^2 * a$a1 * a$a2) + c^2 * a$a2 * beta_u^2
}

# Bona fide SR with slab direction u (u = 1 in the paper's simulations)
SR_est <- function(fit, u = rep(1, length(fit$beta))) {
  fit <- add_eigen(fit)
  z <- drop(crossprod(fit$U, u))
  a <- SR_a(z, fit$lambda)
  beta_u <- sum(u * fit$beta)
  s <- SR_c(SR_mu(fit$sigma2, beta_u), a$a1)
  Sinv_u <- drop(fit$U %*% (z / fit$lambda))
  fit$beta - s * beta_u * Sinv_u
}

# GSR: one slab per eigenvector u_l of X'X with mu_l = sigma^2 / (u_l'beta)^2.
# The shrinkage on u_l is (mu_l / lambda_l) / (1 + mu_l / lambda_l).
GSR_est <- function(fit) {
  fit <- add_eigen(fit)
  bu <- drop(crossprod(fit$U, fit$beta))
  s <- fit$sigma2 / (fit$sigma2 + fit$lambda * bu^2)
  fit$beta - drop(fit$U %*% (s * bu))
}


# Linear shrinkage: LSh -------------------------------------------------------

# LSh works on centered data without the intercept (Sigma_tilde = diag(Sigma)).
# The intercept is put back as mean(y) - xbar'beta so the output has the same
# length as the other estimators.
LSh_est <- function(X, y) {
  Z <- X[, -1, drop = FALSE]
  xbar <- colMeans(Z)
  Zc <- sweep(Z, 2, xbar)
  fit <- ols_parts(Zc, y - mean(y), eigen = FALSE)
  # centering uses one more degree of freedom, so sigma^2 = RSS / (n - p - 1)
  sigma2 <- fit$sigma2 * (nrow(Zc) - ncol(Zc)) / (nrow(Zc) - ncol(Zc) - 1)
  b <- fit$beta
  S <- fit$XtX
  d <- diag(S)
  t1 <- sigma2 * sum(1 / d)
  t2 <- sigma2 * sum(diag(chol2inv(chol(S))))
  # t3 = beta' M'M beta with M = Sigma_tilde^{-1} Sigma - I, eq. (A.9)
  Mb <- drop(S %*% b) / d - b
  t3 <- sum(Mb^2)
  rho <- (t2 - t1) / (t2 - t1 + t3)
  slope <- b + rho * Mb
  c(mean(y) - sum(xbar * slope), slope)
}


# Shrinkage ridge regression: SRR ---------------------------------------------

# H(rho) of Theorem 4, with V_k = (y'X u_k)^2, V_star the mean eigenvalue and
# df = n - p - 1
SRR_H <- function(rho, lambda, V, V_star, yty, df) {
  d <- (1 - rho) * lambda + rho * V_star
  (1 / df) * (yty - 2 * sum(V / d) + sum(lambda * V / d^2)) * sum(lambda / d^2) +
    sum(V * rho^2 * (lambda - V_star)^2 / d^4)
}

SRR_rho <- function(fit) {
  fit <- add_eigen(fit)
  w <- drop(crossprod(fit$U, fit$Xty))
  m <- length(fit$lambda)
  V_star <- (1 / m) * sum(fit$lambda)
  optimize(
    SRR_H, c(0, 1),
    lambda = fit$lambda, V = w^2, V_star = V_star, yty = fit$yty, df = fit$n - m,
    tol = sqrt(.Machine$double.eps)
  )$minimum
}

# SRR = sum_k (y'X u_k) / ((1 - rho*) lambda_k + rho* V_star) u_k, eq. (19)
SRR_est <- function(fit, rho = SRR_rho(fit)) {
  fit <- add_eigen(fit)
  w <- drop(crossprod(fit$U, fit$Xty))
  V_star <- (1 / length(fit$lambda)) * sum(fit$lambda)
  drop(fit$U %*% (w / ((1 - rho) * fit$lambda + rho * V_star)))
}
