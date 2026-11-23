# glm_irls.R
#
# IRLS for a GLM in which the weighted least squares step can be replaced by
# any of the shrinkage estimators. This is the algorithm of glm2::glm.fit2
# (step halving when the deviance is not finite, leaves the valid region, or
# goes up), with the call to lm.fit swapped for wls_step(). If that step fails
# (X'WX singular) or gives non-finite coefficients, IRLS stops at the last
# coefficients with a warning.
#
# x already contains the intercept column. `method` names the estimator used
# for the weighted least squares step (see wls_step below).

# The weighted least squares step for each GLM in Tables D7-D9
wls_step <- function(method, X, y) {
  switch(method,
    OLS = lm.fit(X, y, tol = 1e-9)$coefficients,
    St = St_est(ols_parts(X, y, eigen = FALSE)),
    DSh = DSh_est(ols_parts(X, y, eigen = FALSE)),
    Sh = Sh_est(ols_parts(X, y, eigen = FALSE)),
    SR = SR_est(ols_parts(X, y)),
    GSR = GSR_est(ols_parts(X, y)),
    stop("unknown method ", method)
  )
}

irls_fit <- function(x, y, family, method, start = NULL, maxit = 200, epsilon = 1e-6) {
  nobs <- NROW(y)
  weights <- rep(1, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  dev.resids <- family$dev.resids
  valideta <- if (is.null(family$valideta)) function(eta) TRUE else family$valideta
  validmu <- if (is.null(family$validmu)) function(mu) TRUE else family$validmu
  mustart <- NULL
  eval(family$initialize)

  coefold <- NULL
  if (is.null(start)) {
    eta <- family$linkfun(mustart)
  } else {
    coefold <- start
    eta <- drop(x %*% start)
  }
  mu <- linkinv(eta)
  devold <- sum(dev.resids(y, mu, weights))
  conv <- FALSE
  coef <- coefold

  update <- function(b) {
    eta <<- drop(x %*% b)
    mu <<- linkinv(eta)
    sum(dev.resids(y, mu, weights))
  }

  for (iter in 1:maxit) {
    d_eta <- mu.eta(eta)
    good <- d_eta != 0
    z <- eta[good] + (y - mu)[good] / d_eta[good]
    w <- sqrt(d_eta[good]^2 / variance(mu)[good])
    start <- tryCatch(wls_step(method, x[good, , drop = FALSE] * w, z * w), error = function(e) {
      warning(method, " step failed at iteration ", iter, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(start)) break
    start[is.na(start)] <- 0
    if (any(!is.finite(start))) {
      warning("non-finite coefficients at iteration ", iter)
      break
    }
    dev <- update(start)

    if (!is.finite(dev)) {
      if (is.null(coefold)) stop("no valid coefficients found; please supply starting values")
      ii <- 1
      while (!is.finite(dev)) {
        if (ii > maxit) stop("cannot correct step size")
        ii <- ii + 1
        start <- (start + coefold) / 2
        dev <- update(start)
      }
    }
    if (!(valideta(eta) && validmu(mu))) {
      if (is.null(coefold)) stop("no valid coefficients found; please supply starting values")
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > maxit) stop("cannot correct step size")
        ii <- ii + 1
        start <- (start + coefold) / 2
        dev <- update(start)
      }
    }
    if (iter > 1 && (dev - devold) / (0.1 + abs(dev)) >= epsilon) {
      ii <- 1
      while ((dev - devold) / (0.1 + abs(dev)) > -epsilon) {
        if (ii > maxit) break
        ii <- ii + 1
        start <- (start + coefold) / 2
        dev <- update(start)
      }
    }
    if (abs(dev - devold) / (0.1 + abs(dev)) < epsilon) {
      conv <- TRUE
      coef <- start
      break
    }
    devold <- dev
    coef <- coefold <- start
  }
  if (!conv) warning("IRLS did not converge in ", maxit, " iterations")
  coef
}
