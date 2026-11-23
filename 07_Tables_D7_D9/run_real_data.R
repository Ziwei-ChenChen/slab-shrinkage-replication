# Tables D7-D9: logistic and Poisson GLMs for the cyber-sickness data
# (Appendix D). Each GLM is fitted by IRLS, with the weighted least squares
# step done by OLS or by one of St, DSh, Sh, SR and GSR. RR, EN and LASSO are
# fitted with cv.glmnet. Every replication splits the data 70/30 at random and
# records the test MSE of the nine models.
#
#   Table D7  logistic GLM, raw features, FMS 0 against FMS 1, 2, 3 and 6
#   Table D8  the same with standardized features
#   Table D9  Poisson GLM for the four severity classes, raw and standardized
#
# M = 1 uses the 13 features at lag t-1 and M = 10 the lags t-1, ..., t-10
# (130 features). In both data files the first column is the intercept and
# the last one is the FMS score at time t.
#
# Run from this folder:  Rscript run_real_data.R

source("../R/All_eight_estimators.R")
source("../R/benchmarks.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")
source("glm_irls.R")

N <- 50
test_frac <- 0.3
fms_pairs <- c(1, 2, 3, 6)
models <- c("OLS", "RR", "EN", "LASSO", "St", "DSh", "Sh", "SR", "GSR")

dir.create("results", showWarnings = FALSE)
start_workers()

# Features are standardized over the whole file, before any subsetting. The
# 67 scores of 2.5 and the 74 of 3.5 are rounded down to 2 and 3, as described
# in Appendix D. All other scores are integers.
read_cyber <- function(M, standardize) {
  d <- read.csv(sprintf("data/data_for_regression%d.csv", M), check.names = FALSE)
  X <- as.matrix(d[, -ncol(d)])
  if (standardize) X[, -1] <- scale(X[, -1])
  list(X = X, fms = floor(d[[ncol(d)]]))
}

# Severity classes for the Poisson GLM: FMS 0 -> 0, 1 to 3 -> 1,
# 4 to 6 -> 2 and 7 to 9 -> 3
severity_class <- function(fms) findInterval(fms, c(1, 4, 7))

# Predicted label from the linear predictor
predict_label <- function(eta, family) {
  if (family == "binomial") as.integer(plogis(eta) >= 0.5) else pmin(floor(exp(eta)), 3)
}

# One train/test split. The logistic fits start IRLS from the least squares
# fit of y on X, and the Poisson fits start from mu = y + 0.1. glmnet
# standardizes the features internally (its default), as in the simulations.
glm_rep <- function(r, X, y, family) {
  seed_rep(123 + r, kind = "Mersenne-Twister")
  idx <- sample.int(nrow(X))
  cut <- floor(nrow(X) * (1 - test_frac))
  train <- idx[1:cut]
  test <- idx[(cut + 1):length(idx)]
  Xtr <- X[train, ]
  ytr <- y[train]
  Xte <- X[test, ]
  yte <- y[test]

  start <- NULL
  if (family$family == "binomial") start <- drop(solve(crossprod(Xtr), crossprod(Xtr, ytr)))
  irls <- c("OLS", "St", "DSh", "Sh", "SR", "GSR")
  eta <- lapply(setNames(irls, irls), function(m) drop(Xte %*% irls_fit(Xtr, ytr, family, m, start)))

  seed_rep(10000 + r, kind = "Mersenne-Twister")
  foldid <- make_folds(length(ytr))
  for (m in c("RR", "EN", "LASSO")) {
    alpha <- c(RR = 0, EN = 0.5, LASSO = 1)[[m]]
    cv <- cv.glmnet(Xtr[, -1], ytr, family = family$family, alpha = alpha, foldid = foldid,
                    nlambda = 100)
    eta[[m]] <- drop(predict(cv, Xte[, -1], s = "lambda.min", type = "link"))
  }

  mse <- sapply(eta[models], function(e) mean((yte - predict_label(e, family$family))^2))
  data.frame(replication = r, t(mse))
}

for (input in c("raw", "standardized")) {
  for (M in c(1, 10)) {
    dat <- read_cyber(M, standardize = input == "standardized")

    for (fms in fms_pairs) {
      keep <- dat$fms %in% c(0, fms)
      X <- dat$X[keep, ]
      y <- as.integer(dat$fms[keep] == fms)
      file <- sprintf("results/logit_%s_M%d_FMS_0_vs_%d.csv", input, M, fms)
      run_case(file, N, function(r) glm_rep(r, X, y, binomial()))
    }

    file <- sprintf("results/poisson_%s_M%d.csv", input, M)
    run_case(file, N, function(r) glm_rep(r, dat$X, severity_class(dat$fms), poisson()))
  }
}

plan(sequential)
