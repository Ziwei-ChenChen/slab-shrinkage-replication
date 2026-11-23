# Functions for Table 4: estimators, population and traits, flashfm output,
# one loop of flashfmZero, and power and FDR.

# flashfmZero commit 009fe64 (13 February 2025)
flashfmZero_commit <- "009fe649e9ecf6020671f375c62139d0613f29ca"

# Causal variants (Appendix C): trait 1 has C1 and C2, trait 2 has C1 and C3
causal <- list(
  T1 = c("rs61839660", "rs62626317"),
  T2 = c("rs61839660", "rs11594656")
)

# flashfm-savvySh calls selreg(x, y) with the genotypes only and expects the
# coefficients with the intercept first
selreg_from <- function(est) {
  force(est)
  function(x, y) est(ols_parts(cbind(1, as.matrix(x)), y))
}

estimators <- list(
  OLS = selreg_from(OLS_est),
  St = selreg_from(St_est),
  DSh = selreg_from(DSh_est),
  GSR = selreg_from(GSR_est)
)
methods <- names(estimators)

# Population matrix of genotype scores (0, 1, 2) with rsIDs as column names.
# Variants with MAF <= 0.005 are dropped, as in Appendix C.
read_population <- function(file) {
  G <- as.matrix(readRDS(file))
  storage.mode(G) <- "double"
  if (is.null(colnames(G))) stop("the population matrix needs rsIDs as column names")
  missing <- setdiff(unlist(causal), colnames(G))
  if (length(missing)) stop("causal variants not in the population: ", paste(missing, collapse = ", "))
  af <- colMeans(G) / 2
  G[, pmin(af, 1 - af) > 0.005, drop = FALSE]
}

# n individuals drawn from the population, causal effects from U(0.15, 0.4)
# for each trait, and bivariate normal errors with variance 0.20 and
# correlation 0.40. The traits have no intercept.
draw_traits <- function(G, n = 5000, effect = c(0.15, 0.4), var_e = 0.2, cor_e = 0.4) {
  X <- G[sample.int(nrow(G), n), , drop = FALSE]
  rownames(X) <- NULL
  beta <- lapply(causal, function(snps) setNames(runif(length(snps), effect[1], effect[2]), snps))
  Sigma <- var_e * matrix(c(1, cor_e, cor_e, 1), 2)
  E <- matrix(rnorm(2 * n), n) %*% chol(Sigma)
  Y <- sapply(names(causal), function(t) drop(X[, causal[[t]]] %*% beta[[t]])) + E
  colnames(Y) <- names(causal)
  list(X = X, Y = Y, beta = beta)
}

# mpp.pp$MPP has one matrix per trait, with the single-trait and multi-trait
# MPP of each listed variant. Variants that are not listed have MPP 0.
mpp_long <- function(mpp, snps, beta) {
  do.call(rbind, lapply(names(mpp), function(t) {
    m <- as.matrix(mpp[[t]])
    m <- m[rownames(m) %in% snps, , drop = FALSE]
    data.frame(
      trait = rep(t, nrow(m)), snp = rownames(m),
      mpp_single = unname(m[, 1]), mpp_multi = unname(m[, 2]),
      causal = rownames(m) %in% causal[[t]],
      effect = unname(beta[[t]][rownames(m)]),
      stringsAsFactors = FALSE
    )
  }))
}

# merger() in groupmultiU() of flashfmZero loops over the pairs of groups to
# merge with for (k in 1:nrow(wh)). This sets the loop to seq_len(nrow(wh)),
# which is the same when there is a pair to merge and skipped when there is
# none. The installed package is not changed, and the variant MPPs
# (mpp.pp$MPP) do not depend on the groups.
adapt_flashfmZero <- function() {
  f <- get("groupmultiU", envir = asNamespace("flashfmZero"))
  n <- 0L
  swap <- function(e) {
    if (!is.call(e)) return(e)
    if (identical(e[[1]], as.name("for")) && identical(e[[3]], quote(1:nrow(wh)))) {
      e[[3]] <- quote(seq_len(nrow(wh)))
      n <<- n + 1L
    }
    for (i in seq_along(e)) if (is.call(e[[i]])) e[[i]] <- swap(e[[i]])
    e
  }
  new_body <- swap(body(f))
  if (n == 0L) {
    if (!any(grepl("for (k in seq_len(nrow(wh)))", deparse(f), fixed = TRUE))) {
      stop("groupmultiU() of flashfmZero is not the expected version")
    }
    return(invisible(FALSE))
  }
  if (n > 1L) stop("groupmultiU() of flashfmZero is not the expected version")
  body(f) <- new_body
  utils::assignInNamespace("groupmultiU", f, ns = "flashfmZero")
  invisible(TRUE)
}

# A variant is selected when its MPP is above the threshold. In each
# replication, power is the share of the two causal variants selected and FDR
# the share of selected variants that are not causal (0 if none is selected).
# Both are averaged over the replications. runs lists every replication,
# method and trait, so a run with no variant listed counts as 0.
power_fdr <- function(res, runs, threshold = 0.9) {
  out <- NULL
  for (type in c("single", "multi")) {
    res$n_sel <- res[[paste0("mpp_", type)]] > threshold
    res$tp <- res$n_sel & res$causal
    per_rep <- aggregate(cbind(n_sel, tp) ~ replication + method + trait, data = res, FUN = sum)
    per_rep <- merge(runs, per_rep, all.x = TRUE)
    per_rep[is.na(per_rep)] <- 0
    per_rep$power <- per_rep$tp / lengths(causal[per_rep$trait])
    per_rep$fdr <- ifelse(per_rep$n_sel > 0, (per_rep$n_sel - per_rep$tp) / per_rep$n_sel, 0)
    avg <- aggregate(cbind(power, fdr) ~ method + trait, data = per_rep, FUN = mean)
    avg$replications <- aggregate(power ~ method + trait, data = per_rep, FUN = length)$power
    avg$type <- type
    out <- rbind(out, avg)
  }
  out
}
