# Table B5 from the timing files in results/. Each entry is the mean time of
# a procedure divided by the mean time of OLS for the same setting. The lowest
# value in each column is in red and the second lowest is underlined (ranked
# on the unrounded ratios).
#
# savvySh returns St, DSh and Sh from one call, so the table has one row for
# them (St/DSh/Sh, timed as T_St_DSh_Sh).
#
# Output: tables/Table_B5.tex and tables/Table_B5_ratios.csv

source("../R/latex_tables.R")

n <- 500
sigmas <- c(1, 5)
p_n_ratios <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
procedures <- c(
  RR = "T_RR", EN = "T_EN", LASSO = "T_LASSO", "$\\text{POET}_{99\\%}$" = "T_POET",
  "St/DSh/Sh" = "T_St_DSh_Sh", "SR/GSR" = "T_SR_GSR", SRR = "T_SRR"
)

# Printed as 5.40, 83.9, 930.8, and 1.1e3 or 2e3 for values of 1,000 or more
fmt_ratio <- function(x) {
  ifelse(
    x < 10, sprintf("%.2f", x),
    ifelse(x < 1000, sprintf("%.1f", x), {
      e <- floor(log10(x))
      paste0(sub("\\.0$", "", sprintf("%.1f", x / 10^e)), "e", e)
    })
  )
}

block_ratios <- function(sigma, ratios) {
  out <- NULL
  for (ratio in ratios) for (rho in rho_vals) {
    file <- sprintf("results/runtime_sigma_%g_rho_%g_n_%d_p_%d.csv", sigma, rho, n, round(n * ratio))
    if (!file.exists(file)) {
      warning("missing ", file)
      out <- cbind(out, rep(NA, length(procedures)))
      next
    }
    d <- read.csv(file)
    out <- cbind(out, colMeans(d[procedures]) / mean(d$T_OLS))
  }
  rownames(out) <- names(procedures)
  out
}

rank_cells <- function(x) {
  cells <- matrix(fmt_ratio(x), nrow(x))
  for (j in seq_len(ncol(x))) {
    if (all(is.na(x[, j]))) next
    o <- order(x[, j])
    cells[o[1], j] <- red(cells[o[1], j])
    cells[o[2], j] <- uline(cells[o[2], j])
  }
  cells
}

rho_head <- sapply(seq_along(rep(rho_vals, 3)), function(j) {
  bar <- if (j %in% c(7, 14)) "c|" else "c"
  paste0("\\multicolumn{1}{", bar, "}{\\textbf{", rep(rho_vals, 3)[j], "}}")
})

lines <- c(
  "\\begin{sidewaystable}",
  "\\centering",
  "\\renewcommand{\\arraystretch}{1.15}",
  "\\setlength{\\tabcolsep}{1.3pt}",
  "\\caption{\\textbf{Computational time relative to OLS for regression models with $n=500$}}",
  "\\label{tab:average_CT}",
  "\\scriptsize{",
  "\\begin{tabular}{c|ccccccc|ccccccc|ccccccc}",
  "\\hline",
  paste0("\\multicolumn{1}{c|}{\\textbf{$\\rho$}} & ", paste(rho_head, collapse = " & "), " \\\\"),
  "\\hline"
)
csv <- NULL
for (i in seq_along(sigmas)) {
  lines <- c(
    lines, "",
    paste0("\\multicolumn{22}{l}{\\textbf{Panel ", LETTERS[i], ": Normal Distribution with $\\sigma=", sigmas[i], "$}} \\\\"),
    "\\hline"
  )
  for (b in 1:2) {
    ratios <- p_n_ratios[3 * (b - 1) + 1:3]
    x <- block_ratios(sigmas[i], ratios)
    cells <- rank_cells(x)
    lines <- c(
      lines,
      paste0(
        "\\multicolumn{1}{c|}{\\textbf{$p/n$}} & \\multicolumn{7}{c|}{\\textbf{", pct(ratios[1]),
        "}} & \\multicolumn{7}{c|}{\\textbf{", pct(ratios[2]),
        "}} & \\multicolumn{7}{c}{\\textbf{", pct(ratios[3]), "}} \\\\"
      ),
      "\\hline",
      tex_row("%OLS", rep("1.00", 21), gap = FALSE),
      sapply(seq_len(nrow(x)), function(r) tex_row(rownames(x)[r], cells[r, ], gap = FALSE)),
      "\\hline"
    )
    key <- expand.grid(rho = rho_vals, p_n = ratios)
    csv <- rbind(csv, data.frame(sigma = sigmas[i], p_n = key$p_n, rho = key$rho, t(x), check.names = FALSE))
  }
}
# the sentence on compact notation only when a value of 1,000 or more is printed
compact <- if (any(csv[names(procedures)] >= 1000, na.rm = TRUE)) {
  "Large values are written in compact notation; for example, 4e3 denotes $4\\times10^3$ times the OLS computational time."
}
lines <- c(
  lines,
  "\\end{tabular}%",
  "}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  paste(collapse = " ", c(
    "\\item \\emph{Notes:} This table reports computational time relative to OLS across $N=250$ replications with $n=500$.",
    "For each setting, the average computational time of each estimator is divided by the average computational time of OLS.",
    "A value of $25$, for example, indicates that the estimator requires 25 times the OLS computational time.",
    "RR, EN, and LASSO include 10-fold cross-validation.",
    "St/DSh/Sh denotes the common multiplicative shrinkage call, which returns all three estimators,",
    "and SR/GSR denotes the common slab regression call.",
    "The lowest value in each setting is highlighted in \\textcolor{red}{red} and the second lowest is \\underline{underlined}.",
    compact
  )),
  "\\end{tablenotes}",
  "\\end{sidewaystable}"
)

dir.create("tables", showWarnings = FALSE)
write_tex(lines, "tables/Table_B5.tex")
write.csv(csv, "tables/Table_B5_ratios.csv", row.names = FALSE)
