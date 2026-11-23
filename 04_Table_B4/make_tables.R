# Table B4 from the replication files in results/. Each cell counts the
# replications in which the estimator has the smallest L2 error; the largest
# count in each column is in red. RR, EN and LASSO are NA for p = 1.
#
# Output: tables/Table_B4.tex and tables/Table_B4_counts.csv

source("../R/latex_tables.R")

sigmas <- c(1, 5)
n_vals <- c(1000, 2500, 5000)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
rows <- c("OLS", "RR", "EN", "LASSO", "St", "DSh", "Sh", "SR", "GSR")
columns <- rbind(data.frame(p = 1, rho = NA), expand.grid(rho = rho_vals, p = c(2, 5))[c("p", "rho")])

panel_counts <- function(sigma, n) {
  out <- sapply(seq_len(nrow(columns)), function(j) {
    file <- sprintf("results/binomial_sigma_%g_rho_%s_n_%d_p_%d.csv", sigma, columns$rho[j], n, columns$p[j])
    if (!file.exists(file)) {
      warning("missing ", file)
      return(rep(NA, length(rows)))
    }
    as.numeric(table(factor(read.csv(file)$winner, levels = rows))[rows])
  })
  rownames(out) <- rows
  out[c("RR", "EN", "LASSO"), columns$p == 1] <- NA
  out
}

lines <- c(
  "\\begin{sidewaystable}",
  "\\centering",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\renewcommand{\\arraystretch}{0.6}",
  "\\caption{\\textbf{Best performance regression model}}",
  "\\label{tab:simulated_generic_binomial}",
  "\\footnotesize{",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}c|c|ccccccc|ccccccc@{}}",
  "\\hline",
  paste0("$\\bm\\rho$ & \\textbf{NA} & ", bold_cells(rep(rho_vals, 2)), " \\\\"),
  "\\hline"
)
csv <- NULL
for (sigma in sigmas) {
  lines <- c(lines, paste0("\\multicolumn{16}{c}{$\\bm\\sigma=\\textbf{", sigma, "}$}\\\\"), "\\hline", "")
  for (i in seq_along(n_vals)) {
    counts <- panel_counts(sigma, n_vals[i])
    cells <- matrix(ifelse(is.na(counts), "NA", counts), nrow(counts))
    cells <- highlight(cells, col_is_max(counts))
    lines <- c(
      lines,
      paste0(
        "\\multicolumn{16}{l}{\\textbf{Panel ", LETTERS[i], ": $\\bm n=\\textbf{",
        format(n_vals[i], big.mark = ","), "}$}}\\\\"
      ),
      "\\hline",
      "$\\bm p$ & \\textbf{1} & \\multicolumn{7}{c|}{\\textbf{2}} & \\multicolumn{7}{c}{\\textbf{5}} \\\\",
      "\\hline",
      sapply(seq_along(rows), function(r) tex_row(rows[r], cells[r, ], gap = FALSE)),
      "\\hline",
      ""
    )
    csv <- rbind(csv, data.frame(sigma = sigma, n = n_vals[i], columns, t(counts)))
  }
}
lines <- c(
  lines[-length(lines)],
  "\\end{tabular*}",
  "}",
  "",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  paste(
    "\\item \\emph{Notes:} Performance comparisons across different estimators are shown for Gaussian dependent variables",
    "with $\\sigma\\in\\{1,5\\}$ and multivariate binomial covariates as in \\eqref{dgp_y_1} and \\eqref{dgp_x_2}, respectively.",
    "Estimators are ranked according to the $L_2$ distance from the true regression parameters across $N=250$ replications.",
    "Each entry reports the number of replications in which the corresponding estimator achieves the lowest $L_2$ distance,",
    "with the highest frequency for each scenario highlighted in \\textcolor{red}{red}.",
    "Results are presented for different sample sizes $n$, numbers of covariates $p$, and correlation parameters $\\rho$.",
    "For $p=1$, $\\rho$ is not applicable and is reported as NA."
  ),
  "\\end{tablenotes}",
  "\\end{sidewaystable}"
)

dir.create("tables", showWarnings = FALSE)
write_tex(lines, "tables/Table_B4.tex")
write.csv(csv, "tables/Table_B4_counts.csv", row.names = FALSE)
