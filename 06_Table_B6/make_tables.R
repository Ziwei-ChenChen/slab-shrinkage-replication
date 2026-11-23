# Table B6 from the replication files in results/. Each entry is the mean of
# the squared error ||beta_hat - beta||^2 over the N = 500 replications, three
# decimals. The smallest value in each column is in red.
#
# Output: tables/Table_B6.tex and tables/Table_B6_mse.csv

source("../R/latex_tables.R")

n <- 500
sigmas <- c(1, 5)
p_n_ratios <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95)
rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
rows <- c(ORR = "SE2_ORR", "OSR$(1)$" = "SE2_OSR1", "OSR$(\\mathcal{U})$" = "SE2_OSRU")

# Rounded to four decimals (half up) and then to three
fmt3 <- function(x) sprintf("%.3f", round(round(x, 4) + 1e-9, 3))

block_mse <- function(sigma, ratios) {
  out <- NULL
  for (ratio in ratios) for (rho in rho_vals) {
    file <- sprintf("results/oracle_sigma_%g_rho_%g_n_%d_p_%d.csv", sigma, rho, n, round(n * ratio))
    if (!file.exists(file)) {
      warning("missing ", file)
      out <- cbind(out, rep(NA, length(rows)))
      next
    }
    out <- cbind(out, colMeans(read.csv(file)[rows]))
  }
  rownames(out) <- names(rows)
  out
}

lines <- c(
  "\\begin{sidewaystable}",
  "\\setlength{\\tabcolsep}{1pt} % Adjust column spacing",
  "\\renewcommand{\\arraystretch}{1.5} % Adjust row spacing",
  "\\caption{\\textbf{Oracle MSE comparison of Ridge Regression and Slab Regression}}",
  "\\label{tab:oracle_comparison}",
  "\\huge",
  "\\scriptsize{",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}cc|ccccccc|ccccccc@{}}",
  "\\hline",
  paste0("\\multicolumn{1}{c}{\\textbf{}} & \\multicolumn{1}{c|}{\\textbf{$\\rho$}} & ", bold_cells(rep(rho_vals, 2)), " \\\\"),
  "\\hline"
)
csv <- NULL
for (i in seq_along(sigmas)) {
  lines <- c(
    lines,
    paste0("\\multicolumn{16}{l}{\\textbf{Panel ", LETTERS[i], ": Normal Distribution with $\\sigma = ", sigmas[i], "$}} \\\\"),
    "\\hline"
  )
  for (b in 1:3) {
    ratios <- p_n_ratios[2 * (b - 1) + 1:2]
    x <- block_mse(sigmas[i], ratios)
    cells <- highlight(matrix(fmt3(x), nrow(x)), col_is_min(x))
    lines <- c(
      lines,
      paste0(
        "\\multicolumn{1}{c}{\\textbf{}} & \\multicolumn{1}{c|}{\\textbf{$p/n$}} & \\multicolumn{7}{c|}{\\textbf{",
        pct(ratios[1]), "}} & \\multicolumn{7}{c}{\\textbf{", pct(ratios[2]), "}} \\\\"
      ),
      "\\hline",
      sapply(seq_len(nrow(x)), function(r) tex_row(rownames(x)[r], cells[r, ])),
      "\\hline"
    )
    key <- expand.grid(rho = rho_vals, p_n = ratios)
    csv <- rbind(csv, data.frame(sigma = sigmas[i], p_n = key$p_n, rho = key$rho, t(x), check.names = FALSE))
  }
}
lines <- c(
  lines,
  "\\end{tabular*}",
  "}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  paste(
    "\\item \\emph{Notes:} This table reports simulation-based estimates of the parameter MSE,",
    "$N^{-1}\\sum_{r=1}^{N}\\|\\widehat{\\bm{\\beta}}^{(r)}-\\bm{\\beta}\\|_2^2$, for Oracle Ridge Regression (ORR),",
    "Oracle Slab Regression with the fixed direction $\\mathbf{u}=\\mathbf{1}$, denoted by OSR$(1)$,",
    "and Oracle Slab Regression with the direction selected over the finite candidate set $\\mathcal{U}$, denoted by OSR$(\\mathcal{U})$.",
    "The covariates follow the Gaussian Toeplitz design in \\eqref{dgp_x_1}, with $n=500$ and $\\rho\\in\\{-0.75,-0.50,-0.25,0,0.25,0.50,0.75\\}$.",
    "Each configuration is based on $N=500$ replications.",
    "Panel A considers Gaussian errors with $\\sigma=1$ and Panel B considers Gaussian errors with $\\sigma=5$.",
    "The smallest MSE in each configuration is highlighted in \\textcolor{red}{red}."
  ),
  "\\end{tablenotes}",
  "\\end{sidewaystable}"
)

dir.create("tables", showWarnings = FALSE)
write_tex(lines, "tables/Table_B6.tex")
write.csv(csv, "tables/Table_B6_mse.csv", row.names = FALSE)
