# Tables B1, B2 and B3 from the replication files in results/.
# Each cell counts the replications in which the estimator has the smallest
# L2 error (ties go to the first estimator in the order OLS, St, DSh, Sh, SR,
# GSR, RR, EN, LASSO). The largest count in each column is in red.
#
# Output: tables/Table_B1.tex, Table_B2.tex, Table_B3.tex and the counts as CSV.

source("../R/latex_tables.R")

rho_vals <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
rows <- c("OLS", "RR", "EN", "LASSO", "St", "DSh", "Sh", "SR", "GSR")
panels <- list(
  list(n = 500, ratios = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95), head = "\\multicolumn{23}{l}{\\textbf{Panel A: $\\bm n=\\textbf{500}$}}\\\\"),
  list(n = 1000, ratios = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95), head = "\\multicolumn{23}{l}{\\textbf{Panel B: $\\bm{n}=\\textbf{1,000}$}} \\\\"),
  list(n = 2500, ratios = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95), head = "\\multicolumn{23}{l}{\\textbf{Panel C: $\\bm{n}=\\textbf{2,500}$}} \\\\")
)

dgp_note <- paste(
  "Estimators are ranked in $L_2$ distance (from the ``true'' regression parameters) across $N = 250$ replications,",
  "and the ``best'' estimator is highlighted in \\textcolor{red}{red} for various choices of $n$, $p/n$, and $\\rho$."
)
tables <- list(
  B1 = list(
    tag = "normal_sigma_1", label = "tab:simulated_normal_sigma1",
    title = "\\multicolumn{23}{c}{\\textbf{Normal Distribution}: $\\bm\\sigma=\\textbf{1}$}\\\\",
    note = paste(
      "Performance comparisons across different estimators are shown for Gaussian dependent variable with $\\sigma=1$",
      "and Gaussian covariates as in \\eqref{dgp_y_1} and \\eqref{dgp_x_1}, respectively.", dgp_note
    )
  ),
  B2 = list(
    tag = "normal_sigma_5", label = "tab:simulated_normal_sigma5",
    title = "\\multicolumn{23}{c}{\\textbf{Normal Distribution}: $\\bm\\sigma=\\textbf{5}$}\\\\",
    note = paste(
      "Performance comparisons across different estimators are shown for Gaussian dependent variable with $\\sigma=5$",
      "and Gaussian covariates as in \\eqref{dgp_y_1} and \\eqref{dgp_x_1}, respectively.", dgp_note
    )
  ),
  B3 = list(
    tag = "t_df_50over24", label = "tab:simulated_t",
    title = "\\multicolumn{23}{c}{\\textbf{$\\bm t$-distribution}: $\\bm v=\\textbf{50/24}$}\\\\",
    note = paste(
      "Performance comparisons across different estimators are shown for $t$-distribution dependent variable with $v=50/24$",
      "and Gaussian covariates as in \\eqref{dgp_y_2} and \\eqref{dgp_x_1}, respectively.", dgp_note
    )
  )
)

# counts for one p/n block: rows x (3 ratios x 7 rho)
block_counts <- function(tag, n, ratios) {
  out <- NULL
  for (ratio in ratios) for (rho in rho_vals) {
    p <- round(n * ratio)
    file <- sprintf("results/%s_rho_%g_n_%d_p_%d.csv", tag, rho, n, p)
    counts <- setNames(rep(NA, length(rows)), rows)
    if (file.exists(file)) {
      winner <- read.csv(file)$winner
      counts <- table(factor(winner, levels = rows))[rows]
    } else {
      warning("missing ", file)
    }
    out <- cbind(out, as.numeric(counts))
  }
  rownames(out) <- rows
  out
}

dir.create("tables", showWarnings = FALSE)
for (name in names(tables)) {
  tb <- tables[[name]]
  lines <- c(
    "\\begin{sidewaystable}",
    "\\setlength{\\tabcolsep}{1pt} % Adjust column spacing",
    "\\renewcommand{\\arraystretch}{0.7} % Adjust row spacing",
    "\\caption{\\textbf{Best performance regression model}}",
    paste0("\\label{", tb$label, "}"),
    "\\Huge",
    "\\footnotesize{",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}cc|ccccccc|ccccccc|ccccccc@{}}",
    "\\hline",
    paste0("& $\\bm\\rho$ & ", bold_cells(rep(rho_vals, 3)), " \\\\"),
    "\\hline",
    tb$title,
    "\\hline"
  )
  csv <- NULL
  for (panel in panels) {
    lines <- c(lines, panel$head, "\\hline")
    for (b in seq_len(length(panel$ratios) / 3)) {
      ratios <- panel$ratios[3 * (b - 1) + 1:3]
      counts <- block_counts(tb$tag, panel$n, ratios)
      cells <- highlight(matrix(as.character(counts), nrow(counts)), col_is_max(counts))
      lines <- c(
        lines,
        paste0(
          "\\multirow{2}{*}{} & $\\bm{p/n}$ & \\multicolumn{7}{c|}{\\textbf{", pct(ratios[1]),
          "}} & \\multicolumn{7}{c|}{\\textbf{", pct(ratios[2]),
          "}} & \\multicolumn{7}{c}{\\textbf{", pct(ratios[3]), "}} \\\\"
        ),
        "\\hline",
        sapply(seq_along(rows), function(i) tex_row(rows[i], cells[i, ])),
        "\\hline"
      )
      key <- expand.grid(rho = rho_vals, ratio = ratios)
      csv <- rbind(csv, data.frame(n = panel$n, p_n = key$ratio, rho = key$rho, t(counts)))
    }
  }
  lines <- c(
    lines,
    "\\end{tabular*}%",
    "}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    paste0("\\item \\emph{Notes:} ", tb$note),
    "\\end{tablenotes}",
    "\\end{sidewaystable}"
  )
  write_tex(lines, sprintf("tables/Table_%s.tex", name))
  write.csv(csv, sprintf("tables/Table_%s_counts.csv", name), row.names = FALSE)
}
