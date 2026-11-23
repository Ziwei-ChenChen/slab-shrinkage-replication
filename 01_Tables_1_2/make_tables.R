# Tables 1 and 2 from the replication files in results/.
#
# Table 1: number of replications in which each estimator has the smallest
#          L2 error (every estimator within 1e-10 of the minimum is counted).
#          The largest count in each column is in red.
# Table 2: mean over replications of L2(method) / L2(SRR), three decimals.
#          Ratios above 1 are in red.
#
# Output: tables/Table_1.tex, tables/Table_2.tex and the same numbers as CSV.

source("../R/latex_tables.R")

n <- 1000
sigmas <- c(1, 5)
p_n_ratios <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.95)
f_p_ratios <- c(0.25, 0.50, 0.75)
models <- c("OLS", "RR", "EN", "LASSO", "SRR", "POET_99", "POET_f")
labels <- c(
  OLS = "\\textbf{OLS}", RR = "\\textbf{RR}", EN = "\\textbf{EN}", LASSO = "\\textbf{LASSO}",
  SRR = "\\textbf{SRR}", POET_99 = "\\textbf{$\\text{POET}_{\\bm{99\\%}}$}",
  POET_f = "\\textbf{$\\text{POET}_{\\bm{f}}$}"
)
tol <- 1e-10

settings <- expand.grid(f_p = f_p_ratios, p_n = p_n_ratios, sigma = sigmas)
settings$p <- pmax(1, round(n * settings$p_n))
settings$f <- ceiling(settings$p * settings$f_p)

counts <- ratios <- matrix(NA, length(models), nrow(settings), dimnames = list(models, NULL))
for (j in seq_len(nrow(settings))) {
  s <- settings[j, ]
  file <- sprintf("results/latent_sigma_%g_n_%d_p_%d_f_%d.csv", s$sigma, n, s$p, s$f)
  if (!file.exists(file)) {
    warning("missing ", file)
    next
  }
  L2 <- as.matrix(read.csv(file)[paste0("L2_", models)])
  colnames(L2) <- models
  best <- abs(L2 - apply(L2, 1, min, na.rm = TRUE)) <= tol
  counts[, j] <- colSums(best, na.rm = TRUE)
  ratios[, j] <- colMeans(L2 / L2[, "SRR"], na.rm = TRUE)
}

out <- cbind(settings[c("sigma", "p_n", "f_p", "p", "f")], t(counts))
dir.create("tables", showWarnings = FALSE)
write.csv(out, "tables/Table_1_counts.csv", row.names = FALSE)
out <- cbind(settings[c("sigma", "p_n", "f_p", "p", "f")], t(ratios))
write.csv(out, "tables/Table_2_ratios.csv", row.names = FALSE)

block_header <- function(p_n) {
  c(
    paste0(
      "\\multirow{2}{*}{\\textbf{Model}} & $\\bm{p/n}$ & ",
      paste0("\\multicolumn{3}{c}{\\textbf{", pct(p_n), "}}", collapse = " & "), " \\\\"
    ),
    "\\cmidrule{3-5}",
    "\\cmidrule{6-8}",
    "\\cmidrule{9-11}",
    paste0("\\multirow{2}{*}{} & $\\bm{f/p}$ & ", bold_cells(rep(pct(f_p_ratios), 3)), " \\\\"),
    "\\midrule"
  )
}

body_lines <- function(cells, rows) {
  sapply(rows, function(m) tex_row(labels[m], cells[m, ]))
}

# cells: formatted character matrix, models x settings
table_lines <- function(cells, rows) {
  lines <- character(0)
  for (sigma in sigmas) {
    lines <- c(
      lines,
      paste0("\\multicolumn{11}{c}{\\textbf{Normal Distribution: $\\bm\\sigma=\\textbf{", sigma, "}$}}\\\\"),
      "\\midrule"
    )
    for (b in 1:2) {
      take <- which(settings$sigma == sigma & settings$p_n %in% p_n_ratios[3 * (b - 1) + 1:3])
      lines <- c(lines, block_header(p_n_ratios[3 * (b - 1) + 1:3]), body_lines(cells[, take], rows))
      if (!(sigma == max(sigmas) && b == 2)) lines <- c(lines, "\\midrule")
    }
  }
  lines
}

wrap <- function(body, stretch, caption, label, note) {
  c(
    "\\begin{table}[ht!]",
    "\\setlength{\\tabcolsep}{0pt}",
    paste0("\\renewcommand{\\arraystretch}{", stretch, "}"),
    "\\Large",
    paste0("\\caption{\\textbf{", caption, "}}"),
    paste0("\\label{", label, "}"),
    "\\footnotesize{",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}ccccccccccc@{}}",
    "\\toprule",
    body,
    "\\bottomrule",
    "\\end{tabular*}%",
    "}",
    "\\begin{tablenotes}[flushleft]",
    "\\footnotesize",
    paste0("\\item \\emph{Note:} ", note),
    "\\end{tablenotes}",
    "\\end{table}"
  )
}

# Table 1
cells <- matrix(as.character(counts), nrow(counts), dimnames = dimnames(counts))
flag <- matrix(FALSE, nrow(counts), ncol(counts), dimnames = dimnames(counts))
flag[] <- col_is_max(counts)
cells <- highlight(cells, flag)
note1 <- paste(
  "We tabulate counts of how many times each estimator (OLS, RR, EN, LASSO, SRR, $\\text{POET}_{\\bm{99\\%}}$, and $\\text{POET}_{\\bm{f}}$)",
  "achieves the lowest $L_2$-distance (from the ``true'' regression parameters) across $N = 250$ replications of samples of size $n = 1,000$",
  "for various choices of $p/n$ and $f/p$ under $\\sigma=1, 5$. The best-performing method is highlighted in \\textcolor{red}{red}. "
)
write_tex(
  wrap(table_lines(cells, models), "0.8", "Counts of models achieving the minimum $L_2$ error", "counts_mean_l2_error_SRR", note1),
  "tables/Table_1.tex"
)

# Table 2
rows <- setdiff(models, "SRR")
cells <- matrix(sprintf("%.3f", ratios), nrow(ratios), dimnames = dimnames(ratios))
cells <- highlight(cells, ratios > 1)
note2 <- paste(
  "We report the mean $L_2$ error ratios of each estimator (OLS, RR, EN, LASSO, $\\text{POET}_{\\bm{99\\%}}$, and $\\text{POET}_{\\bm{f}}$)",
  "relative to SRR across $N = 250$ simulations with sample size $n = 1,000$, and $\\sigma =\\{1,5\\}$.",
  "A ratio greater/lower than $1$ indicates ``worse''/``better'' performance than SRR.",
  "Instances highlighted in \\textcolor{red}{red} indicate cases where SRR outperforms the comparator (i.e., those with ratio larger than $1$)."
)
write_tex(
  wrap(table_lines(cells, rows), "1", "Mean $L_2$ error ratios of estimators relative to SRR", "mean_l2_error_ratios_SRR", note2),
  "tables/Table_2.tex"
)
