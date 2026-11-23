# Tables D7-D9 from the files in results/. Each entry is the test MSE averaged
# over the N = 50 splits, four decimals. The smallest value in each column,
# after rounding, is in red (ties included).
#
# Output: tables/Table_D7.tex, Table_D8.tex and Table_D9.tex, with the
# unrounded means in tables/Table_D7_mse.csv etc.

source("../R/latex_tables.R")

models <- c("OLS", "RR", "EN", "LASSO", "St", "DSh", "Sh", "SR", "GSR")
fms_pairs <- c(1, 2, 3, 6)

mean_mse <- function(file) {
  if (!file.exists(file)) {
    warning("missing ", file)
    return(setNames(rep(NA_real_, length(models)), models))
  }
  colMeans(read.csv(file)[models])
}

body_rows <- function(x) {
  xr <- round(x, 4)
  cells <- highlight(matrix(sprintf("%.4f", xr), nrow(xr)), col_is_min(xr))
  paste0("\\textbf{", models, "} & ", apply(cells, 1, paste, collapse = " & "), " \\\\")
}

save_means <- function(x, file) {
  write.csv(data.frame(Model = models, x, check.names = FALSE), file, row.names = FALSE)
}

table_end <- c("\\bottomrule", "\\end{tabular*}", "}", "\\end{table}")
dir.create("tables", showWarnings = FALSE)


# Tables D7 and D8 ------------------------------------------------------------

logit_means <- function(input) {
  key <- expand.grid(M = c(1, 10), fms = fms_pairs)
  x <- sapply(sprintf("results/logit_%s_M%d_FMS_0_vs_%d.csv", input, key$M, key$fms), mean_mse)
  colnames(x) <- sprintf("FMS_0_vs_%d_M%d", key$fms, key$M)
  x
}

logit_table <- function(x, caption, label) {
  c(
    "\\begin{table}[ht!]",
    "\\setlength{\\tabcolsep}{0pt} % Adjust column spacing",
    "\\renewcommand{\\arraystretch}{1}",
    paste0("\\caption{\\textbf{", caption, "}}"),
    paste0("\\label{", label, "}"),
    "\\footnotesize{",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}ccccccccc@{}}",
    "\\toprule",
    paste0(" & ", paste0("\\multicolumn{2}{c}{\\textbf{FMS ($0$ vs. $", fms_pairs, "$)}}", collapse = " & "), "\\\\"),
    " \\cmidrule(r){2-3} \\cmidrule(r){4-5} \\cmidrule(r){6-7} \\cmidrule(r){8-9}",
    paste0("\\multirow{-2}{*}{\\textbf{Model}} & ", paste(rep(c("$M=1$", "$M=10$"), 4), collapse = " & "), "\\\\"),
    "\\midrule",
    body_rows(x),
    table_end
  )
}

x <- logit_means("raw")
write_tex(
  logit_table(x, "Estimated MSE for binary classifications with raw input data",
              "tab:results_Binomial_log_regression_raw_data"),
  "tables/Table_D7.tex"
)
save_means(x, "tables/Table_D7_mse.csv")

x <- logit_means("standardized")
write_tex(
  logit_table(x, "Estimated MSE for binary classifications with standardized input data",
              "tab:results_Binomial_log_regression_standardized_data"),
  "tables/Table_D8.tex"
)
save_means(x, "tables/Table_D8_mse.csv")


# Table D9 --------------------------------------------------------------------

key <- expand.grid(M = c(1, 10), input = c("raw", "standardized"), stringsAsFactors = FALSE)
x <- sapply(sprintf("results/poisson_%s_M%d.csv", key$input, key$M), mean_mse)
colnames(x) <- sprintf("%s_M%d", key$input, key$M)

write_tex(
  c(
    "\\begin{table}[htbp!]",
    "\\centering",
    "\\caption{\\textbf{Estimated MSE for Poisson GLM with raw and standardized input data}}",
    "\\label{tab:results_Poisson_log_regression_raw_data}",
    "\\setlength{\\tabcolsep}{8pt}",
    "\\renewcommand{\\arraystretch}{1}",
    "\\footnotesize{",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}ccccc@{}}",
    "\\toprule",
    " & \\multicolumn{2}{c}{\\textbf{Raw Data}} & \\multicolumn{2}{c}{\\textbf{Standardized Data}}\\\\",
    " \\cmidrule{2-3} \\cmidrule{4-5}",
    "\\multirow{-2}{*}{\\textbf{Model}} & $M=1$ & $M=10$ & $M=1$ & $M=10$ \\\\",
    "\\midrule",
    body_rows(x),
    table_end
  ),
  "tables/Table_D9.tex"
)
save_means(x, "tables/Table_D9_mse.csv")
