# Table 4 from the replication files in results/: power and FDR at the MPP
# threshold 0.9 for single-trait and multi-trait fine-mapping. The highest
# power for a trait among the single rows, and among the multi rows, is in red
# when it is above that of OLS.
#
# Output: tables/Table_4.tex and tables/Table_4.csv

source("../R/All_eight_estimators.R")
source("../R/latex_tables.R")
source("table4_functions.R")

N <- 100
traits <- names(causal)

files <- list.files("results", pattern = "^rep_[0-9]+\\.csv$", full.names = TRUE)
if (!length(files)) stop("no replications in results/")
res <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
res <- res[res$method %in% methods, ]
reps <- sort(unique(res$replication))
if (length(reps) < N) warning(length(reps), " of ", N, " replications in results/")

runs <- merge(unique(res[c("replication", "method")]), data.frame(trait = traits))
pf <- power_fdr(res, runs, threshold = 0.9)

# One row per method and type, in the order of the paper
rows <- expand.grid(type = c("single", "multi"), method = methods, stringsAsFactors = FALSE)[c("method", "type")]
value <- function(what, r, t) {
  pf[[what]][pf$method == rows$method[r] & pf$type == rows$type[r] & pf$trait == t]
}
power <- sapply(traits, function(t) sapply(seq_len(nrow(rows)), value, what = "power", t = t))
fdr <- sapply(traits, function(t) sapply(seq_len(nrow(rows)), value, what = "fdr", t = t))

best <- matrix(FALSE, nrow(rows), length(traits))
for (type in c("single", "multi")) {
  i <- which(rows$type == type)
  ols <- i[rows$method[i] == "OLS"]
  for (j in seq_along(traits)) {
    top <- i[power[i, j] == max(power[i, j])]
    best[top, j] <- power[top, j] > power[ols, j]
  }
}

num <- function(x) as.character(round(x, 3))
cells <- cbind(highlight(matrix(num(power), nrow(rows)), best), matrix(num(fdr), nrow(rows)))
labels <- paste0("\\textbf{", rows$type, "-", rows$method, "}")

lines <- c(
  "\\begin{table}[htbp!]",
  "\\caption{\\textbf{Power and FDR (false discovery rate) comparison for single and multi-trait fine-mapping based on four different estimators.}}",
  "\\label{tab:genetics_power}",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\renewcommand{\\arraystretch}{1}",
  "\\footnotesize{",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}ccccc@{}}",
  "\\toprule",
  " & \\multicolumn{2}{c}{\\textbf{Power}} & \\multicolumn{2}{c}{\\textbf{FDR}} \\\\",
  " \\cmidrule{2-3} \\cmidrule{4-5}",
  "\\textbf{Method} & \\textbf{Trait 1} & \\textbf{Trait 2} & \\textbf{Trait 1} & \\textbf{Trait 2} \\\\",
  "\\midrule",
  sapply(seq_len(nrow(rows)), function(r) tex_row(labels[r], cells[r, ], gap = FALSE)),
  "\\bottomrule",
  "\\end{tabular*}",
  "}",
  "\\end{table}"
)

dir.create("tables", showWarnings = FALSE)
write_tex(lines, "tables/Table_4.tex")
csv <- data.frame(
  rows,
  setNames(as.data.frame(power), paste0("power_", traits)),
  setNames(as.data.frame(fdr), paste0("fdr_", traits)),
  replications = length(reps)
)
write.csv(csv, "tables/Table_4.csv", row.names = FALSE)
print(csv, row.names = FALSE)
