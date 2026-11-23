# latex_tables.R
#
# Small helpers shared by the make_tables.R scripts. The layout of each table
# is written out in the make_tables.R of its folder.

red <- function(x) paste0("\\textcolor{red}{", x, "}")
uline <- function(x) paste0("\\underline{", x, "}")

# Highlight cells of a character matrix where `flag` is TRUE
highlight <- function(cells, flag, how = red) {
  cells[flag & !is.na(flag)] <- how(cells[flag & !is.na(flag)])
  cells
}

# Column-wise maximum or minimum, all ties flagged
col_is_max <- function(x) {
  apply(x, 2, function(v) !is.na(v) & v == max(v, na.rm = TRUE))
}
col_is_min <- function(x) {
  apply(x, 2, function(v) !is.na(v) & v == min(v, na.rm = TRUE))
}

# One body row: "label &  & a & b & ... \\" (gap = TRUE leaves the empty
# second column used by most of the tables)
tex_row <- function(label, cells, gap = TRUE) {
  paste0(label, if (gap) " &  & " else " & ", paste(cells, collapse = " & "), " \\\\")
}

# Bold column headings such as \textbf{-0.75} & \textbf{-0.5} & ...
bold_cells <- function(x) paste0("\\textbf{", x, "}", collapse = " & ")

pct <- function(x) paste0(round(100 * x), "\\%")

write_tex <- function(lines, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, file)
  cat("written:", file, "\n")
}
