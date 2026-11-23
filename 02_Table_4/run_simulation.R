# Table 4: single-trait (JAM expanded) and multi-trait (flashfm) fine-mapping
# in the IL2RA region with OLS, St, DSh and GSR estimates. In each of N = 100
# replications, 5,000 of the 10,000 individuals are drawn, two traits are
# simulated, and both are fine-mapped with each estimator.
# Needs Java and data/IL2RA_population.rds (see README.md).
#
# Run from this folder:
#   Rscript run_simulation.R            all replications
#   Rscript run_simulation.R 1-20 45    replications 1 to 20 and 45

source("../R/All_eight_estimators.R")
source("../R/data_generation.R")
source("../R/simulation_tools.R")
source("table4_functions.R")

N <- 100
n <- 5000
population_file <- "data/IL2RA_population.rds"
results_dir <- "results"

args <- commandArgs(trailingOnly = TRUE)
reps <- 1:N
if (length(args)) {
  reps <- sort(unique(unlist(lapply(strsplit(args, "-", fixed = TRUE), function(r) {
    r <- as.integer(r)
    if (length(r) == 2) r[1]:r[2] else r
  }))))
  if (anyNA(reps) || any(reps < 1 | reps > N)) stop("replications must be between 1 and ", N)
}

if (!nzchar(Sys.which("java"))) stop("JAM needs Java on the PATH")
if (!file.exists(population_file)) stop(population_file, " is missing (see README.md)")

# flashfmZero and flashfm-savvySh.R of Jennifer Asimit, at fixed commits
flashfmZero_lib <- "flashfmZero_lib"
if (!nzchar(system.file(package = "flashfmZero", lib.loc = flashfmZero_lib))) {
  dir.create(flashfmZero_lib, showWarnings = FALSE)
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github(
    paste0("jennasimit/flashfmZero@", flashfmZero_commit),
    lib = flashfmZero_lib, upgrade = "never", repos = "https://cloud.r-project.org"
  )
}
flashfmZero_lib <- normalizePath(flashfmZero_lib)
.libPaths(c(flashfmZero_lib, .libPaths()))
if (!nzchar(system.file(package = "savvySh"))) install.packages("savvySh", repos = "https://cloud.r-project.org")

flashfm_file <- "flashfm-savvySh.R"
flashfm_url <- paste0(
  "https://raw.githubusercontent.com/jennasimit/flashfm-savvySh/",
  "ad7a101b296687dd5112f73ef42c9c0acece2783/flashfm-savvySh.R"
)
if (!file.exists(flashfm_file)) download.file(flashfm_url, flashfm_file, mode = "wb")
source(flashfm_file)
adapt_flashfmZero()

# The same code is loaded on every worker
code_files <- normalizePath(c(
  "../R/All_eight_estimators.R", "../R/data_generation.R", "table4_functions.R", flashfm_file
))
load_code <- function() {
  .libPaths(c(flashfmZero_lib, .libPaths()))
  for (f in code_files) sys.source(f, envir = globalenv())
  adapt_flashfmZero()
}

G <- read_population(population_file)
cat("population:", nrow(G), "individuals and", ncol(G), "variants with MAF > 0.005\n")

dir.create(results_dir, showWarnings = FALSE)
results_dir <- normalizePath(results_dir)

# Settings of the example in flashfm-savvySh.R. JAM takes its MCMC seed from
# R's generator, which is reset before each estimator.
fine_map <- function(d, method, k) {
  path <- file.path(tempdir(), sprintf("jam_rep_%03d_%s", k, method))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(path, recursive = TRUE))
  seed_rep(100000 + k)
  fm <- tryCatch(
    FLASHFMwithJAMd_rawX_selreg(
      selreg = estimators[[method]], Y = d$Y, X = d$X, TOdds = 1, cpp = 0.99,
      NCORES = 1, save.path = path, maxcv = 3, maxcv_stop = 10, jam.nM.iter = 2
    ),
    error = function(e) stop(method, ": ", conditionMessage(e), call. = FALSE)
  )
  fm$mpp.pp$MPP
}

one_rep <- function(k) {
  load_code()
  file <- file.path(results_dir, sprintf("rep_%03d.csv", k))
  error_file <- sub("\\.csv$", "_error.txt", file)
  seed_rep(123 + k)
  d <- draw_traits(G, n)
  out <- tryCatch(
    do.call(rbind, lapply(methods, function(m) {
      start <- Sys.time()
      mpp <- fine_map(d, m, k)
      seconds <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 1)
      data.frame(method = m, seconds = seconds, mpp_long(mpp, colnames(d$X), d$beta))
    })),
    error = function(e) conditionMessage(e)
  )
  if (is.character(out)) {
    writeLines(out, error_file)
    return(FALSE)
  }
  out <- data.frame(replication = k, seed = 123 + k, out)
  write.csv(out, paste0(file, ".tmp"), row.names = FALSE)
  file.rename(paste0(file, ".tmp"), file)
  if (file.exists(error_file)) file.remove(error_file)
  TRUE
}

todo <- Filter(function(k) !file.exists(file.path(results_dir, sprintf("rep_%03d.csv", k))), reps)
if (length(todo) == 0) {
  cat("these replications are all in", results_dir, "\n")
} else {
  cat("running", length(todo), "of the", length(reps), "replications, results in", results_dir, "\n")
  start <- Sys.time()
  start_workers(min(max_workers, length(todo)))
  ok <- unlist(future_lapply(todo, one_rep, future.seed = TRUE, future.chunk.size = 1))
  plan(sequential)
  cat("finished:", sum(ok), "replications saved in", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 1), "minutes")
  if (any(!ok)) cat(";", sum(!ok), "failed, see rep_*_error.txt")
  cat("\n")
}
