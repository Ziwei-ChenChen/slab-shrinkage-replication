# simulation_tools.R
#
# Running one simulation setting in parallel batches, with the results saved
# after every batch so an interrupted run carries on where it stopped.

library(future)
library(future.apply)

# Largest number of workers, 10 or fewer if the machine has fewer cores.
# Raise the 10 on a bigger machine. availableCores() counts only the cores
# given to the job on a cluster.
max_workers <- min(10, parallelly::availableCores(omit = 1))

start_workers <- function(n_workers = max_workers) {
  if (n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)
  invisible(n_workers)
}

# As many workers as fit in memory_gb when each needs per_worker GB
fit_workers <- function(per_worker, memory_gb) {
  w <- as.integer(max(1, min(max_workers, floor(memory_gb / per_worker))))
  if (nbrOfWorkers() != w) start_workers(w)
  w
}

read_case <- function(file) {
  if (!file.exists(file)) return(NULL)
  out <- read.csv(file, stringsAsFactors = FALSE)
  out[!duplicated(out$replication), , drop = FALSE]
}

# TRUE when all N replications of a setting are saved
case_done <- function(file, N) {
  done <- read_case(file)
  !is.null(done) && all(seq_len(N) %in% done$replication)
}

# one_rep(k) returns a one-row data frame with a column `replication`
run_case <- function(file, N, one_rep, batch_size = 50) {
  done <- read_case(file)
  todo <- setdiff(1:N, done$replication)
  if (length(todo) == 0) {
    cat("done already:", basename(file), "\n")
    return(invisible(done))
  }
  cat("running", basename(file), "-", length(todo), "replications left\n")
  for (reps in split(todo, ceiling(seq_along(todo) / batch_size))) {
    # The seed is set inside one_rep, so future.seed = TRUE only avoids a
    # warning from future. Each worker gets one chunk of the batch.
    out <- future_lapply(reps, one_rep, future.seed = TRUE)
    done <- rbind(done, do.call(rbind, out))
    done <- done[order(done$replication), , drop = FALSE]
    tmp <- paste0(file, ".tmp")
    write.csv(done, tmp, row.names = FALSE)
    file.rename(tmp, file)
  }
  invisible(done)
}
