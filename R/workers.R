# workers.R -- the future plan for one mvImpute() call (spec section 5.1)
#
# Purpose: decide how many workers a call runs on and put the matching plan in place
#   for the duration of the call, restoring the caller's plan on exit.
# Estimand: not an estimator. Results do not depend on the worker count: random numbers
#   come from pre-generated streams (rng.R), never from the worker that happens to run
#   a task.
# Usage (internal):
#   n <- local_workers(workers)    # in mvImpute(), after the no-missing early return;
#                                  # returns workers in effect
#   check_workers(workers)         # validation alone (mvImpute() runs it before the return)
# Dependencies: future (plan, multisession, nbrOfWorkers, availableCores).
#
# Rules (workers argument):
#   NULL, active plan not sequential -> use it, leave it in place (the user chose it)
#   NULL, plan sequential            -> multisession with k = min(8, availableCores() - 1)
#                                       workers if k >= 2, else stay sequential
#   1                                -> sequential for the call; an active worker pool is shut
#                                       down for it and restarted when the call exits
#   k > 1                            -> multisession with k workers for the call
# future cannot tell "never set" from "explicitly sequential", so users who want serial
# pass workers = 1. Multisession on every OS: forking is unavailable on Windows and
# unsafe in RStudio.
#
# Why k < 2 stays sequential: inside a future (mvImpute() called from a worker)
# availableCores() is 1, so k = 0. A one-worker multisession would only add start-up
# cost and nest a worker pool inside a worker.

MAX_AUTO_WORKERS <- 8L

# BLOCKING check, shared with mvImpute() so a bad `workers` is refused even when the call
# returns early and never starts a pool. The upper bound keeps as.integer() in
# local_workers() from turning e.g. 3e9 into NA with a warning.
check_workers <- function(workers) {
  if (!is.null(workers) && !(is_whole(workers, 1) && workers <= .Machine$integer.max))
    stop("workers must be NULL or a whole number >= 1", call. = FALSE)
  invisible(workers)
}

local_workers <- function(workers, envir = parent.frame()) {
  check_workers(workers)

  if (is.null(workers)) {
    # The user's plan is neither replaced nor restored: touching it would restart their workers.
    if (!inherits(future::plan(), "sequential")) return(as.integer(future::nbrOfWorkers()))
    k <- min(MAX_AUTO_WORKERS, as.integer(future::availableCores()) - 1L)
    if (k < 2L) return(1L)
    return(set_multisession(k, envir))
  }

  workers <- as.integer(workers)
  if (workers == 1L) {
    old <- future::plan(future::sequential)
    withr::defer(future::plan(old), envir)
    return(1L)
  }
  set_multisession(workers, envir)
}

# Install a multisession plan with k workers; the previous plan comes back when `envir`
# exits (also on error). future::plan() returns the previous strategy stack, which
# plan() accepts back.
set_multisession <- function(k, envir) {
  old <- future::plan(future::multisession, workers = k)
  withr::defer(future::plan(old), envir)
  k
}
