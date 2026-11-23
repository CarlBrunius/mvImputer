# stage.R -- the stage engine: iterate model fits until imputed cells stop moving
#   (spec sections 5.4 and 5.5)
#
# Purpose: one imputation stage. Each iteration refits every target column that still has
#   unfrozen missing cells, from the previous iteration's values, writes the predictions
#   into those cells, and freezes the cells the stopping rule calls converged.
# Estimand: not an estimator. The returned W holds, at missing cells, the model's estimate
#   of each missing value given the other columns (conditional mean for a numeric target,
#   most likely level for a factor target) at the point the stage stopped.
# Usage (internal; mvImpute() calls it once per stage s):
#   scales <- col_scales(enc$W, enc$miss, enc$type)      # once, before stage 1
#   streams <- make_streams(seed, n_stage, ncol(W))
#   res <- run_stage(W, enc, scales, streams[[s]], tol[s], maxIter[s], method,
#                    ntree, mtry, nComp, forceZero)
#   res$W; res$iterations; res$not_converged; res$oob
# RNG: run_stage() advances the caller's .Random.seed (future_lapply does, even with
#   pre-generated seeds), so its caller must hold local_rng_guard().
# Dependencies: future.apply (future_lapply with a list of per-job seeds: the only way to
#   give each job its own pre-generated stream whatever the worker count); R/encode.R,
#   R/converge.R, R/rng.R, R/backends.R.
#
# Two deliberate departures from mvImp:
# 1. Training rows (Carl, 2026-11-23). A column's model is trained only on the rows where
#    that column was ORIGINALLY observed. Old mvImp trained on every row not still
#    iterating, so its own converged imputations became training targets for it.
# 2. A factor is fitted as one target (spec F3). Its own indicator (or code) columns are
#    removed from its predictors through dm$map, so it cannot predict itself from its
#    current guess; that self-reference made old mvImp oscillate.
#
# Seeding (spec 5.5). Job k runs under stream_for_iteration(streams[[col]], iter, crn),
# keyed by the job's W column, never by its position in the job list: as columns freeze
# and drop out of the list, no other column's random numbers shift. Jobs are fitted from
# one design matrix built before any of them runs and written back only after all have
# returned, so the result does not depend on job order or the number of workers.
#
# What crosses to the workers: stage_fit_job() is a namespace-level function and gets
# the design matrix D, its job and the backend settings as explicit arguments. `enc` is
# never an argument or a global of the worker, so enc$X (excluded columns, possibly an
# outcome) cannot reach a model.

# --- stage_jobs -------------------------------------------------------------------

# One job per target column with at least one unfrozen missing cell. `train` is the
# originally observed rows whatever is frozen (departure 1 above); `pred` the unfrozen
# missing rows. unname(): miss carries row names, and which() would carry them along.
stage_jobs <- function(miss, frozen, targets) {
  stopifnot(is.matrix(miss), is.logical(miss), is.matrix(frozen), is.logical(frozen),  # BLOCKING
            identical(dim(miss), dim(frozen)), !anyNA(miss), !anyNA(frozen),
            is.numeric(targets), all(targets %in% seq_len(ncol(miss))))
  jobs <- lapply(as.integer(targets), function(col) {
    pred <- unname(which(miss[, col] & !frozen[, col]))
    if (!length(pred)) return(NULL)
    list(col = col, train = unname(which(!miss[, col])), pred = pred)
  })
  jobs[!vapply(jobs, is.null, TRUE)]
}

# --- stage_fit_job (runs on a worker) ---------------------------------------------

# `job` carries col/train/pred plus what fit_predict() needs about the target: `y` (its
# originally observed values), `drop` (its own columns in D), `is_factor`, `n_levels`,
# `scale`. Everything else the worker sees is in the argument list.
stage_fit_job <- function(job, D, method, ntree, mtry, nComp) {
  x_train <- D[job$train, -job$drop, drop = FALSE]
  x_pred  <- D[job$pred, -job$drop, drop = FALSE]
  fit_predict(method, job$y, x_train, x_pred, job$is_factor, job$n_levels, job$scale,
              ntree, mtry, nComp)
}

# --- run_stage ----------------------------------------------------------------------

run_stage <- function(W, enc, scales, streams, tol, maxIter, method, ntree, mtry, nComp,
                      forceZero, crn = TRUE, verbose = FALSE) {
  miss <- enc$miss
  # BLOCKING: a misaligned W, scale or stream vector would pair a column with another
  # column's scale or random numbers, silently.
  stopifnot(is.matrix(W), identical(dim(W), dim(miss)), !anyNA(W),
            length(scales) == ncol(W), is.list(streams), length(streams) == ncol(W),
            is.numeric(tol), length(tol) == 1L, is.finite(tol), tol > 0,
            is_whole(maxIter, 1),
            is.logical(forceZero), length(forceZero) == 1L, !is.na(forceZero),
            is.logical(verbose), length(verbose) == 1L, !is.na(verbose))

  is_factor <- enc$type != TYPE_NUMERIC
  constant  <- enc$colnames %in% enc$constant
  targets   <- which(colSums(miss) > 0L & !constant)
  # Constant columns' cells start frozen: first_guess() has filled them with the
  # constant, and their scale is 0, so the stopping rule must never see them.
  frozen <- matrix(FALSE, nrow(W), ncol(W))
  frozen[, constant] <- TRUE
  # Spec 5.4: every stage starts with all (non-constant) missing cells unfrozen and no
  # history; a factor cell needs two unchanged predictions within THIS stage.
  prev_same <- matrix(FALSE, nrow(W), ncol(W))
  oob <- stats::setNames(rep(NA_real_, length(targets)), enc$colnames[targets])

  iterations <- 0L
  for (iter in seq_len(maxIter)) {
    jobs <- stage_jobs(miss, frozen, targets)
    if (!length(jobs)) break
    dm <- design_matrix(W, enc)
    payload <- lapply(jobs, function(job) {
      col <- job$col
      c(job, list(y = unname(W[job$train, col]), drop = dm$map[[col]],
                  is_factor = is_factor[col],
                  n_levels = if (is_factor[col]) length(enc$levels[[enc$colnames[col]]]) else NA,
                  scale = scales[col]))
    })
    seeds <- lapply(jobs, function(job) stream_for_iteration(streams[[job$col]], iter, crn))
    fits <- future.apply::future_lapply(payload, stage_fit_job, D = dm$D, method = method,
                                        ntree = ntree, mtry = mtry, nComp = nComp,
                                        future.seed = seeds)

    # Write back only after every job has returned (see header: order-free result).
    for (k in seq_along(jobs)) {
      col <- jobs[[k]]$col; rows <- jobs[[k]]$pred
      new <- fits[[k]]$pred
      if (forceZero && !is_factor[col]) new <- pmax(new, 0)
      res <- converge_cells(old = unname(W[rows, col]), new = new,
                            prev_same = prev_same[rows, col],
                            scale = rep(scales[col], length(rows)),
                            is_factor = rep(is_factor[col], length(rows)), tol = tol)
      W[rows, col] <- new
      frozen[rows, col] <- res$converged
      prev_same[rows, col] <- res$same
      oob[[enc$colnames[col]]] <- fits[[k]]$oob
    }
    iterations <- iter
    if (verbose) message("Iteration ", iter, ": ", sum(miss & !frozen), " cells still iterating")
  }

  left <- which(miss & !frozen, arr.ind = TRUE)
  list(W = W, iterations = iterations,
       not_converged = data.frame(row = unname(left[, "row"]),
                                  column = enc$colnames[left[, "col"]]),
       oob = oob)
}
