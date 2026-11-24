# mvImpute.R -- the package's only export: multivariate imputation of a mixed table
#   (spec sections 5.1, 5.4, 5.5)
#
# Purpose: wire the internal steps together: validate, resolve the seed, encode, first
#   guess, one run_stage() per stage, decode, attach diagnostics.
# Estimand: not an estimator. At each missing cell the returned value is the model's
#   estimate of that value given the other (non-excluded) columns of the SAME table:
#   conditional mean (numeric) or most likely level (factor-like), at the point each
#   stage stopped. One deterministic completed table; no imputation uncertainty.
# Usage:
#   out <- mvImpute(X, seed = 1)          # X a data.frame or numeric matrix
#   attr(out, "mvImpute")                 # diagnostics list
# Dependencies: every internal file of the package (validate, encode, rng, workers, converge,
#   stage, backends). No new external dependency.
#
# Order of the calls below is load-bearing (see also R/rng.R):
#   check_input()   before encode_input(): checks must see the user's own types and values.
#   resolve_seed()  before local_rng_guard(): the ONE draw seed = NULL costs from the
#                   caller's RNG must survive the guard's restore.
#   local_rng_guard() before local_workers(): withr defers run last-in first-out, so the
#                   future plan is torn down INSIDE the guarded window and the guard puts
#                   .Random.seed back last.
#   local_workers() after the "nothing missing" return: no worker pool is started for nothing.
# No plain on.exit() (without add = TRUE) may appear in this function: it would silently
# drop the withr restores above.
#
# Verdict classes: stop() and stopifnot() are BLOCKING (the run cannot be trusted). The
# only DIAGNOSTIC outputs are message() (constant columns, progress, nothing to impute) and
# the final-stage not-converged warning(); the run still returns its table.
#
# The exclude boundary: excluded columns are read only by check_input() (type and name
# checks) and rows_without_observed(); encode_input() drops them before anything numeric
# is built, and decode_output() hands them back from a copy of X.

#' Impute missing values in a numeric or mixed-type table
#'
#' Iterative random-forest (or PLS) imputation. Each column with gaps is predicted from
#' all the other columns, the predictions replace the gaps, and this is repeated until the
#' imputed values stop moving. Numeric, integer, factor, ordered factor, character and
#' logical columns are supported. Factors are imputed as one target each.
#'
#' The algorithm runs in stages (`length(tol)` of them). Stage 1 is deliberately loose;
#' every later stage restarts from the previous result with all missing cells open again.
#' Within a stage, a numeric cell is frozen once its change between two iterations is
#' below `tol[stage]` times the standard deviation of that column's observed values; a
#' factor cell is frozen once its predicted level has been unchanged for two consecutive
#' iterations. A column's model is trained only on the rows where that column was
#' originally observed.
#'
#' Results are identical for a given `seed` whatever `workers` is (same machine, same R
#' and package versions). The caller's `.Random.seed` and `RNGkind()` are left as they
#' were, except that `seed = NULL` draws exactly one integer from the caller's generator.
#'
#' @section Assumptions:
#' * **Pooled imputation.** Every row, including rows that a later cross-validation fold
#'   holds out, informs every imputed value. Downstream cross-validated or test-set
#'   performance is therefore optimistic, more so if the outcome is a column of `X`. Leave
#'   the outcome out of `X`, or pass it in `exclude`: excluded columns are returned
#'   unchanged (gaps stay gaps), are never predictors and are never imputed.
#' * **Missing at random.** The models assume that values are missing at random given the
#'   other columns. Values missing because they are below a detection limit (missing not
#'   at random) are imputed too high: forests cannot predict below the observed range.
#' * **Single imputation.** The result is one deterministic completed table. Downstream
#'   standard errors that treat imputed values as observed are too small.
#' * **Linear dependencies.** Exact linear dependencies among numeric columns that you
#'   supply (your own one-hot columns, or parts and their total) can make the iterations
#'   oscillate. Recode a one-hot block as one factor column; factor columns are handled
#'   safely.
#'
#' @section Choosing `guess`:
#' `"minVar"` (the default) starts each numeric gap at the column minimum, which suits
#' omics intensities with detection-limit gaps. For a covariate table (clinical or
#' questionnaire variables) use `guess = "median"`. `"minTotal"` starts every numeric gap
#' at the smallest observed value over all numeric columns. A matrix or data frame of the
#' shape of `X` supplies a start per cell; it is read only at the missing cells of numeric
#' columns, where it must be a finite number. Factor-like columns always start at their
#' most frequent observed level (ties: the first level).
#'
#' @section Other notes:
#' * An integer matrix with gaps is returned as a double matrix (a matrix has one storage
#'   type). An integer column of a data frame comes back as double if it had gaps and as
#'   integer if it had none; imputed values are not rounded.
#' * Column names must be unique (duplicates are refused). Matrices without column names
#'   are reported as `V1`, `V2`, ... in messages and diagnostics, and returned without them.
#' * A column whose observed values are all equal is filled with that value (with a
#'   message) and is not a target.
#' * With `method = "pls"`, "most likely level" means the largest predicted indicator
#'   score, not a probability.
#' * `workers = 1` runs serially for the call; if a parallel plan is active it is shut
#'   down for the call and restarted when the call ends. `workers = NULL` uses an active
#'   non-sequential plan as it is and leaves it in place.
#' * `oob_error` is estimated on observed rows, so when the gaps are not missing at
#'   random it overstates how well the imputed values can be predicted.
#' * Known limitation: on tables with factor-like columns the final stage often reaches
#'   `maxIter`, because predicted factor levels can cycle between iterations, so the "not
#'   converged" warning is common there. The cells concerned are listed in the
#'   diagnostics' `not_converged`.
#'
#' @param X A data frame or a numeric matrix. The same class is returned.
#' @param method `"ranger"` (default), `"randomForest"`, or `"pls"` (needs the
#'   `mixOmics` package).
#' @param exclude Character vector of column names that are passed through untouched
#'   (identifiers, the outcome). The checks on column values and types skip them; the
#'   checks that column names are unique, that `X` is a data frame or numeric matrix and
#'   that every name in `exclude` exists cover all columns. At least two columns must
#'   remain.
#' @param seed `NULL` (draw one integer from the caller's random-number generator) or a
#'   single whole number within the integer range. The seed used is recorded in the
#'   diagnostics.
#' @param workers `NULL`: use the active [future::plan()] if it is not sequential,
#'   otherwise a multisession plan with `min(8, availableCores() - 1)` workers if that
#'   is at least 2, otherwise serial. `1`: serial. `k > 1`: a multisession plan with `k`
#'   workers. The previous plan is restored on exit.
#' @param tol Numeric vector, one value per stage: the convergence tolerance as a
#'   fraction of the standard deviation of the observed values of each numeric column.
#'   The defaults were chosen by simulation.
#' @param maxIter Whole numbers, one per stage: the iteration cap of that stage.
#'   `length(maxIter)` must equal `length(tol)`.
#' @param guess `"minVar"`, `"median"`, `"minTotal"`, or a matrix or data frame with the
#'   dimensions of `X` (see "Choosing `guess`").
#' @param forceZero `TRUE` floors numeric imputations at zero.
#' @param ntree Number of trees per forest.
#' @param mtry Number of predictors tried per split (a single whole number; capped at the
#'   number of predictors).
#' @param nComp Number of PLS components (capped at the rank of the training predictors).
#' @param verbose `TRUE` reports each stage and iteration through [message()].
#'
#' @return `X` with its missing values filled, in the same class, with the original
#'   column types (factor levels including unused ones, ordered factors, character,
#'   logical), row names and column names. Observed values are returned exactly as given.
#'   If `X` has no missing values in the columns to impute it is returned unchanged with
#'   a message; its diagnostics then have `workers = NA` and empty per-stage fields.
#'
#'   The diagnostics are in `attr(out, "mvImpute")`, a list with
#'   * `seed`: the seed used; `rng_kind`: `"L'Ecuyer-CMRG"`; `method`;
#'   * `settings`: list of the other arguments as given;
#'   * `workers`: number of workers used (integer);
#'   * `iterations`: iterations run per stage (integer vector);
#'   * `not_converged`: per stage, a data frame with columns `row` (row number in `X`)
#'     and `column` (column name) of the cells still iterating at the cap, zero rows if
#'     none;
#'   * `time`: elapsed seconds per stage;
#'   * `oob_error`: named numeric, one value per imputed column, from the last fit in
#'     the final stage. Numeric column: out-of-bag RMSE divided by the column's
#'     standard deviation. Factor-like column: out-of-bag misclassification rate. `NA`
#'     for `method = "pls"`;
#'   * `rows_no_observed`: row numbers with no observed value in the non-excluded
#'     columns (integer; such rows also give a warning).
#'
#' A `warning()` is given only if the final stage ends with unconverged cells; earlier
#' stages record theirs in `not_converged`.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(a = rnorm(40), b = rnorm(40))
#' d$c <- d$a + d$b + rnorm(40, sd = 0.1)
#' d$g <- factor(ifelse(d$a > 0, "up", "down"))
#' d$a[c(3, 9)] <- NA
#' d$g[c(5, 20)] <- NA
#' out <- mvImpute(d, seed = 1, workers = 1, ntree = 50, maxIter = c(3, 5))
#' attr(out, "mvImpute")$iterations
#' @export
mvImpute <- function(X, method = c("ranger", "randomForest", "pls"), exclude = NULL,
                     seed = NULL, workers = NULL, tol = c(0.05, 0.025),
                     maxIter = c(15, 60), guess = "minVar", forceZero = FALSE,
                     ntree = 100, mtry = 5, nComp = 2, verbose = FALSE) {
  method <- match.arg(method)
  # BLOCKING, and before check_input(): that function lets mtry = NULL through ("not
  # supplied") but the backends need a number, so mvImpute() has its own single message.
  if (!is_pos_whole(mtry)) stop("mtry must be a single positive whole number", call. = FALSE)
  flag_ok <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!flag_ok(forceZero) || !flag_ok(verbose))                                  # BLOCKING
    stop("forceZero and verbose must each be TRUE or FALSE", call. = FALSE)
  check_input(X, exclude = exclude, tol = tol, maxIter = maxIter, method = method,
              nComp = nComp, ntree = ntree, mtry = mtry, guess = guess)
  check_workers(workers)   # validates only; the pool starts later, after the early return

  seed <- resolve_seed(seed)          # before the guard: the one draw must survive it
  local_rng_guard()
  enc <- encode_input(X, exclude)

  settings <- list(exclude = exclude, tol = tol, maxIter = maxIter, guess = guess,
                   forceZero = forceZero, ntree = ntree, mtry = mtry, nComp = nComp,
                   verbose = verbose)
  diagnostics <- list(seed = seed, rng_kind = STREAM_RNG_KIND, method = method,
                      settings = settings, workers = NA_integer_, iterations = integer(),
                      not_converged = list(), time = numeric(), oob_error = numeric(),
                      rows_no_observed = rows_without_observed(X, exclude))
  if (!any(enc$miss)) {
    message("X has no missing values in the columns to impute; returned unchanged.")
    attr(X, "mvImpute") <- diagnostics
    return(X)
  }
  diagnostics$workers <- local_workers(workers)

  # Constant columns are filled by first_guess() and never iterate (R/stage.R); say so once.
  for (nm in enc$colnames[enc$colnames %in% enc$constant & colSums(enc$miss) > 0L])
    message("Column '", nm, "' is constant; its missing values were filled with that value.")

  W       <- first_guess(enc$W, enc, guess)
  scales  <- col_scales(enc$W, enc$miss, enc$type)   # SD of OBSERVED values, fixed for all stages
  n_stage <- length(tol)
  streams <- make_streams(seed, n_stage, ncol(W))

  diagnostics$iterations    <- integer(n_stage)
  diagnostics$not_converged <- vector("list", n_stage)
  diagnostics$time          <- numeric(n_stage)
  for (s in seq_len(n_stage)) {
    if (verbose) message("Stage ", s, " of ", n_stage, ": tol = ", tol[s],
                         ", maxIter = ", maxIter[s])
    t0  <- proc.time()[["elapsed"]]
    res <- run_stage(W, enc, scales, streams[[s]], tol[s], maxIter[s], method, ntree, mtry,
                     nComp, forceZero, verbose = verbose)
    diagnostics$time[s]            <- proc.time()[["elapsed"]] - t0
    diagnostics$iterations[s]      <- res$iterations
    diagnostics$not_converged[[s]] <- data.frame(row = as.integer(res$not_converged$row),
                                                 column = as.character(res$not_converged$column))
    W <- res$W                       # every stage starts from the previous stage's table
  }
  diagnostics$oob_error <- res$oob

  # BLOCKING: decode_output() passes NA through silently, which would hand back a table
  # that still has gaps without any sign.
  if (anyNA(W[enc$miss])) stop("internal error: missing cells remain after the last stage",
                               call. = FALSE)
  n_left <- nrow(diagnostics$not_converged[[n_stage]])
  if (n_left > 0L)                                                              # DIAGNOSTIC
    warning(n_left, " imputed cells not converged after maxIter in the final stage.",
            call. = FALSE)

  out <- decode_output(W, enc)
  attr(out, "mvImpute") <- diagnostics
  out
}
