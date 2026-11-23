# backends.R -- one model fit and prediction per imputation target (spec section 5.3)
#
# Purpose: `fit_predict()` fits ONE target (a numeric column, or the integer codes of a
#   factor) on the training rows with ranger, randomForest or mixOmics PLS, and predicts
#   the rows to impute. The stage engine calls it inside workers, each job with its own
#   pre-seeded L'Ecuyer stream already installed as the R RNG.
# Estimand: not an estimator. The prediction is the forest's / PLS regression's estimate
#   of the conditional mean (numeric target) or of the most likely class (factor target)
#   of the target given the other columns, for the rows to impute. `oob` estimates how
#   well the target can be predicted from the other columns (out-of-bag RMSE / scale, or
#   out-of-bag misclassification rate); it is a diagnostic and never feeds the imputation.
# Usage (internal):
#   r <- fit_predict("ranger", y, x_train, x_pred, is_factor, n_levels, scale,
#                    ntree, mtry, nComp)
#   r$pred   # unnamed numeric vector, or integer codes in 1..n_levels (factor target)
#   r$oob    # numeric scalar; NA_real_ for pls and for the no-predictor fallback
# Dependencies: ranger, randomForest (both Imports); mixOmics only for method "pls"
#   (Suggests, checked at call time). All calls namespace-qualified.
#
# RNG contract. All randomness comes from the current R RNG, which the caller has set;
# nothing here calls set.seed() or touches RNGkind(). ranger gets `seed = NULL`, so it
# draws its own seed from R's RNG (verified: the draw advances .Random.seed and the same
# state gives identical forests); ranger always runs with `num.threads = 1`, because
# parallelism comes from `future` and the thread count must not change a result.
# randomForest uses R's RNG directly; PLS is deterministic. Class prediction is an argmax
# with ties going to the first level (`max.col(ties.method = "first")`), which makes no
# draw; randomForest's own "response" prediction breaks ties at random, so its vote
# fractions (type = "prob") are used instead.
#
# Verified finding (spec section 5.1 open hypothesis; ranger 0.18.0): for a probability
# forest `fit$prediction.error` is a Brier-type score, mean((1 - p_true)^2) over rows,
# where p_true is the out-of-bag probability of the observed class. It equals that
# quantity to machine precision and differs from the out-of-bag misclassification rate
# (toy 3-class example: 0.0727 vs 0.0933). So the factor `oob` below is computed from
# `fit$predictions` (the out-of-bag class probabilities), never from `prediction.error`.
# For a regression forest `prediction.error` is the out-of-bag MSE, so
# sqrt(prediction.error) / scale is the scaled out-of-bag RMSE.

BACKENDS <- c("ranger", "randomForest", "pls")

# randomForest's regression warning for a response with five or fewer unique values.
# 0/1 and small-integer numeric columns are legitimate targets here, so this one message
# is muffled; any other warning still reaches the caller.
RF_FEW_VALUES_MSG <- "five or fewer unique values"

# mixOmics warns "At least one study has less than 5 samples, mean centering might not do
# as expected" when a PLS is fitted or predicted on 2-4 rows. Centring 2-4 rows is still
# ordinary centring and spec 5.1 allows such training sets, so this one advisory message
# is muffled; any other mixOmics warning still reaches the caller.
PLS_FEW_SAMPLES_MSG <- "less than 5 samples, mean centering"

# --- fit_predict ----------------------------------------------------------------

#' Fit one target and predict the rows to impute
#'
#' Order of operations per fit (each step is a named object a reader can print):
#' validate; drop zero-variance predictors; no predictor left -> mean / mode fallback;
#' cap `mtry` and `nComp` (PLS: at the training rank); factor target -> relabel to the observed levels; fit and
#' predict with the backend; map class scores back to codes in `1..n_levels`.
#'
#' Failure states: a non-finite prediction is an error (never "converged"); a missing
#' mixOmics for `"pls"` is an error; a factor target with a single observed level takes
#' the fallback (no model can be fitted on one class).
#'
#' @param method `"ranger"`, `"randomForest"` or `"pls"`.
#' @param y Numeric vector (numeric target) or integer codes in `1..n_levels` (factor).
#' @param x_train,x_pred Numeric matrices, same columns; callers subset with `drop = FALSE`.
#' @param is_factor Logical scalar.
#' @param n_levels Number of declared levels (factor target); ignored otherwise.
#' @param scale Stopping-rule scale of the target (numeric target); ignored otherwise.
#' @param ntree,mtry,nComp Forest size, predictors tried per split, PLS components.
#' @return List: `pred` (unnamed) and `oob` (numeric scalar).
#' @noRd
fit_predict <- function(method, y, x_train, x_pred, is_factor, n_levels, scale,
                        ntree, mtry, nComp) {
  # --- BLOCKING input checks: a wrong shape here would silently mis-align predictors.
  stopifnot(is.character(method), length(method) == 1L, method %in% BACKENDS,
            is.matrix(x_train), is.numeric(x_train), is.matrix(x_pred),
            is.numeric(x_pred), ncol(x_train) == ncol(x_pred),
            nrow(x_train) == length(y), nrow(x_train) >= 1L, nrow(x_pred) >= 1L,
            is.logical(is_factor), length(is_factor) == 1L, !is.na(is_factor),
            is.numeric(ntree), length(ntree) == 1L, ntree >= 1,
            is.numeric(mtry), length(mtry) == 1L, mtry >= 1,
            is.numeric(nComp), length(nComp) == 1L, nComp >= 1,
            all(is.finite(x_train)), all(is.finite(x_pred)))
  if (is_factor) {
    stopifnot(is.numeric(y), all(is.finite(y)), all(y == round(y)),
              is.numeric(n_levels), length(n_levels) == 1L, !is.na(n_levels),
              all(y >= 1 & y <= n_levels))
    y <- as.integer(y)
  } else {
    stopifnot(is.numeric(y), all(is.finite(y)),
              is.numeric(scale), length(scale) == 1L, is.finite(scale), scale > 0)
  }
  if (method == "pls" && !requireNamespace("mixOmics", quietly = TRUE))
    stop("method = \"pls\" needs the 'mixOmics' package, which is not installed")

  # Zero-variance predictors are dropped from the training rows' point of view, and the
  # same positions from x_pred: mixOmics warns on them and randomForest can stop.
  keep <- vapply(seq_len(ncol(x_train)), function(j) {
    r <- range(x_train[, j])
    r[2L] > r[1L]
  }, logical(1))
  x_train <- x_train[, keep, drop = FALSE]
  x_pred  <- x_pred[, keep, drop = FALSE]
  if (ncol(x_train) == 0L)
    return(fallback_prediction(y, nrow(x_pred), is_factor, n_levels))

  # Positional names: the backends never see user names, so duplicated, empty or odd
  # column names cannot break ranger/mixOmics name matching between fit and predict.
  colnames(x_train) <- colnames(x_pred) <- paste0("p", seq_len(ncol(x_train)))
  mtry  <- min(mtry, ncol(x_train))
  nComp <- min(nComp, ncol(x_train))
  if (method == "pls") {
    # PLS cannot extract more components than the rank of the CENTRED training
    # predictors: mixOmics 6.36 stops with "system is computationally singular" in
    # predict() otherwise. Rank-deficiency arises with two complementary indicator
    # columns of a 2-level factor, a duplicated predictor, or n_train <= ncol (rank is
    # at most n_train - 1 after centring). Capping at the rank loses no information:
    # components beyond it are numerically empty.
    rank <- qr(sweep(x_train, 2L, colMeans(x_train)))$rank
    # Rank 0 (predictors constant up to numerical tolerance): no model to fit, so the
    # training mean / mode is predicted, as when no predictor survives the variance drop.
    if (rank == 0L) return(fallback_prediction(y, nrow(x_pred), is_factor, n_levels))
    nComp <- max(1L, min(nComp, rank))
  }

  if (!is_factor) {
    res <- switch(method,
      ranger       = ranger_numeric(y, x_train, x_pred, scale, ntree, mtry),
      randomForest = rf_numeric(y, x_train, x_pred, scale, ntree, mtry),
      pls          = pls_numeric(y, x_train, x_pred, nComp))
    # BLOCKING: a non-finite prediction must never be treated as converged.
    if (!all(is.finite(res$pred))) stop("non-finite prediction")
    return(list(pred = as.numeric(unname(res$pred)), oob = res$oob))
  }

  # Factor target: fit on the levels observed in the training rows, then map back.
  observed <- sort(unique(y))
  if (length(observed) == 1L)
    return(fallback_prediction(y, nrow(x_pred), is_factor, n_levels))
  y_obs <- match(y, observed)        # 1..length(observed), the relabelled response
  res <- switch(method,
    ranger       = ranger_factor(y_obs, x_train, x_pred, ntree, mtry),
    randomForest = rf_factor(y_obs, x_train, x_pred, ntree, mtry),
    pls          = pls_factor(y_obs, x_train, x_pred, nComp))
  idx <- argmax_first(res$scores)
  list(pred = observed[idx], oob = res$oob)
}

# --- shared helpers ---------------------------------------------------------------

#' Training mean (numeric) or most frequent code (factor) when no model can be fitted
#'
#' `which.max(tabulate())` returns the first maximum, so ties go to the lowest code,
#' matching the spec's "ties: first level". `oob` is NA: nothing was predicted by a model.
#' @noRd
fallback_prediction <- function(y, n_pred, is_factor, n_levels) {
  if (is_factor) {
    value <- which.max(tabulate(y, nbins = n_levels))
  } else {
    value <- mean(y)
  }
  list(pred = rep(value, n_pred), oob = NA_real_)
}

#' Row-wise argmax of a class-score matrix; ties go to the first column, no RNG draw
#'
#' Errors on any non-finite score: a NaN score would otherwise become NA or a silently
#' arbitrary class.
#' @noRd
argmax_first <- function(scores) {
  stopifnot(is.matrix(scores))
  if (!all(is.finite(scores))) stop("non-finite prediction")
  max.col(scores, ties.method = "first")
}

#' Out-of-bag misclassification rate from out-of-bag class scores
#'
#' Rows that were never out of bag carry NaN scores (few trees) and are left out.
#' NA_real_ when no row has an out-of-bag score. DIAGNOSTIC value, never blocks.
#' @noRd
oob_misclass <- function(scores, y_obs) {
  ok <- is.finite(rowSums(scores))
  if (!any(ok)) return(NA_real_)
  mean(max.col(scores[ok, , drop = FALSE], ties.method = "first") != y_obs[ok])
}

# --- ranger -------------------------------------------------------------------------

ranger_numeric <- function(y, x_train, x_pred, scale, ntree, mtry) {
  fit <- ranger::ranger(x = x_train, y = y, num.trees = ntree, mtry = mtry,
                        num.threads = 1, seed = NULL, verbose = FALSE)
  pred <- stats::predict(fit, data = x_pred, num.threads = 1)$predictions
  # prediction.error of a regression forest is the out-of-bag MSE (see header).
  list(pred = pred, oob = finite_or_na(sqrt(fit$prediction.error) / scale))
}

ranger_factor <- function(y_obs, x_train, x_pred, ntree, mtry) {
  lev <- seq_len(max(y_obs))
  fit <- ranger::ranger(x = x_train, y = factor(y_obs, levels = lev),
                        num.trees = ntree, mtry = mtry, probability = TRUE,
                        num.threads = 1, seed = NULL, verbose = FALSE)
  scores <- stats::predict(fit, data = x_pred, num.threads = 1)$predictions
  # BLOCKING: score columns must follow the level order, or argmax picks the wrong class.
  stopifnot(identical(colnames(scores), as.character(lev)))
  # Misclassification from the stored out-of-bag probabilities, not prediction.error
  # (a Brier-type score for probability forests, see header).
  list(scores = scores, oob = oob_misclass(fit$predictions, y_obs))
}

# --- randomForest -------------------------------------------------------------------

# Only the "five or fewer unique values" warning is muffled; others pass through.
rf_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(RF_FEW_VALUES_MSG, conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

rf_numeric <- function(y, x_train, x_pred, scale, ntree, mtry) {
  fit <- rf_quiet(randomForest::randomForest(x = x_train, y = y, ntree = ntree,
                                             mtry = mtry))
  pred <- stats::predict(fit, newdata = x_pred)
  list(pred = pred, oob = finite_or_na(sqrt(fit$mse[ntree]) / scale))
}

rf_factor <- function(y_obs, x_train, x_pred, ntree, mtry) {
  lev <- seq_len(max(y_obs))
  fit <- randomForest::randomForest(x = x_train, y = factor(y_obs, levels = lev),
                                    ntree = ntree, mtry = mtry)
  # Vote fractions, argmax taken by argmax_first(): the class with the most votes, ties
  # to the first level, no RNG draw.
  scores <- stats::predict(fit, newdata = x_pred, type = "prob")
  # BLOCKING: score columns must follow the level order, or argmax picks the wrong class.
  stopifnot(identical(colnames(scores), as.character(lev)))
  list(scores = unname(scores),
       oob = finite_or_na(unname(fit$err.rate[ntree, "OOB"])))
}

# --- mixOmics PLS ---------------------------------------------------------------------

# mixOmics::predict() returns an array [rows, response columns, components]. Taking the
# slice for `nComp` components with drop = FALSE keeps the row dimension when only one
# row is predicted (otherwise mixOmics collapses it to a vector).
# Only the few-samples advisory is muffled (see PLS_FEW_SAMPLES_MSG).
pls_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(PLS_FEW_SAMPLES_MSG, conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

pls_scores <- function(fit, x_pred, nComp) {
  arr <- pls_quiet(stats::predict(fit, newdata = x_pred))$predict
  matrix(arr[, , nComp, drop = FALSE], nrow = dim(arr)[1L])
}

pls_numeric <- function(y, x_train, x_pred, nComp) {
  fit <- pls_quiet(mixOmics::pls(X = x_train,
                                 Y = matrix(y, ncol = 1L, dimnames = list(NULL, "y")),
                                 ncomp = nComp, mode = "regression", scale = TRUE))
  list(pred = pls_scores(fit, x_pred, nComp)[, 1L], oob = NA_real_)
}

pls_factor <- function(y_obs, x_train, x_pred, nComp) {
  K <- max(y_obs)
  # Indicator block of the observed levels (PLS-DA); every column has variance because
  # at least two levels are observed.
  Y <- outer(y_obs, seq_len(K), "==") + 0
  colnames(Y) <- paste0("c", seq_len(K))
  fit <- pls_quiet(mixOmics::pls(X = x_train, Y = Y, ncomp = nComp,
                                 mode = "regression", scale = TRUE))
  list(scores = pls_scores(fit, x_pred, nComp), oob = NA_real_)
}

#' NaN / Inf out-of-bag error (e.g. no row out of bag with very few trees) becomes NA
#' @noRd
finite_or_na <- function(v) if (is.finite(v)) as.numeric(v) else NA_real_
