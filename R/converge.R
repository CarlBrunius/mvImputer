# converge.R -- stopping rule for the iterative imputation stages (spec section 5.4)
#
# Purpose: pure functions the stage engine calls once per iteration to decide which
#   imputed cells have stopped moving. No state, no RNG, no I/O.
# Estimand: not an estimator. `converge_cells()` returns, per cell, whether the
#   change between two successive predictions is small relative to the spread of the
#   feature (numeric) or whether the predicted level has been stable (factor-like).
# Usage (internal):
#   s   <- col_scales(W, miss, type)                    # once, before stage 1
#   res <- converge_cells(old, new, prev_same, s[col], is_factor[col], tol)
#   still_left <- !res$converged                        # exact complement
# Dependencies: stats (sd) only.
#
# Why the old rule was replaced (spec F2). The old engine used
#   abs(old - new) / old < tolerance
# which (i) passes a cell moving from -1 to -2 because the ratio is negative, (ii) gives
# Inf for any move away from 0 so 0/1 columns ran to the iteration cap, (iii) leaves a
# ratio exactly equal to the tolerance in neither "converged" nor "still left", and
# (iv) is not shift-invariant. Dividing by the SD of the observed values of the feature
# removes all four: the denominator is fixed, positive and location-free, and
# "converged" / "still left" are complements by construction.
#
# Why factor-like cells need two unchanged predictions. A cell whose first prediction
# merely equals the modal first guess would otherwise freeze on a model that was fitted
# to first-guess predictors. Requiring the level to be unchanged for two consecutive
# iterations means at least one prediction came from a model fitted on updated values.

# --- col_scales ---------------------------------------------------------------

#' Per-column scale for the stopping rule
#'
#' SD of the observed values of each numeric column, computed once before stage 1.
#' Returned unnamed, one value per column of `W`.
#'
#' Failure states: a numeric column with fewer than two observed values, or a constant
#' one, gets scale 0 (the engine removes constant columns beforehand, section 5.1, so
#' 0 is not expected in use; `converge_cells()` treats a 0 scale as "any change is
#' not converged" rather than dividing by it silently). Factor-like columns get `NA`
#' and are never scaled.
#'
#' @param W Numeric matrix (cells to impute hold any value; only observed ones are used).
#' @param miss Logical matrix, same dimensions as `W`, `TRUE` where the cell is missing.
#' @param type Character vector of length `ncol(W)`; `"numeric"` marks numeric columns,
#'   anything else is factor-like.
#' @return Unnamed numeric vector of length `ncol(W)`.
#' @noRd
col_scales <- function(W, miss, type) {
  # BLOCKING: shapes must agree, or scales silently attach to the wrong column.
  stopifnot(is.matrix(W), is.matrix(miss), identical(dim(W), dim(miss)),
            length(type) == ncol(W))
  out <- vapply(seq_len(ncol(W)), function(j) {
    if (!identical(type[[j]], "numeric")) return(NA_real_)
    x <- W[!miss[, j], j]
    # stats::sd() is NA for a single value; a single observation has no spread, so 0.
    if (length(x) < 2L) 0 else stats::sd(x)
  }, numeric(1))
  unname(out)
}

# --- converge_cells -----------------------------------------------------------

#' Decide which cells have converged
#'
#' Numeric cell: converged when `old == new` or `abs(new - old) / scale < tol`
#' (strict inequality: a change exactly at `tol` is not converged). Factor-like cell
#' (`is_factor`): converged when the prediction is unchanged now and was unchanged in
#' the previous iteration (`same & prev_same`); its `scale` is ignored (it is `NA`).
#'
#' @param old,new Numeric vectors: previous and current predictions (level codes for
#'   factor-like cells).
#' @param prev_same Logical vector: was the prediction unchanged in the previous
#'   iteration? Use `FALSE` in the first iteration.
#' @param scale Numeric vector from `col_scales()`, expanded to one value per cell;
#'   `NA` for factor-like cells.
#' @param is_factor Logical vector, `TRUE` for factor-like cells.
#' @param tol Scalar tolerance, a fraction of the feature SD.
#' @return List with `converged` (logical, length of `new`) and `same` (logical,
#'   `new == old`, to be passed back as `prev_same` next iteration).
#' @noRd
converge_cells <- function(old, new, prev_same, scale, is_factor, tol) {
  n <- length(new)
  # BLOCKING: equal lengths; recycling here would pair cells with the wrong scale.
  stopifnot(length(old) == n, length(prev_same) == n, length(scale) == n,
            length(is_factor) == n, length(tol) == 1L, is.finite(tol), tol > 0,
            !anyNA(is_factor), !anyNA(prev_same), !anyNA(old))
  # BLOCKING: a NaN/Inf prediction would make the comparison NA, and NA is neither
  # converged nor still left. Failing loudly is the only safe outcome.
  if (!all(is.finite(new))) {
    stop("non-finite prediction: ", sum(!is.finite(new)), " of ", n,
         " cells; the learner returned NA/NaN/Inf", call. = FALSE)
  }
  # BLOCKING: every numeric cell needs a usable scale; an NA scale would turn
  # `converged` into NA for that cell (factor cells never use their scale).
  if (anyNA(scale[!is_factor])) {
    stop("numeric cell with missing scale; col_scales() must give every numeric ",
         "column a finite scale", call. = FALSE)
  }

  same <- new == old
  converged <- logical(n)
  # Factor-like cells: two consecutive unchanged predictions.
  converged[is_factor] <- same[is_factor] & prev_same[is_factor]
  # Numeric cells: SD-scaled change. `same` is checked first so that a zero scale with
  # no change gives converged rather than 0/0 = NaN; a zero scale with a change gives
  # Inf, which is not < tol, i.e. not converged.
  num <- !is_factor
  change <- abs(new[num] - old[num]) / scale[num]
  converged[num] <- same[num] | (!is.nan(change) & change < tol)
  list(converged = converged, same = same)
}
