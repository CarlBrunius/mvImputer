# encode.R -- the internal data model: encode a mixed table, build the design matrix,
#   make the first guess, decode back (spec sections 5.1, 5.2)
#
# Purpose: turn a data.frame or numeric matrix into one numeric matrix W (numeric columns
#   as values, factor-like columns as integer level codes), and turn W back into the
#   user's types. Everything the imputation stages do happens on W in between.
# Estimand: none. Nothing here estimates a population quantity. first_guess() only
#   produces starting values that the iteration replaces; encode/decode must be lossless.
# Usage (internal; check_input() has already run):
#   enc <- encode_input(X, exclude)         # list, see "Data model" below
#   W   <- first_guess(enc$W, enc, guess)   # every missing cell filled
#   dm  <- design_matrix(W, enc)            # predictors: list(D, map)
#   out <- decode_output(W, enc)            # same class and dimnames as X
# Dependencies: stats (median). R/validate.R for col_labels() and GUESS_CHOICES.
#
# Data model (`enc`, one entry per NON-EXCLUDED column, in original order, p of them):
#   W        n x p double. Numeric/integer: the values. Factor-like: level codes 1..K.
#            NA where missing (NaN is turned into NA).
#   miss     n x p logical, TRUE where the original value was missing.
#   type     "numeric" | "nominal" | "ordered" (R/converge.R keys on these labels).
#   class    original class: numeric, integer, factor, ordered, character, logical.
#   levels   named list (factor-like columns only) of the full level set, as character.
#   constant names of columns whose observed values are all equal.
#   exclude  excluded names.   X  the input, kept for decoding.
#   colnames internal names (V1..Vp for an unnamed matrix); had_names: logical.
#   keep     positions in X of the p non-excluded columns. Added to the plan's model:
#            first_guess() needs the position, not the name, to read a user guess table,
#            because names can repeat in a data.frame and are absent in a matrix.
#
# The exclude boundary: excluded columns are never read after the `keep` selection in
# encode_input(). They are not in W, miss or D, and decode_output() starts from a copy
# of X, so they come back untouched, gaps included.
#
# Names and row names: column names must be unique (check_input() refuses duplicates), so
# enc$constant and enc$levels are looked up by name. W, miss and D carry ROW names (X's own,
# else "1".."n"): a matrix with only column names returns NAMED scalars from W[i, "col"].
# Traps: W[rows, k] is a named vector (use unname() before identical()); which(arr.ind = TRUE)
# on miss returns row names in its dimnames.
#
# Observed values: decode_output() takes them from enc$W, never from the W it is given,
# so no stage can alter an observed cell by mistake; only missing cells come from W.
#
# Verdict classes: every stop() here is BLOCKING (a malformed W, enc or guess would
# give a wrong table without any other sign). There are no DIAGNOSTIC checks.

TYPE_NUMERIC <- "numeric"
TYPE_NOMINAL <- "nominal"
TYPE_ORDERED <- "ordered"
LOGICAL_LEVELS <- c("FALSE", "TRUE")   # fixed, so a logical column with only TRUE observed
                                       # still decodes from a two-level code space

# --- encode_input ---------------------------------------------------------------
encode_input <- function(X, exclude = NULL) {
  labels <- col_labels(X)
  keep   <- which(!(labels %in% exclude))
  cols   <- lapply(keep, function(j) encode_column(if (is.data.frame(X)) X[[j]] else X[, j]))
  names(cols) <- labels[keep]

  W <- do.call(cbind, lapply(cols, `[[`, "w"))
  # Row names are set (X's own, else "1".."n") because a matrix with only column names
  # returns NAMED scalars from W[i, "col"], which breaks exact comparisons downstream.
  rn <- rownames(X)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(X)))
  dimnames(W) <- list(rn, labels[keep])
  miss <- is.na(W)
  type <- vapply(cols, `[[`, "", "type")
  lvls <- lapply(cols[type != TYPE_NUMERIC], `[[`, "levels")
  # Constant = at least one observed value and all observed values equal. A column with no
  # observed value is not constant (check_input() has already refused those).
  is_const <- vapply(seq_along(cols), function(k) {
    length(unique(W[!miss[, k], k])) == 1L
  }, TRUE)

  list(W = W, miss = miss, type = unname(type),
       class = unname(vapply(cols, `[[`, "", "class")),
       levels = lvls, constant = labels[keep][is_const],
       exclude = if (is.null(exclude)) character(0) else exclude,
       X = X, colnames = labels[keep],
       had_names = !is.null(colnames(X)), keep = keep)
}

# One column -> list(w = numeric vector, class, type, levels). The level set is fixed here
# and nowhere else: logical is always FALSE/TRUE; character is the sorted observed values
# with method = "radix" (C-locale, so it does not depend on the user's locale); a factor
# keeps its declared levels, unused ones included.
encode_column <- function(x) {
  if (is.ordered(x))
    return(list(w = as.numeric(as.integer(x)), class = "ordered", type = TYPE_ORDERED,
                levels = levels(x)))
  if (is.factor(x))
    return(list(w = as.numeric(as.integer(x)), class = "factor", type = TYPE_NOMINAL,
                levels = levels(x)))
  if (is.character(x)) {
    lv <- sort(unique(x[!is.na(x)]), method = "radix")
    return(list(w = as.numeric(match(x, lv)), class = "character", type = TYPE_NOMINAL,
                levels = lv))
  }
  if (is.logical(x))
    return(list(w = as.numeric(x) + 1, class = "logical", type = TYPE_NOMINAL,
                levels = LOGICAL_LEVELS))
  if (is.numeric(x)) {
    w <- as.numeric(x)
    w[is.nan(w)] <- NA_real_
    return(list(w = w, class = if (is.integer(x)) "integer" else "numeric",
                type = TYPE_NUMERIC, levels = NULL))
  }
  stop("unsupported column type ", class(x)[1], " reached encode_input()",       # BLOCKING
       call. = FALSE)
}

# --- design_matrix --------------------------------------------------------------
# Predictor view of W. Nominal -> one 0/1 column per DECLARED level (so an unused level
# is an all-zero column; the backends decide what to do with it), ordered -> its integer
# code, numeric -> as is. A factor-like cell must already be filled: a hard 0/1 or code
# cannot be built from NA, and silently coding NA as "no level" would bias the models.
design_matrix <- function(W, enc) {
  if (!is.matrix(W) || !identical(dim(W), dim(enc$W)))
    stop("W must be a matrix with the dimensions of enc$W", call. = FALSE)         # BLOCKING
  blocks <- lapply(seq_len(ncol(W)), function(k) {
    nm <- enc$colnames[k]
    if (enc$type[k] == TYPE_NUMERIC) {
      m <- W[, k, drop = FALSE]
      colnames(m) <- nm
      return(m)
    }
    codes <- W[, k]
    lv <- enc$levels[[nm]]
    if (anyNA(codes))
      stop("column '", nm, "' has missing level codes; fill them (first_guess()) ",  # BLOCKING
           "before design_matrix()", call. = FALSE)
    # A code outside 1..K would give an all-zero one-hot row, i.e. "no level".
    if (any(codes != round(codes) | codes < 1 | codes > length(lv)))
      stop("column '", nm, "' holds values that are not valid level codes",       # BLOCKING
           call. = FALSE)
    if (enc$type[k] == TYPE_ORDERED) {
      m <- matrix(codes, ncol = 1L, dimnames = list(NULL, nm))
      return(m)
    }
    m <- outer(codes, seq_along(lv), "==") * 1
    colnames(m) <- paste0(nm, "_", lv)
    m
  })
  ends <- cumsum(vapply(blocks, ncol, 1L))
  starts <- ends - vapply(blocks, ncol, 1L) + 1L
  map <- Map(seq.int, starts, ends)
  names(map) <- enc$colnames
  D <- do.call(cbind, blocks)
  # D's names are internal labels and can collide (numeric column "f_a" next to factor f with
  # level "a"); map holds the indices, so unique names lose nothing.
  colnames(D) <- make.unique(colnames(D))
  list(D = D, map = map)
}

# --- first_guess ----------------------------------------------------------------
# Fills every missing cell of W. Observed values are read from enc$W and never changed.
#   constant column     -> its constant, whatever `guess` is
#   factor-like column  -> most frequent observed code; ties go to the LOWEST code
#                          (no RNG: the first guess must not consume or depend on a seed)
#   numeric column      -> minVar: column minimum | median: column median |
#                          minTotal: minimum over all numeric columns' observed values |
#                          a user table, read at this column's missing cells only
# A user table has X's full shape (it may carry excluded columns); the column of X is
# located through enc$keep. It is never read at observed cells, constant, factor-like
# or excluded columns, so NA or non-numeric content there is harmless.
first_guess <- function(W, enc, guess = "minVar") {
  if (!is.matrix(W) || !identical(dim(W), dim(enc$W)))
    stop("W must be a matrix with the dimensions of enc$W", call. = FALSE)         # BLOCKING
  is_table <- is.matrix(guess) || is.data.frame(guess)
  if (is_table && !identical(dim(guess), dim(enc$X)))
    stop("guess must have the same dimensions as X", call. = FALSE)                # BLOCKING
  if (!is_table && !(is.character(guess) && length(guess) == 1L && guess %in% GUESS_CHOICES))
    stop("guess must be one of ", paste0("'", GUESS_CHOICES, "'", collapse = ", "), # BLOCKING
         " or a matrix or data.frame with the dimensions of X", call. = FALSE)

  total_min <- if (identical(guess, "minTotal")) numeric_total_min(enc) else NA_real_
  const <- enc$colnames %in% enc$constant
  for (k in which(colSums(enc$miss) > 0L)) {
    miss_k <- enc$miss[, k]
    obs <- enc$W[!miss_k, k]
    if (!length(obs))
      stop("column '", enc$colnames[k], "' has no observed value to start from",   # BLOCKING
           call. = FALSE)
    W[miss_k, k] <- if (const[k]) obs[1L]
      else if (enc$type[k] != TYPE_NUMERIC) mode_code(obs)
      else if (is_table) table_guess(guess, enc, k)
      else switch(guess, minVar = min(obs), median = stats::median(obs), minTotal = total_min)
  }
  W
}

# Minimum of the observed values of every numeric-type column, constants included (a
# constant column's value is an observed value). NA when there is no numeric column.
numeric_total_min <- function(enc) {
  vals <- unlist(lapply(which(enc$type == TYPE_NUMERIC),
                        function(k) enc$W[!enc$miss[, k], k]))
  if (length(vals)) min(vals) else NA_real_
}

# Most frequent code; which.max() returns the first maximum, i.e. the lowest code on ties.
mode_code <- function(codes) {
  which.max(tabulate(codes, nbins = max(codes)))
}

# The user's guess for the missing cells of W column k. Errors name the column as the user
# sees it (V<k> for an unnamed matrix, via enc$colnames).
table_guess <- function(guess, enc, k) {
  j <- enc$keep[k]
  g <- if (is.data.frame(guess)) guess[[j]] else guess[, j]
  g <- g[enc$miss[, k]]
  if (!is.numeric(g))
    stop("guess is not numeric where column '", enc$colnames[k], "' is missing ",   # BLOCKING
         "(a guess table is read at the missing cells of numeric columns)", call. = FALSE)
  if (anyNA(g))
    stop("guess is NA at a missing cell of column '", enc$colnames[k], "'; a guess ", # BLOCKING
         "table must hold a value wherever X is missing in a numeric column",
         call. = FALSE)
  if (any(is.infinite(g)))
    stop("guess is Inf or -Inf at a missing cell of column '", enc$colnames[k], "'", # BLOCKING
         call. = FALSE)
  g
}

# --- decode_output --------------------------------------------------------------
# Back to the user's class, dimnames and column types. Starts from a copy of enc$X, so
# table-level attributes (dimnames, row names, class) and excluded columns are untouched;
# column attributes of DECODED columns are not kept, except on a double column. Per decoded column the value is
# enc$W at observed cells (exact copy of the input), W at missing cells, and the constant
# for a constant column. Cells of W that are still NA stay NA (a decode before
# first_guess(), used by the round-trip tests).
# Type rules (spec 5.1): an integer column with any gap returns double, unrounded; one
# without gaps stays integer. For an integer MATRIX the rule applies to the whole matrix
# (one storage type), so it is integer only if no cell was missing.
decode_output <- function(W, enc) {
  if (!is.matrix(W) || !identical(dim(W), dim(enc$W)))
    stop("W must be a matrix with the dimensions of enc$W", call. = FALSE)         # BLOCKING
  out <- enc$X
  was_integer_matrix <- is.matrix(out) && is.integer(out)
  if (was_integer_matrix) storage.mode(out) <- "double"
  const <- enc$colnames %in% enc$constant

  for (k in seq_len(ncol(W))) {
    miss_k <- enc$miss[, k]
    v <- unname(enc$W[, k])
    v[miss_k] <- W[miss_k, k]
    # Only missing cells: observed cells must come from enc$W alone (0 and -0 are one
    # constant but different bits).
    if (const[k]) v[miss_k] <- enc$W[!miss_k, k][1L]
    j <- enc$keep[k]
    orig <- if (is.data.frame(out)) out[[j]] else NULL
    col <- decode_column(v, enc$class[k], enc$levels[[enc$colnames[k]]], orig,
                         any_gap = any(miss_k), label = enc$colnames[k])
    if (is.data.frame(out)) out[[j]] <- col else out[, j] <- col
  }
  if (was_integer_matrix && !any(enc$miss)) storage.mode(out) <- "integer"
  out
}

# One decoded column. `orig` is the input column (data.frame only; NULL for a matrix).
decode_column <- function(v, cls, lv, orig, any_gap, label) {
  if (cls %in% c("factor", "ordered", "character", "logical")) {
    ok <- is.na(v) | (v == round(v) & v >= 1 & v <= length(lv))
    if (!all(ok))
      stop("column '", label, "' holds values that are not valid level codes",     # BLOCKING
           call. = FALSE)
    codes <- as.integer(v)
    return(switch(cls,
      factor    = , ordered = structure(codes, levels = lv, class = class(orig)),
      character = lv[codes],
      logical   = as.logical(lv[codes])))
  }
  if (cls == "integer" && !any_gap) return(as.integer(v))
  if (is.double(orig)) { orig[] <- v; return(orig) }
  v
}
