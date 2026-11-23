# validate.R -- input checks that run before encoding (spec section 5.1)
#
# Purpose: make bad input fail loudly instead of producing a plausible-looking but wrong
#   imputation. check_input() raises the errors and warnings of spec section 5.1;
#   rows_without_observed() is shared with mvImpute()'s diagnostics.
# Estimand: none. Nothing here estimates anything; it only validates arguments.
# Usage (internal):
#   check_input(X, exclude, tol, maxIter, method, nComp, ntree, mtry, guess)  # NULL, or stop/warn
#   rows_without_observed(X, exclude)                  # integer row numbers
# Dependencies: stats (complete.cases). mixOmics is only probed, via has_pkg().
#
# Scope. Every check except "unknown column" looks at NON-EXCLUDED columns only, so an
# identifier or an outcome with gaps can be passed through with `exclude`. All checks run
# on the raw table, before any encoding: they must see the user's own types and values.
#
# Verdict classes: every stop() below is BLOCKING (the run cannot be trusted or cannot
# start); the two warning() calls are DIAGNOSTIC (the run proceeds; the user is told
# about a property that makes the result weaker than it looks).
#
# Not checked here, on purpose: `method` is matched by mvImpute() with match.arg() before
# this runs; `seed` is validated by resolve_seed() (R/rng.R).

GUESS_CHOICES      <- c("minVar", "median", "minTotal")
IDENTIFIER_MIN_OBS <- 3L   # below this, "all distinct" is chance (2 observed values)
NAMES_SHOWN        <- 5L   # columns / rows listed in a message before "and N more"

# Wrapper so a test can mock "mixOmics is not installed" (local_mocked_bindings).
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# Column labels used in every message; V1, V2, ... for a matrix without colnames, so
# messages and the output's (absent) dimnames agree.
col_labels <- function(X) {
  nm <- colnames(X)
  if (is.null(nm)) paste0("V", seq_len(ncol(X))) else nm
}

fmt_names <- function(x, max = NAMES_SHOWN) {
  shown <- paste0("'", x[seq_len(min(length(x), max))], "'", collapse = ", ")
  if (length(x) > max) paste0(shown, " and ", length(x) - max, " more") else shown
}

# --- rows_without_observed ------------------------------------------------------
# Row numbers where every non-excluded column is NA (NaN counts as NA, as is.na() does).
# Used by check_input() for the warning and by mvImpute() for the diagnostics list, so
# both always agree on what "no observed value" means.
rows_without_observed <- function(X, exclude = NULL) {
  keep <- !(col_labels(X) %in% exclude)
  observed <- !is.na(X)                     # logical matrix for matrix and data.frame
  unname(which(rowSums(observed[, keep, drop = FALSE]) == 0))
}

# --- check_input ----------------------------------------------------------------
# Order matters: structural errors first (so later checks can assume a table of at
# least 2 supported columns), argument errors next, warnings last (an error pre-empts
# a warning). `mtry` NULL means "not supplied"; the others are always checked, because
# a malformed nComp that happens to be unused by this method is still a malformed call.
check_input <- function(X, exclude, tol, maxIter, method, nComp,
                        ntree = 100L, mtry = NULL, guess = "minVar") {
  if (!(is.data.frame(X) || is.matrix(X)))
    stop("X must be a data.frame or matrix, not ", class(X)[1], call. = FALSE)  # BLOCKING
  if (is.matrix(X) && !is.numeric(X))
    stop("A matrix X must be a numeric matrix (it is ", typeof(X),              # BLOCKING
         "); use a data.frame for mixed types", call. = FALSE)

  # Names are the key of the internal data model (levels, constants, exclude); duplicates
  # would silently mix up columns. Excluded columns count: exclude by name is ambiguous.
  nm <- colnames(X)
  if (!is.null(nm) && anyDuplicated(nm))
    stop("Column names must be unique (duplicated: ", fmt_names(unique(nm[duplicated(nm)])), # BLOCKING
         "); rename them, e.g. with make.unique().", call. = FALSE)

  labels <- col_labels(X)
  check_exclude(exclude, colnames(X))   # real names only: V1.. are labels, not names
  use <- which(!(labels %in% exclude))
  if (length(use) < 2L)
    stop("Need at least 2 columns to impute from each other after exclude (got ",  # BLOCKING
         length(use), ")", call. = FALSE)

  cols <- lapply(use, function(j) if (is.data.frame(X)) X[[j]] else X[, j])
  names(cols) <- labels[use]
  check_columns(cols)
  check_tuning(tol, maxIter, method, nComp, ntree, mtry)
  check_guess(guess, X)

  rows <- rows_without_observed(X, exclude)
  if (length(rows))                                                             # DIAGNOSTIC
    warning(length(rows), " row(s) have no observed value in the columns to impute (",
            "row numbers ", paste(rows[seq_len(min(length(rows), NAMES_SHOWN))], collapse = ", "),
            if (length(rows) > NAMES_SHOWN) ", ...", "); they are imputed from first guesses and ",
            "other imputations only", call. = FALSE)
  for (block in find_one_hot_blocks(cols))                                      # DIAGNOSTIC
    warning("Columns ", fmt_names(block, 20L), " are numeric 0/1 columns that sum to 1 ",
            "in every complete row (a one-hot block). Imputed separately they can give ",
            "rows with no or several 1s; recode them as one factor column",
            call. = FALSE)
  invisible(NULL)
}

check_exclude <- function(exclude, real_names) {
  if (is.null(exclude)) return(invisible(NULL))
  # Not coerced: exclude = 1 would otherwise be read as the name "1", or ignored.
  if (!is.character(exclude) || anyNA(exclude))
    stop("exclude must be a character vector of column names", call. = FALSE)   # BLOCKING
  unknown <- setdiff(exclude, real_names)
  if (length(unknown))
    stop("exclude names an unknown column: ", fmt_names(unknown),               # BLOCKING
         if (is.null(real_names)) " (X has no column names)"
         else paste0(" (columns are ", fmt_names(real_names), ")"), call. = FALSE)
  invisible(NULL)
}

# `cols`: named list of the non-excluded columns. Each check is one pass over the list,
# and a check only runs once the earlier ones hold (e.g. Inf is looked for only in
# columns already known to be a supported type).
check_columns <- function(cols) {
  supported <- vapply(cols, function(x) is.numeric(x) || is.factor(x) ||
                        is.character(x) || is.logical(x), TRUE)
  if (!all(supported))
    stop("unsupported column type in ", fmt_names(names(cols)[!supported]), " (",    # BLOCKING
         paste(vapply(cols[!supported], function(x) class(x)[1], ""), collapse = ", "),
         "); supported are numeric, integer, factor, ordered, character and logical",
         call. = FALSE)

  has_inf <- vapply(cols, function(x) is.numeric(x) && any(is.infinite(x)), TRUE)
  if (any(has_inf))
    stop("Inf or -Inf values in ", fmt_names(names(cols)[has_inf]),                  # BLOCKING
         ": set them to NA if they mean missing (they are not treated as missing)",
         call. = FALSE)

  n_obs   <- vapply(cols, function(x) sum(!is.na(x)), 1L)
  too_few <- n_obs < 2L & n_obs < lengths(cols)
  if (any(too_few))
    stop("fewer than 2 observed values in ", fmt_names(names(cols)[too_few]),        # BLOCKING
         " (a column with gaps needs at least 2 observed values to be imputed)",
         call. = FALSE)

  # All-distinct text would become one class per row and be "imputed" as noise.
  is_id <- vapply(cols, function(x) {
    if (!(is.character(x) || is.factor(x))) return(FALSE)
    obs <- x[!is.na(x)]
    length(obs) >= IDENTIFIER_MIN_OBS && !anyDuplicated(obs)
  }, TRUE)
  if (any(is_id))
    stop("Column(s) ", fmt_names(names(cols)[is_id]), " look like identifiers (every ", # BLOCKING
         "observed value is distinct); remove them from the table or pass them in exclude",
         call. = FALSE)
  invisible(NULL)
}

is_pos_whole <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == round(x)
}

check_tuning <- function(tol, maxIter, method, nComp, ntree, mtry) {
  if (length(tol) != length(maxIter) || length(tol) < 1L)
    stop("tol and maxIter must give one value per stage: equal length, at least 1 ",  # BLOCKING
         "(got lengths ", length(tol), " and ", length(maxIter), ")", call. = FALSE)
  # Non-finite values are refused too: NA would crash the comparison later, Inf would
  # make a stage "converge" (tol) or never stop (maxIter) without saying so.
  if (!is.numeric(tol) || any(!is.finite(tol)) || any(tol <= 0))
    stop("tol must be finite numbers > 0 (a fraction of each feature's SD)",         # BLOCKING
         call. = FALSE)
  if (!is.numeric(maxIter) || any(!is.finite(maxIter)) || any(maxIter < 1) ||
      any(maxIter != round(maxIter)))
    stop("maxIter must be whole numbers >= 1", call. = FALSE)                        # BLOCKING
  if (identical(method, "pls") && !has_pkg("mixOmics"))
    stop("method = \"pls\" needs the 'mixOmics' package, which is not installed",    # BLOCKING
         call. = FALSE)
  if (!is_pos_whole(nComp))
    stop("nComp must be a single positive whole number", call. = FALSE)              # BLOCKING
  if (!is_pos_whole(ntree))
    stop("ntree must be a single positive whole number", call. = FALSE)              # BLOCKING
  if (!is.null(mtry) && !is_pos_whole(mtry))
    stop("mtry must be NULL or a single positive whole number", call. = FALSE)       # BLOCKING
  invisible(NULL)
}

# Only the shape is checked here. Non-NA at missing cells is checked by first_guess() (R/encode.R).
check_guess <- function(guess, X) {
  ok_name  <- is.character(guess) && length(guess) == 1L && guess %in% GUESS_CHOICES
  ok_table <- (is.matrix(guess) || is.data.frame(guess)) && identical(dim(guess), dim(X))
  if (!(ok_name || ok_table))
    stop("guess must be one of ", fmt_names(GUESS_CHOICES), " or a matrix or ",      # BLOCKING
         "data.frame with the same dimensions as X (", paste(dim(X), collapse = " x "),
         ")", call. = FALSE)
  invisible(NULL)
}

# --- one-hot detection ----------------------------------------------------------
# Brief's approach: among numeric columns whose observed values are exactly {0, 1}, two
# columns are linked when never both 1 in rows where both are observed; within each
# connected group of >= 2 columns, report the group if its row sum is 1 in every row
# complete for the group (and at least one such row exists).
# Design call beyond the brief: a column with only one distinct observed value (all 0
# or all 1) is not a candidate. It is "never both 1" with everything, so an all-zero
# column would link two unrelated one-hot blocks into one group whose sums are 2.
# Known limit: an unrelated rare 0/1 indicator that never co-occurs with a block's
# columns joins its group and hides the block (false negative; the warning is advisory).
find_one_hot_blocks <- function(cols) {
  is_bin <- vapply(cols, function(x) {
    obs <- x[!is.na(x)]
    is.numeric(x) && length(obs) > 0L && all(obs %in% c(0, 1)) && length(unique(obs)) == 2L
  }, TRUE)
  if (sum(is_bin) < 2L) return(list())
  B <- do.call(cbind, cols[is_bin])
  is_one <- matrix(as.numeric(!is.na(B) & B == 1), nrow(B))
  linked <- crossprod(is_one) == 0               # TRUE: never both 1
  diag(linked) <- FALSE
  grp <- connected_groups(linked)

  blocks <- list()
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) < 2L) next
    sub  <- B[, idx, drop = FALSE]
    done <- stats::complete.cases(sub)
    if (any(done) && all(rowSums(sub[done, , drop = FALSE]) == 1))
      blocks[[length(blocks) + 1L]] <- colnames(B)[idx]
  }
  blocks
}

# Connected components of a symmetric logical adjacency matrix: integer group labels.
connected_groups <- function(linked) {
  k <- nrow(linked)
  grp <- integer(k)
  g <- 0L
  for (i in seq_len(k)) {
    if (grp[i] != 0L) next
    g <- g + 1L
    grp[i] <- g
    todo <- i
    while (length(todo)) {
      cur  <- todo[1L]
      todo <- todo[-1L]
      nb   <- which(linked[cur, ] & grp == 0L)
      grp[nb] <- g
      todo <- c(todo, nb)
    }
  }
  grp
}
