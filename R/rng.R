# rng.R -- seeds, RNG isolation and per-column random-number streams (spec section 5.5)
#
# Purpose: turn the user's `seed` into reproducible L'Ecuyer-CMRG streams, one per
#   (stage, column of W), without leaving any trace in the caller's RNG state.
# Estimand: not an estimator. Every function is a pure function of its arguments
#   except resolve_seed(NULL), which draws exactly one integer from the caller's RNG.
# Usage (internal):
#   seed    <- resolve_seed(seed)                        # NULL -> one caller draw
#   streams <- make_streams(seed, n_stage, n_col)        # streams[[stage]][[col]]
#   st      <- stream_for_iteration(streams[[s]][[j]], iter, crn)
#   local_rng_guard()    # at the top of any function that touches the RNG
# Caller order (mvImpute()): seed <- resolve_seed(seed); local_rng_guard(); local_workers(w).
#   resolve_seed() comes first so the one draw seed = NULL costs is made before the guard
#   records the state, i.e. survives the restore. local_workers() comes after the guard so
#   the plan teardown (a deferred call) runs inside the guarded window and the guard puts
#   .Random.seed back last. Never a plain on.exit() (without add = TRUE) in the same
#   function: it silently replaces the withr defers and drops these restores.
# Dependencies: parallel (nextRNGStream / nextRNGSubStream, the L'Ecuyer machinery);
#   withr (scoped restore of .Random.seed, including its absence).
#
# Why a guard is required, not optional. future.apply::future_lapply() advances the
# caller's .Random.seed even when it is given a list of pre-generated seeds (verified
# 2026-11-23, future.apply 1.20.2). Without a guard, calling mvImpute() would change
# the caller's later random numbers. local_rng_guard() restores .Random.seed and the
# three RNG kinds when the calling function exits, normally or by error.
#
# Why RNGkind() is restored beside withr::local_preserve_seed(). withr restores the
# caller's kind only when a .Random.seed exists. A caller with NO .Random.seed (fresh
# session, or after rm(.Random.seed)) gets the seed removed again, but the kind set
# inside the scope would stay: after make_streams() the session would silently run on
# L'Ecuyer-CMRG. The explicit RNGkind() restore closes exactly that case (tested in
# test-rng.R, "an absent .Random.seed stays absent").

# L'Ecuyer-CMRG is the only generator with the stream/sub-stream machinery that gives
# non-overlapping streams to parallel workers; the other two kinds are pinned so the
# streams do not depend on the caller's normal.kind or sample.kind.
STREAM_RNG_KIND        <- "L'Ecuyer-CMRG"
STREAM_RNG_NORMAL_KIND <- "Inversion"
STREAM_RNG_SAMPLE_KIND <- "Rejection"

# --- resolve_seed ---------------------------------------------------------------

# `seed` NULL: one draw from the caller's RNG, so set.seed(1); mvImpute(X) reproduces
# and this is the only touch of the caller's RNG state. Otherwise a single finite whole
# number within integer range, returned as integer. Failure state: anything else
# (NA, non-numeric, fractional, out of range, length != 1) is a BLOCKING error;
# as.integer() would otherwise turn 3e9 into NA with only a warning.
resolve_seed <- function(seed) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max, 1L))
  check_seed(seed)
  as.integer(seed)
}

check_seed <- function(seed) {
  ok <- is.numeric(seed) && length(seed) == 1L && is.finite(seed) &&
    seed == round(seed) && abs(seed) <= .Machine$integer.max
  if (!ok) stop("seed must be NULL or a single finite whole number", call. = FALSE)  # BLOCKING
  invisible(seed)
}

# --- local_rng_guard ------------------------------------------------------------

# Restore .Random.seed (including its absence) and the three RNG kinds when `envir`
# exits. Order matters: withr runs deferred calls last-in first-out, so the RNGkind()
# restore (registered second) runs FIRST; RNGkind() rescrambles the seed, and
# preserve_seed then puts the caller's seed back (or removes it again).
# suppressWarnings(): restoring sample.kind = "Rounding" warns, and that is the
# caller's own setting being put back, not a problem of ours.
local_rng_guard <- function(envir = parent.frame()) {
  kind_before <- RNGkind()
  withr::local_preserve_seed(.local_envir = envir)
  withr::defer(suppressWarnings(RNGkind(kind_before[1], kind_before[2], kind_before[3])),
               envir)
  invisible(NULL)
}

# --- make_streams ---------------------------------------------------------------

# Returns `n_stage` lists of `n_col` L'Ecuyer-CMRG states (integer vectors of length 7).
# `n_col` counts the columns of W (non-excluded): a stream is keyed by position in W, so
# adding or removing excluded columns cannot change any imputation (spec section 5.5).
# The first stream is the state set.seed(seed) produces; every further stream is
# nextRNGStream() of the previous one, in stage-major order (column fastest). Streams are
# therefore distinct. Stream (s, j) is chain position (s - 1) * n_col + j, so only stage-1
# streams are independent of n_col; later stages depend on it. That is harmless: n_col is
# ncol(W), fixed for a data set (excluded columns are not counted).
make_streams <- function(seed, n_stage, n_col) {
  check_seed(seed)
  stopifnot("n_stage must be a whole number >= 1" = is_whole(n_stage, 1),   # BLOCKING
            "n_col must be a whole number >= 1" = is_whole(n_col, 1))
  local_rng_guard()
  RNGkind(STREAM_RNG_KIND, STREAM_RNG_NORMAL_KIND, STREAM_RNG_SAMPLE_KIND)
  set.seed(as.integer(seed))
  current <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n_stage)
  for (i in seq_len(n_stage)) {
    out[[i]] <- vector("list", n_col)
    for (j in seq_len(n_col)) {
      out[[i]][[j]] <- current
      current <- parallel::nextRNGStream(current)
    }
  }
  out
}

# TRUE for a single finite whole number >= `min`. Shared by the argument checks here
# and in workers.R.
is_whole <- function(x, min) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) && x >= min
}

# --- stream_for_iteration -------------------------------------------------------

# crn = TRUE (common random numbers): every iteration of a stage reuses the stream, so a
# change between iterations reflects the data, not new random draws. crn = FALSE: iteration
# `iter` uses the (iter - 1)-th sub-stream, i.e. iteration 1 is the stream itself.
stream_for_iteration <- function(stream, iter, crn) {
  stopifnot("iter must be a whole number >= 1" = is_whole(iter, 1),        # BLOCKING
            "crn must be TRUE or FALSE" = is.logical(crn) && length(crn) == 1L && !is.na(crn))
  if (crn) return(stream)
  for (k in seq_len(iter - 1L)) stream <- parallel::nextRNGSubStream(stream)
  stream
}
