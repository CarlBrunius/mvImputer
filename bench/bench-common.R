# bench-common.R -- helpers shared by bench/calibrate_tol.R and bench/benchmark.R
#
# Purpose: build the simulated scenarios, run the old engine (StatTools::mvImpWrap) with
#   its iteration log read back, and score an imputed table against the simulated truth.
#   These were moved out of calibrate_tol.R unchanged (apart from the arguments named
#   below) so that both scripts score and run the old engine in exactly the same way.
# Estimand: none of its own; score_run() returns, for one imputed table, the RMSE of the
#   imputed numeric cells against the truth (raw, and per column scaled by the SD of the
#   column's observed values) and the share of imputed factor-like cells with the true label.
# Usage: source("bench/bench-common.R") from the package root, after
#   source("tests/testthat/helper-sim.R") (sim_numeric, sim_mixed, punch_holes).
# Dependencies: StatTools (the old engine, only inside old_engine() and
#   assert_old_engine()); withr (a scoped seed for the old engine's call); stats, utils.

HOLE_SEED_OFFSET <- 1000L                     # holes use a seed apart from the data's
HOLE_FRAC  <- 0.05

# mvImpWrap defaults written out; only the worker count varies between callers.
old_args <- function(workers)
  list(method = "RF", rfMeth = "ranger", tol1 = 0.05, n1 = 15, tol2 = 0.025, n2 = 60,
       nCore = workers)

# BLOCKING: the reference must be the pinned old engine, or the comparison is not the one
# the spec asks for.
assert_old_engine <- function() {
  stopifnot(requireNamespace("StatTools", quietly = TRUE),
            as.character(utils::packageVersion("StatTools")) == "0.0.917")
}

# --- scenarios ------------------------------------------------------------------------
# truth: complete table; holes: the table with gaps given to mvImpute(); old_input: the
# numeric matrix the old engine gets; decode_old: old engine output -> table like truth.
# p_numeric: number of columns of the numeric scenarios (the mixed one has a fixed layout).
make_scenario <- function(id, rep, p_numeric) {
  hole_seed <- HOLE_SEED_OFFSET + rep
  if (id %in% c("i", "ii")) {
    truth <- sim_numeric(200, p_numeric, seed = rep)
    holes <- punch_holes(truth, HOLE_FRAC, if (id == "i") "mcar" else "left", seed = hole_seed)
    return(list(truth = truth, holes = holes, old_input = holes, decode_old = identity))
  }
  truth <- sim_mixed(900, seed = rep, extra = 10)
  holes <- punch_holes(truth, HOLE_FRAC, "mcar", seed = hole_seed)
  onehot <- one_hot(holes)
  list(truth = truth, holes = holes, old_input = onehot$M,
       decode_old = function(M) decode_one_hot(M, onehot$blocks, holes))
}

# Factor-like columns (factor, ordered, character, logical) -> K indicator columns each;
# numeric columns pass through. model.frame(na.action = na.pass): the default drops NA rows.
one_hot <- function(d) {
  blocks <- list(); parts <- list()
  for (nm in names(d)) {
    x <- d[[nm]]
    if (is.numeric(x)) { parts[[nm]] <- matrix(as.numeric(x), ncol = 1, dimnames = list(NULL, nm)); next }
    f <- if (is.factor(x)) x else factor(x)
    m <- stats::model.matrix(~ f - 1, stats::model.frame(~ f, na.action = stats::na.pass))
    colnames(m) <- make.names(paste0(nm, "_", levels(f)))
    stopifnot(ncol(m) == nlevels(f), nrow(m) == nrow(d))          # BLOCKING: full K, no row dropped
    blocks[[nm]] <- list(cols = colnames(m), levels = levels(f))
    parts[[nm]] <- m
  }
  M <- do.call(cbind, unname(parts))
  attributes(M)[c("assign", "contrasts")] <- NULL
  list(M = M, blocks = blocks)
}

# Largest indicator wins; max.col(ties.method = "first") sends ties to the first level.
decode_one_hot <- function(M, blocks, template) {
  out <- template
  for (nm in names(template)) {
    if (is.null(blocks[[nm]])) { out[[nm]] <- M[, nm]; next }
    out[[nm]] <- blocks[[nm]]$levels[max.col(M[, blocks[[nm]]$cols, drop = FALSE],
                                             ties.method = "first")]
  }
  out
}

# --- scoring ----------------------------------------------------------------------------
# Scaled RMSE pooled over imputed cells of numeric columns with spread (a constant column
# is filled exactly and has scale 0); factor-like columns are compared as labels.
score_run <- function(imp, sc) {
  truth <- as.data.frame(sc$truth); holes <- as.data.frame(sc$holes); imp <- as.data.frame(imp)
  num_r <- num_s <- numeric(); hits <- logical()
  for (nm in names(truth)) {
    miss <- is.na(holes[[nm]])
    if (!any(miss)) next
    if (is.numeric(truth[[nm]])) {
      s <- stats::sd(holes[[nm]][!miss])
      if (!(s > 0)) next
      r <- as.numeric(imp[[nm]][miss]) - truth[[nm]][miss]
      num_r <- c(num_r, r); num_s <- c(num_s, r / s)
    } else {
      hits <- c(hits, as.character(imp[[nm]][miss]) == as.character(truth[[nm]][miss]))
    }
  }
  stopifnot(length(num_r) > 0, !anyNA(num_r), !anyNA(hits))     # BLOCKING: no cell left unimputed
  c(rmse_scaled = sqrt(mean(num_s^2)), rmse_raw = sqrt(mean(num_r^2)),
    factor_acc = if (length(hits)) mean(hits) else NA_real_)
}

# --- old engine -------------------------------------------------------------------------
# mvImpWrap on `input` with `workers` PSOCK workers, after a scoped set.seed(seed) (the
# brief asks for it; with_seed restores the session's RNG afterwards, and the kinds are
# pinned to R's defaults so the seed means the same whatever ran before). The seed has no
# effect on the engine's workers, which are not seeded, so its runs do not reproduce.
# Iterations and the cells left at the cap are read from the engine's printed log.
# Returns the raw engine output (decode it with the scenario), the iteration counts and
# the elapsed seconds of the engine call alone (no decoding, no scoring).
old_engine <- function(input, workers, seed) {
  args <- old_args(workers)
  t0 <- proc.time()[["elapsed"]]
  log <- utils::capture.output(
    res <- withr::with_seed(seed, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                            .rng_sample_kind = "Rejection",
                            do.call(StatTools::mvImpWrap, c(list(MAT = input), args))))
  secs <- proc.time()[["elapsed"]] - t0
  is_iter <- grepl("^Iteration", log)
  second  <- which(grepl("SECOND ROUND", log))
  stopifnot(length(second) == 1L, sum(is_iter) >= 2L)             # BLOCKING: log format as read
  it1 <- sum(is_iter[seq_len(second)]); it2 <- sum(is_iter) - it1
  # mvImp prints "Still N NAs" only when a stage stops at its cap with cells left, so the
  # stage-2 count is read from the lines after "SECOND ROUND" only (stage 1 has its own).
  log2  <- log[seq.int(second + 1L, length(log))]
  still <- regmatches(log2, regexpr("Still [0-9]+ NAs", log2))
  stopifnot(length(still) <= 1L)                                   # BLOCKING: log format as read
  left <- if (length(still)) as.numeric(gsub("[^0-9]", "", still)) else 0
  list(res = res, it1 = it1, it2 = it2, left = left, n2 = args$n2, seconds = secs)
}
