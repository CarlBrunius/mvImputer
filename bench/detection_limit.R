# detection_limit.R -- imputing values that went missing below a detection limit (simulated)
#
# Purpose: compare mvImpute() settings (ranger or PLS, log or original scale, forceZero, number
#   of PLS components, guess) with a column-minimum fill, when low values are more likely to be
#   missing than high ones. Writes bench/output/detection-limit.md, the source of every figure
#   quoted in the README section "Choosing settings for your data".
# Estimand: per data set (seed) and method, how far the imputed cells are from the simulated
#   truth: RMSE with each error divided by the SD of the column's observed values, on the log
#   scale and on the original scale; the share of imputed cells above the truth, below the
#   column's observed minimum, and zero or negative. Reported as the mean (and SD) over 5 data
#   sets, plus in how many of the 5 a method beat the minimum fill. No uncertainty beyond that.
# Usage (from the package root, after R CMD INSTALL --no-docs .; about 10 minutes, one process):
#   Rscript bench/detection_limit.R          # every seed not yet cached, then the report
#   Rscript bench/detection_limit.R report   # the report again from the cached runs
#   The published run was started at BelowNormal priority (PowerShell Start-Process
#   -PriorityClass BelowNormal), so the machine stays usable; children inherit it.
# Dependencies: mvImputer (installed), mixOmics (the PLS backend), withr (via helper-sim.R).
#   Helpers from tests/testthat/helper-sim.R (sim_numeric, punch_holes, sim_seeded) and
#   bench/bench-common.R (HOLE_SEED_OFFSET, HOLE_FRAC; its StatTools helpers are not called).
#
# Design (fixed; two earlier exploratory scripts were merged into this one):
# - Data: sim_numeric(200, 100, seed = s), s = 1..5: log-normal, rank-3 signal plus noise.
# - Missingness, 5 % of cells, three patterns, each with its own hole seed:
#     left   exactly the lowest 5 % of every column (hole seed 1000 + s; no randomness).
#     graded a soft limit: in each column 5 % of rows are drawn with weight exp(-2 z), z the
#            standardised log value, so low values are likely but not certain to go missing
#            (hole seed 3000 + s).
#     mixed  the lowest 2.5 % of every column plus 2.5 % of the remaining cells at random
#            (seeds 1000 + s and 2000 + s).
# - mvImpute(..., seed = s, workers = 1): serial, so seconds contain no worker pool.
# - Original-scale and forceZero fits are run for the "left" pattern only. Their log-scale
#   counterparts (ranger, PLS with 2 components) are the "left" rows of the log-scale runs
#   and are not fitted twice (the call is identical, so is the result).
# - The minimum fill sets every gap to its column's smallest observed value.
# BLOCKING checks: positive truth, observed cells untouched, the design constants above.
# DIAGNOSTIC: cells not converged are counted per fit and reported; the package's warning
#   for it is muffled, any other warning surfaces.

stopifnot(file.exists("DESCRIPTION"),                                   # BLOCKING
          read.dcf("DESCRIPTION", fields = "Package")[1, 1] == "mvImputer")
source("tests/testthat/helper-sim.R")
source("bench/bench-common.R")

N_ROW <- 200L; N_COL <- 100L; SEEDS <- 1:5; WORKERS <- 1L
MIXED_SEED_OFFSET  <- 2000L           # random half of the "mixed" pattern
GRADED_SEED_OFFSET <- 3000L
K_WEIGHT           <- 2               # "graded": cell weight exp(-K_WEIGHT * z)
PATTERNS   <- c("left", "graded", "mixed")
OUT_DIR    <- "bench/output"
CACHE_DIR  <- file.path(OUT_DIR, "cache-detection-limit")
REPORT_PATH <- file.path(OUT_DIR, "detection-limit.md")
stopifnot(identical(HOLE_FRAC, 0.05), identical(HOLE_SEED_OFFSET, 1000L))   # BLOCKING: the design

# --- holes -----------------------------------------------------------------------------------
# Returns the table with gaps and `low`, the cells that are lowest-in-column by construction
# (all gaps for left and graded; the lowest half for mixed).
make_holes <- function(X, pattern, seed) {
  hole_seed <- HOLE_SEED_OFFSET + seed
  if (pattern == "left") {
    H <- punch_holes(X, HOLE_FRAC, "left", seed = hole_seed)
    return(list(H = H, low = is.na(H)))
  }
  if (pattern == "mixed") {
    H <- punch_holes(X, HOLE_FRAC / 2, "left", seed = hole_seed)
    low <- is.na(H)
    sim_seeded(MIXED_SEED_OFFSET + seed, {          # `code` runs here: H is assigned in this frame
      free <- which(!is.na(H)); n_add <- round(HOLE_FRAC / 2 * length(H))
      H[sample(free, n_add)] <- NA
    })
    return(list(H = H, low = low))
  }
  sim_seeded(GRADED_SEED_OFFSET + seed, {
    H <- X
    for (j in seq_len(ncol(X))) {
      z <- as.numeric(scale(log(X[, j])))
      H[sample(nrow(X), round(HOLE_FRAC * nrow(X)), prob = exp(-K_WEIGHT * z)), j] <- NA
    }
  })
  list(H = H, low = is.na(H))
}

# --- scoring ---------------------------------------------------------------------------------
# A: completed table on the ORIGINAL scale; X: truth; H: table with gaps; m: imputed cells.
# Each error is divided by the SD of that column's observed values on the scale of the error,
# so columns count equally. The log-scale RMSE is NA when an imputed value is not positive.
score <- function(A, X, H, m, low) {
  by_col <- function(v) matrix(v, nrow(X), ncol(X), byrow = TRUE)
  s_orig  <- by_col(apply(H, 2, stats::sd, na.rm = TRUE))
  s_log   <- by_col(apply(log(H), 2, stats::sd, na.rm = TRUE))
  col_min <- by_col(apply(H, 2, min, na.rm = TRUE))
  rms <- function(e, w) if (any(w)) sqrt(mean(e[w]^2)) else NA_real_
  e_orig <- (A - X) / s_orig
  positive <- all(A[m] > 0)
  e_log <- if (positive) (log(A) - log(X)) / s_log else NA_real_
  c(log_all = if (positive) rms(e_log, m) else NA_real_,
    log_low = if (positive) rms(e_log, m & low) else NA_real_,
    log_rand = if (positive) rms(e_log, m & !low) else NA_real_,
    orig_all = rms(e_orig, m), bias_orig = mean(e_orig[m]),
    above = mean(A[m] > X[m]), below_min = mean(A[m] < col_min[m]), nonpos = mean(A[m] <= 0))
}

# --- fits ------------------------------------------------------------------------------------
# One mvImpute() call; on the log scale the table is logged first and the result back-
# transformed. Returns the table on the original scale, seconds, iterations and the number of
# cells not converged in the final stage (the package's warning for it is muffled here).
fit <- function(H, seed, scale, ...) {
  Y <- if (scale == "log") log(H) else H
  t0 <- proc.time()[["elapsed"]]
  out <- withCallingHandlers(
    suppressMessages(mvImputer::mvImpute(Y, seed = seed, workers = WORKERS, ...)),
    warning = function(w) if (grepl("not converged after maxIter", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  sec <- proc.time()[["elapsed"]] - t0
  d <- attr(out, "mvImpute")
  list(A = if (scale == "log") exp(out) else out, sec = sec, iter = sum(d$iterations),
       left = nrow(d$not_converged[[length(d$not_converged)]]))
}

# Rows of the result: one per (pattern, seed, method).
run_seed <- function(seed) {
  X <- sim_numeric(N_ROW, N_COL, seed = seed)
  stopifnot(all(X > 0))                                                 # BLOCKING: log needs positive truth
  rows <- list()
  for (pattern in PATTERNS) {
    hl <- make_holes(X, pattern, seed); H <- hl$H; m <- is.na(H)
    stopifnot(any(m), all(H[!m] == X[!m]))                              # BLOCKING
    fits <- list("min fill" = list(A = H, sec = 0, iter = 0L, left = 0L))
    for (j in seq_len(ncol(H))) fits[["min fill"]]$A[m[, j], j] <- min(H[, j], na.rm = TRUE)
    # Original scale and forceZero: the "left" pattern only (see the header).
    if (pattern == "left") {
      fits[["ranger, original"]] <- fit(H, seed, "original", method = "ranger", forceZero = FALSE)
      fits[["pls 2 comp, original, forceZero = TRUE"]]  <- fit(H, seed, "original", method = "pls", forceZero = TRUE)
      fits[["pls 2 comp, original, forceZero = FALSE"]] <- fit(H, seed, "original", method = "pls", forceZero = FALSE)
    }
    fits[["ranger, log"]]     <- fit(H, seed, "log", method = "ranger")
    fits[["pls 2 comp, log"]] <- fit(H, seed, "log", method = "pls", nComp = 2)
    fits[["pls 3 comp, log"]] <- fit(H, seed, "log", method = "pls", nComp = 3)
    fits[["pls 3 comp, log, guess median"]] <- fit(H, seed, "log", method = "pls", nComp = 3,
                                                   guess = "median")
    for (nm in names(fits)) {
      f <- fits[[nm]]
      rows[[length(rows) + 1L]] <- data.frame(
        pattern = pattern, seed = seed, method = nm, t(score(f$A, X, H, m, hl$low)),
        sec = f$sec, iter = f$iter, not_converged = f$left, check.names = FALSE)
    }
  }
  do.call(rbind, rows)
}

cache_path <- function(seed) file.path(CACHE_DIR, sprintf("seed-%d.rds", seed))

read_runs <- function() {
  paths <- cache_path(SEEDS)
  if (!all(file.exists(paths))) stop("cache incomplete: ", paste(paths[!file.exists(paths)], collapse = ", "))  # BLOCKING
  do.call(rbind, lapply(paths, readRDS))
}

# --- report ----------------------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "n/a", formatC(x, format = "f", digits = d))
pct <- function(x) ifelse(is.na(x), "n/a", sprintf("%.0f %%", 100 * x))
md_table <- function(head, rows)
  c(paste0("| ", paste(head, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(head)), collapse = "|"), "|"),
    vapply(rows, function(r) paste0("| ", paste(r, collapse = " | "), " |"), ""))

# Mean over seeds of one column for (pattern, method); NA if any seed is NA.
avg <- function(runs, pattern, method, col, f = mean) {
  v <- runs[runs$pattern == pattern & runs$method == method, col]
  stopifnot(length(v) == length(SEEDS))                                 # BLOCKING
  f(v)
}

# Seeds (of 5) where `method` has a smaller log-scale RMSE than the minimum fill.
wins <- function(runs, pattern, method) {
  a <- runs[runs$pattern == pattern & runs$method == method, ]; a <- a[order(a$seed), ]
  b <- runs[runs$pattern == pattern & runs$method == "min fill", ]; b <- b[order(b$seed), ]
  sum(a$log_all < b$log_all)
}

table_a <- function(runs) {
  meth <- c("min fill", "ranger, original", "pls 2 comp, original, forceZero = TRUE",
            "pls 2 comp, original, forceZero = FALSE", "ranger, log", "pls 2 comp, log")
  rows <- lapply(meth, function(mt) {
    g <- function(col, f = mean) avg(runs, "left", mt, col, f)
    log_rmse <- if (is.na(g("log_all"))) "n/a" else sprintf("%s (%s)", fmt(g("log_all")), fmt(g("log_all", stats::sd)))
    c(mt, log_rmse, fmt(g("orig_all")), fmt(g("bias_orig")), pct(g("above")), pct(g("below_min")),
      pct(g("nonpos")), fmt(g("not_converged"), 1), if (mt == "min fill") "-" else fmt(g("sec"), 1))
  })
  md_table(c("Method (pattern left)", "Log-scale RMSE (SD over seeds)", "Original-scale RMSE",
             "Mean error, original scale", "Above truth", "Below observed min.",
             "Zero or negative", "Cells not converged", "Seconds"), rows)
}

table_b <- function(runs) {
  meth <- c("min fill", "ranger, log", "pls 2 comp, log", "pls 3 comp, log",
            "pls 3 comp, log, guess median")
  rows <- list()
  for (p in PATTERNS) for (mt in meth) {
    g <- function(col, f = mean) avg(runs, p, mt, col, f)
    rows[[length(rows) + 1L]] <- c(p, mt, sprintf("%s (%s)", fmt(g("log_all")), fmt(g("log_all", stats::sd))),
      if (p == "mixed") fmt(g("log_low")) else "-", if (p == "mixed") fmt(g("log_rand")) else "-",
      fmt(g("orig_all")), pct(g("above")), pct(g("below_min")),
      if (mt == "min fill") "-" else fmt(g("iter"), 1), fmt(g("not_converged"), 1),
      if (mt == "min fill") "-" else fmt(g("sec"), 1))
  }
  md_table(c("Pattern", "Method (log scale)", "Log-scale RMSE (SD over seeds)", "of which lowest cells",
             "of which random cells", "Original-scale RMSE", "Above truth", "Below observed min.",
             "Iterations", "Cells not converged", "Seconds"), rows)
}

table_wins <- function(runs) {
  meth <- c("ranger, log", "pls 2 comp, log", "pls 3 comp, log", "pls 3 comp, log, guess median")
  rows <- lapply(meth, function(mt) c(mt, vapply(PATTERNS, function(p)
    sprintf("%d of %d", wins(runs, p, mt), length(SEEDS)), "")))
  md_table(c("Method (log scale)", PATTERNS), rows)
}

write_report <- function(runs, path) {
  L <- function(p, mt, col = "log_all") avg(runs, p, mt, col)
  secs_ratio <- vapply(PATTERNS, function(p) L(p, "ranger, log", "sec") / L(p, "pls 2 comp, log", "sec"), 0)
  nz <- function(mt) L("left", mt, "nonpos")
  body <- c(
    "# Detection-limit simulation: which settings impute low values best", "",
    "## Summary", "",
    sprintf("- Simulated data (200 rows x 100 columns, log-normal, 5 data sets), 5 %% of cells missing, imputed on the log scale. **Strictly the lowest values missing (`left`):** the column-minimum fill was closest (%s), then PLS (%s), then ranger (%s); PLS (2 components) beat the minimum fill in %d of %d seeds.",
            fmt(L("left", "min fill")), fmt(L("left", "pls 2 comp, log")), fmt(L("left", "ranger, log")),
            wins(runs, "left", "pls 2 comp, log"), length(SEEDS)),
    sprintf("- **Low values more likely, not only the lowest, missing (`graded`, a soft limit):** PLS was closest (%s, better than the minimum fill in %d of %d seeds); minimum fill %s; ranger %s.",
            fmt(L("graded", "pls 2 comp, log")), wins(runs, "graded", "pls 2 comp, log"), length(SEEDS),
            fmt(L("graded", "min fill")), fmt(L("graded", "ranger, log"))),
    sprintf("- **Half lowest-in-column, half at random (`mixed`):** PLS %s, ranger %s, minimum fill %s. The minimum fill is good on the lowest cells (%s) and poor on the random ones (%s); PLS is the reverse (%s and %s).",
            fmt(L("mixed", "pls 2 comp, log")), fmt(L("mixed", "ranger, log")), fmt(L("mixed", "min fill")),
            fmt(L("mixed", "min fill", "log_low")), fmt(L("mixed", "min fill", "log_rand")),
            fmt(L("mixed", "pls 2 comp, log", "log_low")), fmt(L("mixed", "pls 2 comp, log", "log_rand"))),
    sprintf("- **Scale and forceZero (pattern `left`):** on the original scale ranger gave original-scale RMSE %s against %s on the log scale; PLS on original-scale intensities gave %s with `forceZero = TRUE` (%s of imputed cells set to exactly zero) and %s with `forceZero = FALSE` (%s zero or negative). Log-scale PLS was %s. Log-transforming first lowered the original-scale error of both methods, and PLS on untransformed intensities gave many zero or negative values.",
            fmt(L("left", "ranger, original", "orig_all")), fmt(L("left", "ranger, log", "orig_all")),
            fmt(L("left", "pls 2 comp, original, forceZero = TRUE", "orig_all")), pct(nz("pls 2 comp, original, forceZero = TRUE")),
            fmt(L("left", "pls 2 comp, original, forceZero = FALSE", "orig_all")), pct(nz("pls 2 comp, original, forceZero = FALSE")),
            fmt(L("left", "pls 2 comp, log", "orig_all"))),
    sprintf("- **Forests never go below the observed minimum** (share of imputed cells below it in `left`: %s with ranger on the original scale, %s on the log scale): a forest predicts averages of observed values. PLS on log data did go below it (%s of its imputed cells in `left`).",
            pct(L("left", "ranger, original", "below_min")), pct(L("left", "ranger, log", "below_min")),
            pct(L("left", "pls 2 comp, log", "below_min"))),
    sprintf("- **Speed (one worker):** ranger took %s to %s times as long as PLS (%s s against %s s in `left`).",
            fmt(min(secs_ratio), 1), fmt(max(secs_ratio), 1), fmt(L("left", "ranger, log", "sec"), 1),
            fmt(L("left", "pls 2 comp, log", "sec"), 1)),
    sprintf("- **PLS components and `guess`:** 2 and 3 components differ little (log RMSE %s against %s in `graded`); `guess = \"median\"` against `\"minVar\"` (default) with 3 components: %s against %s in `graded`.",
            fmt(L("graded", "pls 2 comp, log")), fmt(L("graded", "pls 3 comp, log")),
            fmt(L("graded", "pls 3 comp, log, guess median")), fmt(L("graded", "pls 3 comp, log"))),
    "- Only simulated, low-rank data, 5 data sets and the settings listed here: see Caveats. These numbers say what happened in this simulation, not what will happen on a given real data set.", "",
    "## Definitions", "",
    "- **Log-scale RMSE**: root mean squared error over the imputed cells of the log of the imputed value minus the log of the true value, each error divided by the SD of that column's observed log values, pooled over columns. Imputations made on the log scale are exponentiated and logged again, so the same cells are compared. Not defined (n/a) when an imputed value is zero or negative.",
    "- **Original-scale RMSE**: the same with the values themselves, each error divided by the SD of that column's observed values. It is dominated by the largest values, which are rarely the missing ones.",
    "- **Above truth**: share of imputed cells whose imputed value is larger than the true value. **Below observed minimum**: share of imputed cells below the smallest observed value of their column. **Zero or negative**: share of imputed cells at or below zero after back-transformation.",
    "- **Mean error, original scale**: mean of the same scaled errors, with sign (positive: imputed too high).",
    "- **Lowest cells / random cells** (`mixed` only): the log-scale RMSE over the cells that are lowest-in-column by construction, and over the cells removed at random.",
    "- **Cells not converged**: imputed cells still changing after the iteration cap in the final stage (mean over seeds). **Seconds**: elapsed time of the `mvImpute()` call, one worker, no worker pool; mean over seeds.",
    "- **Minimum fill**: each gap set to the smallest observed value of its column.", "",
    "## Setup", "",
    sprintf("- Data: `sim_numeric(%d, %d, seed = s)`, s = %s (log-normal, rank-3 signal plus noise). Gaps: 5 %% of cells. `left`: the lowest 5 %% of every column (deterministic). `graded`: in each column 5 %% of rows drawn with weight exp(-%d z), z the standardised log value (hole seed %d + s). `mixed`: the lowest 2.5 %% of each column plus 2.5 %% of the remaining cells at random (hole seeds %d + s and %d + s).",
            N_ROW, N_COL, paste(range(SEEDS), collapse = " to "), K_WEIGHT, GRADED_SEED_OFFSET, HOLE_SEED_OFFSET, MIXED_SEED_OFFSET),
    "- `mvImpute()` with its defaults (100 trees, mtry 5, `tol = c(0.05, 0.025)`, `maxIter = c(15, 60)`, `guess = \"minVar\"`) except what the row says; `seed = s`, `workers = 1`. Log-scale rows impute `log(table)` and exponentiate.",
    sprintf("- %s; mvImputer %s; ranger %s; mixOmics %s; %d logical cores. One machine; other jobs may have run at the same time, so seconds are indicative.",
            R.version.string, utils::packageVersion("mvImputer"), utils::packageVersion("ranger"),
            utils::packageVersion("mixOmics"), parallel::detectCores()),
    "- Script: `bench/detection_limit.R`. It replaces two exploratory scripts (not committed) whose designs it keeps; their log-scale `left` fits are not run twice.", "",
    "## Experiment 1: scale, engine and forceZero (pattern `left`)", "",
    "Log-scale RMSE is n/a where an imputed value was zero or negative (log undefined); the original-scale RMSE then still applies.", "",
    table_a(runs), "",
    "## Experiment 2: three missingness patterns, log scale", "",
    table_b(runs), "",
    "## Seeds (of 5) where the method beat the minimum fill (log-scale RMSE)", "",
    table_wins(runs), "",
    "## Caveats", "",
    "- The data are simulated, low-rank (rank 3), log-normal. Their columns are strongly linearly related on the log scale, which favours PLS, a linear method. Real data with weaker or non-linear structure may favour the forest.",
    "- `left` (exactly the lowest 5 % of every column) is the most favourable case for a minimum fill: the missing values really are the lowest ones. Real detection limits are rarely this clean.",
    "- Only 2 and 3 PLS components were tried, and ranger is the only forest backend (no `randomForest` runs here). Other tuning values (trees, mtry, tol) were not varied.",
    "- Five data sets per pattern, one simulated structure, one machine. The SD over seeds is shown for the headline error; a difference of a few hundredths is within it.",
    "- Seconds are for one serial run (`workers = 1`) and include no worker pool; with several workers the ratio between methods may differ. Other jobs may have run on the machine.",
    "- The simulation knows the truth, so these RMSEs cannot be computed on real data; there, the check is how imputed values behave (for example whether they fall below the observed range), not a distance to the truth.", "")
  writeLines(body, path, useBytes = TRUE)
  invisible(NULL)
}

# --- main --------------------------------------------------------------------------------------
main <- function(args) {
  mode <- if (length(args)) args[1] else "run"
  stopifnot(mode %in% c("run", "report"))                               # BLOCKING
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  if (mode == "run") for (s in SEEDS) {
    if (file.exists(cache_path(s))) next
    t0 <- proc.time()[["elapsed"]]
    saveRDS(run_seed(s), cache_path(s))
    message(format(Sys.time(), "%H:%M:%S"), " seed ", s, " done in ",
            round(proc.time()[["elapsed"]] - t0), " s")
  }
  write_report(read_runs(), REPORT_PATH)
}

main(commandArgs(trailingOnly = TRUE))
