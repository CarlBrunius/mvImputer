# calibrate_tol.R -- choose the default convergence tolerances of mvImpute() (spec 5.4)
#
# Purpose: run mvImpute() over a grid of (tol1, tol2) on three simulated scenarios and
#   compare each setting with the old engine (StatTools::mvImpWrap with its defaults),
#   then apply the selection rule and write bench/output/calibration.md.
# Estimand: per scenario and setting, the mean over replicate data sets of the RMSE of the
#   imputed numeric cells against the simulated truth, each column scaled by the SD of its
#   observed values (the stopping rule's scale), and the mean total number of iterations
#   (both stages). Selection rule (spec 5.4): per scenario, the setting with the smallest
#   mean scaled RMSE among those whose mean total iterations do not exceed the old engine's.
# Usage (from the package root, after R CMD INSTALL --no-docs .):
#   Rscript bench/calibrate_tol.R probe   # one replicate, default tol + old engine; timings
#   Rscript bench/calibrate_tol.R         # full grid, 5 replicates; writes the report
#   Rscript bench/calibrate_tol.R run i 3 # one slice (scenario, replicate) into the cache
#   Rscript bench/calibrate_tol.R report  # the report from a complete cache
# Dependencies: mvImputer (installed: multisession workers load the installed package);
#   StatTools 0.0.917 (the old engine, the reference); future (one 2-worker pool for all
#   mvImpute() calls); withr (scoped seeds). Test helpers and bench/bench-common.R are source()d by path.
#
# Choices recorded here (Carl to confirm):
# - Selection metric is the SCALED RMSE, pooled over all imputed numeric cells of a data
#   set: raw RMSE on log-normal intensities is dominated by the few high-abundance columns.
#   Raw RMSE and (scenario iii) factor accuracy are reported beside it, not used to select.
# - Both engines start from their default first guess ("minVar", column minimum), also in
#   the mixed scenario; mvImpute() otherwise runs with all its defaults except tol.
# - Old engine, mixed scenario: character and logical columns become factors, each factor
#   becomes K indicator columns (model.matrix, NA rows kept), decoded by the largest
#   indicator, ties to the first level (spec section 9, the route earlier analyses used).
# - The old engine's workers are not seeded, so its runs do not reproduce: only its means
#   over replicates are meaningful. mvImpute() runs are reproducible from their seed.
# - Finished runs are cached in bench/output/cache-calibration/ so an interrupted run
#   resumes; delete that folder after any change to this script or the package.

stopifnot(file.exists("DESCRIPTION"),                                   # BLOCKING
          read.dcf("DESCRIPTION", fields = "Package")[1, 1] == "mvImputer")
source("tests/testthat/helper-sim.R")
source("bench/bench-common.R")     # scenarios, scoring and the old-engine runner, shared with benchmark.R

TOL1_GRID  <- c(0.1, 0.05, 0.025)
TOL2_GRID  <- c(0.05, 0.025, 0.01, 0.005)
MAX_ITER   <- c(15, 60)                       # mvImpute() default, equal to the old n1, n2
REP_SEEDS  <- 1:5                             # data, holes and engine seeds derive from these
WORKERS    <- 2L
P_NUMERIC  <- 500L                            # columns of the numeric scenarios
OUT_DIR    <- "bench/output"
CACHE_DIR  <- file.path(OUT_DIR, "cache-calibration")
# Carl's decision of 2026-11-24, recorded in the report (the grid does not choose it).
DECISION   <- "tol = c(0.05, 0.025)"
# Why the old engine's numbers in the mixed table are not a reference cost (from its code:
# it stops a cell when abs(new - old) / old < tol, old being the current value).
OLD_ENGINE_III_NOTE <- paste(
  "Note on the old engine in this scenario: its \"Cells left\" counts cells of the one-hot",
  "indicator matrix, not table cells, so it is not comparable with the new engine's column.",
  "Its stopping test divides the change by the current value; indicator cells start at 0",
  "(the column minimum) and stay near 0, so the test is rarely met for them, while any",
  "negative current value counts as converged. Its 75 iterations here follow from that rule",
  "and are not a meaningful reference cost. In scenarios i and ii (positive data) its caps",
  "are genuine.")
SCENARIOS  <- c(i = "numeric 200 x 500, 5 % MCAR", ii = "numeric 200 x 500, 5 % left-censored",
                iii = "mixed 900 x 20, 5 % MCAR")

assert_old_engine()                              # BLOCKING (bench-common.R)

# --- engines ------------------------------------------------------------------------------
# The final-stage "not converged" warning is a measured outcome here (stage2_at_cap), so it
# is muffled; any other warning surfaces.
run_new <- function(sc, tol1, tol2, seed) {
  t0 <- proc.time()[["elapsed"]]
  out <- withCallingHandlers(
    suppressMessages(mvImputer::mvImpute(sc$holes, seed = seed, workers = NULL,
                                         tol = c(tol1, tol2), maxIter = MAX_ITER)),
    warning = function(w) if (grepl("not converged after maxIter", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  secs <- proc.time()[["elapsed"]] - t0
  it <- attr(out, "mvImpute")$iterations
  stopifnot(attr(out, "mvImpute")$workers == WORKERS)             # BLOCKING: ran on the pool
  left <- nrow(attr(out, "mvImpute")$not_converged[[2]])
  # At cap = cap reached WITH cells left: a stage that converges on its last iteration
  # has also run maxIter iterations, but it did converge.
  c(score_run(out, sc), iter1 = it[1], iter2 = it[2], iter_total = sum(it),
    stage2_at_cap = it[2] == MAX_ITER[2] && left > 0, left = left, seconds = secs)
}

# The old engine is run and its log read by old_engine() (bench-common.R): set.seed(rep)
# before the call as the brief asks, scoped, and iterations / cells left from its printed log.
run_old <- function(sc, seed) {
  o <- old_engine(sc$old_input, WORKERS, seed)
  c(score_run(sc$decode_old(o$res), sc), iter1 = o$it1, iter2 = o$it2,
    iter_total = o$it1 + o$it2, stage2_at_cap = o$it2 == o$n2 && o$left > 0, left = o$left,
    seconds = o$seconds)
}

# One cached run: the key names engine, scenario, setting and replicate.
cached <- function(key, expr) {
  path <- file.path(CACHE_DIR, paste0(key, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  val <- expr
  saveRDS(val, path)
  val
}

run_design <- function(design) {
  rows <- vector("list", nrow(design))
  for (k in seq_len(nrow(design))) {
    d <- design[k, ]
    sc <- make_scenario(d$scenario, d$rep, P_NUMERIC)
    key <- sprintf("%s_%s_%s_%s_r%d", d$engine, d$scenario, d$tol1, d$tol2, d$rep)
    res <- cached(key, if (d$engine == "old") run_old(sc, d$rep) else
      run_new(sc, d$tol1, d$tol2, d$rep))
    message(format(Sys.time(), "%H:%M:%S"), " ", key, ": ", round(res[["seconds"]], 1), " s")
    rows[[k]] <- cbind(d, as.data.frame(as.list(res)))
  }
  do.call(rbind, rows)
}

# --- design --------------------------------------------------------------------------------
build_design <- function(probe) {
  reps <- if (probe) REP_SEEDS[1] else REP_SEEDS
  grid <- if (probe) data.frame(tol1 = 0.05, tol2 = 0.025) else
    expand.grid(tol1 = TOL1_GRID, tol2 = TOL2_GRID)
  new <- merge(expand.grid(scenario = names(SCENARIOS), rep = reps, stringsAsFactors = FALSE), grid)
  old <- expand.grid(scenario = names(SCENARIOS), rep = reps, tol1 = NA, tol2 = NA,
                     stringsAsFactors = FALSE)
  rbind(cbind(engine = "old", old), cbind(engine = "new", new))
}

# --- summary and selection ---------------------------------------------------------------
summarise_runs <- function(res) {
  keys <- unique(res[c("engine", "scenario", "tol1", "tol2")])
  rows <- lapply(seq_len(nrow(keys)), function(k) {
    r <- merge(keys[k, ], res)
    stopifnot(nrow(r) == length(REP_SEEDS))                       # BLOCKING: every replicate present
    data.frame(keys[k, ], rmse = mean(r$rmse_scaled), rmse_sd = stats::sd(r$rmse_scaled),
               rmse_raw = mean(r$rmse_raw), acc = mean(r$factor_acc),
               iter1 = mean(r$iter1), iter2 = mean(r$iter2), iter = mean(r$iter_total),
               iter_min = min(r$iter_total), iter_max = max(r$iter_total),
               # at cap = cap reached with cells left, derived here so cached runs follow it
               at_cap = sum(r$iter2 == MAX_ITER[2] & r$left > 0), left = mean(r$left),
               seconds = mean(r$seconds))
  })
  do.call(rbind, rows)
}

# Spec 5.4: per scenario, the smallest mean scaled RMSE among settings whose mean total
# iterations do not exceed the old engine's mean. NA when no setting qualifies.
select_setting <- function(s) {
  old_iter <- s$iter[s$engine == "old"]
  new <- s[s$engine == "new", ]
  ok <- new[new$iter <= old_iter, ]
  if (!nrow(ok)) return(NULL)
  ok[which.min(ok$rmse), ]
}

# Paired over replicates (same data in every setting): mean and SD of the per-replicate
# difference in scaled RMSE from the selected setting. Shows whether a gap between
# settings is larger than replicate-to-replicate noise.
paired_diff <- function(res, scen, pick) {
  new <- res[res$engine == "new" & res$scenario == scen, ]
  ref <- new[new$tol1 == pick$tol1 & new$tol2 == pick$tol2, c("rep", "rmse_scaled")]
  m <- merge(new, ref, by = "rep", suffixes = c("", "_ref"))
  d <- m$rmse_scaled - m$rmse_scaled_ref
  agg <- stats::aggregate(d, by = list(tol1 = m$tol1, tol2 = m$tol2),
                          FUN = function(x) c(mean = mean(x), sd = stats::sd(x)))
  data.frame(tol1 = agg$tol1, tol2 = agg$tol2, d_mean = agg$x[, "mean"], d_sd = agg$x[, "sd"])
}

# Builder's aggregations (not the spec rule; reported as such). The rule can only trade
# accuracy against iterations through the old engine's count. "Efficient" = among settings
# whose mean paired scaled-RMSE difference to the rule's pick is within 2 standard errors
# (SD / sqrt(replicates)), the one with the fewest mean total iterations (ties: larger tol).
efficient_setting <- function(s, res, scen, pick) {
  pd <- paired_diff(res, scen, pick)
  ok <- pd[pd$d_mean <= 2 * pd$d_sd / sqrt(length(REP_SEEDS)) | pd$d_mean <= 0, ]
  new <- merge(s[s$engine == "new" & s$scenario == scen, ], ok[c("tol1", "tol2")])
  new[order(new$iter, -new$tol1, -new$tol2), ][1, ]
}

# The setting with fewer mean iterations than `chosen` that came closest to passing the
# 2-SE screen, and by how much it missed (NULL if none): shows how fragile the pick is.
efficient_near_miss <- function(s, res, scen, pick, chosen) {
  pd <- paired_diff(res, scen, pick)
  pd$miss_by <- pd$d_mean - 2 * pd$d_sd / sqrt(length(REP_SEEDS))
  new <- merge(s[s$engine == "new" & s$scenario == scen, c("tol1", "tol2", "iter")], pd)
  out <- new[new$miss_by > 0 & new$d_mean > 0 & new$iter < chosen$iter, ]
  if (!nrow(out)) return(NULL)
  out[which.min(out$miss_by), ]
}

# How unstable the per-scenario ranking is: the number of distinct settings that are best
# in some replicate.
n_rep_best <- function(res, scen) {
  new <- res[res$engine == "new" & res$scenario == scen, ]
  best <- vapply(split(new, new$rep), function(x)
    paste(x$tol1, x$tol2)[which.min(x$rmse_scaled)], "")
  length(unique(best))
}

# One default for all scenarios: the setting whose worst (over scenarios) mean paired
# scaled-RMSE excess over that scenario's pick is smallest (minimax regret).
common_setting <- function(res, picks) {
  pds <- lapply(names(picks), function(sc) {
    d <- paired_diff(res, sc, picks[[sc]]); d$scenario <- sc; d })
  pd <- do.call(rbind, pds)
  worst <- stats::aggregate(d_mean ~ tol1 + tol2, pd, max)
  worst[order(worst$d_mean), ][1:2, ]          # best and runner-up, for the margin
}

# --- report ------------------------------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

scenario_table <- function(s, res, scen, pick) {
  s <- s[s$scenario == scen, ]
  s <- s[order(s$engine != "old", -s$tol1, -s$tol2), ]
  pd <- if (!is.null(pick)) paired_diff(res, scen, pick) else NULL
  has_acc <- scen == "iii"
  head <- c("Engine", "tol1", "tol2", "Scaled RMSE (SD)", "Diff. to selected (SD)", "Raw RMSE",
            if (has_acc) "Factor accuracy", "Iter. stage 1 + 2", "Total iter. (range)",
            "Stage 2 at cap", "Cells left (final stage)", "Seconds")
  lines <- c(paste0("| ", paste(head, collapse = " | "), " |"),
             paste0("|", paste(rep("---", length(head)), collapse = "|"), "|"))
  for (k in seq_len(nrow(s))) {
    r <- s[k, ]; is_old <- r$engine == "old"
    sel <- !is.null(pick) && !is_old && r$tol1 == pick$tol1 && r$tol2 == pick$tol2
    d <- if (is_old || is.null(pd)) "" else {
      p <- pd[pd$tol1 == r$tol1 & pd$tol2 == r$tol2, ]
      sprintf("%s (%s)", fmt(p$d_mean, 4), fmt(p$d_sd, 4))
    }
    cells <- c(if (is_old) "old" else if (sel) "**new (selected)**" else "new",
               if (is_old) "0.05*" else r$tol1, if (is_old) "0.025*" else r$tol2,
               sprintf("%s (%s)", fmt(r$rmse), fmt(r$rmse_sd)), d, fmt(r$rmse_raw),
               if (has_acc) fmt(r$acc), sprintf("%s + %s", fmt(r$iter1, 1), fmt(r$iter2, 1)),
               sprintf("%s (%d-%d)", fmt(r$iter, 1), r$iter_min, r$iter_max),
               sprintf("%d / %d", r$at_cap, length(REP_SEEDS)), fmt(r$left, 1), fmt(r$seconds, 0))
    lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  lines
}

write_report <- function(res, path) {
  s <- summarise_runs(res)
  picks <- lapply(names(SCENARIOS), function(sc) select_setting(s[s$scenario == sc, ]))
  names(picks) <- names(SCENARIOS)
  stopifnot(!vapply(picks, is.null, TRUE))   # BLOCKING: the proposal below needs a pick per scenario
  eff <- lapply(names(SCENARIOS), function(sc) efficient_setting(s, res, sc, picks[[sc]]))
  names(eff) <- names(SCENARIOS)
  common2 <- common_setting(res, picks)
  common <- common2[1, ]; runner <- common2[2, ]
  near <- lapply(names(SCENARIOS), function(sc)
    efficient_near_miss(s, res, sc, picks[[sc]], eff[[sc]]))
  names(near) <- names(SCENARIOS)
  near <- near[!vapply(near, is.null, TRUE)]
  near_txt <- if (!length(near)) "" else paste0(" Sensitivity: ", paste(vapply(names(near), function(sc)
    sprintf("in (%s) tol = c(%s, %s) has fewer iterations and misses the 2-SE screen by only %s",
            sc, near[[sc]]$tol1, near[[sc]]$tol2, formatC(near[[sc]]$miss_by, format = "e", digits = 1)), ""),
    collapse = "; "), ".")
  rep_best <- vapply(names(SCENARIOS), function(sc) n_rep_best(res, sc), 1L)
  setting <- function(r) sprintf("tol = c(%s, %s)", r$tol1, r$tol2)
  per_scen <- function(lst) paste(sprintf("(%s) %s", names(lst), vapply(lst, setting, "")),
                                  collapse = "; ")
  old_iter <- stats::setNames(s$iter[s$engine == "old"], s$scenario[s$engine == "old"])
  n_excluded <- vapply(names(SCENARIOS), function(sc)
    sum(s$engine == "new" & s$scenario == sc & s$iter > old_iter[[sc]]), 1L)
  spread <- vapply(names(SCENARIOS), function(sc) max(paired_diff(res, sc, picks[[sc]])$d_mean), 0)
  rmse_level <- vapply(names(SCENARIOS), function(sc) mean(s$rmse[s$scenario == sc & s$engine == "new"]), 0)
  old_rmse <- stats::setNames(s$rmse[s$engine == "old"], s$scenario[s$engine == "old"])
  cap_new <- vapply(names(SCENARIOS), function(sc) sum(s$at_cap[s$engine == "new" & s$scenario == sc]), 0)
  n_new <- length(TOL1_GRID) * length(TOL2_GRID) * length(REP_SEEDS)
  desc <- utils::packageDescription("StatTools")
  body <- c(
    "# Calibration of the default tolerances (`tol`)", "",
    "## Summary", "",
    sprintf("- Compared `mvImpute()` over tol1 in {%s} x tol2 in {%s} with the old engine (`StatTools::mvImpWrap()`, its defaults) on three simulated scenarios, %d replicate data sets each. The old engine ran both stages to their caps (mean total iterations %s of a possible %d), so the rule's iteration condition excluded %s of %d settings and the rule reduced to the smallest mean scaled RMSE.",
            paste(TOL1_GRID, collapse = ", "), paste(TOL2_GRID, collapse = ", "), length(REP_SEEDS),
            paste(sprintf("(%s) %s", names(old_iter), fmt(old_iter, 1)), collapse = ", "), sum(MAX_ITER),
            paste(sprintf("(%s) %d", names(n_excluded), n_excluded), collapse = ", "),
            length(TOL1_GRID) * length(TOL2_GRID)),
    sprintf("- Accuracy is flat over the grid: the largest mean paired difference between any setting and the scenario's pick is %s scaled-RMSE units, against scaled RMSEs of about %s (old engine %s). Stage 2 reached its cap in %s of %d new-engine runs; in the mixed table a smaller tol does not remove the cells left.",
            paste(sprintf("(%s) %s", names(spread), fmt(spread, 4)), collapse = ", "),
            paste(sprintf("(%s) %s", names(rmse_level), fmt(rmse_level, 2)), collapse = ", "),
            paste(sprintf("(%s) %s", names(old_rmse), fmt(old_rmse, 3)), collapse = ", "),
            paste(sprintf("(%s) %d", names(cap_new), cap_new), collapse = ", "), n_new),
    sprintf("- **Decision (Carl, 2026-11-24): keep %s, the current defaults.** Accuracy is flat over the grid (largest mean paired gap %s scaled-RMSE units) and the iteration condition never binds, so no setting in the grid is supported over the current one.",
            DECISION, fmt(max(spread), 3)),
    sprintf("- Spec rule picks per scenario: %s. Builder's aggregations (not the spec rule; see Definitions): common setting %s, ahead of the runner-up %s by %s in worst-case excess, while the best setting differs between replicates (%s distinct best settings in %d replicates); efficient setting %s.%s",
            per_scen(picks), setting(common), setting(runner), fmt(runner$d_mean - common$d_mean, 4),
            paste(sprintf("(%s) %d", names(rep_best), rep_best), collapse = ", "), length(REP_SEEDS),
            per_scen(eff), near_txt),
    sprintf("- Uncertainty: %d replicates per setting, and the gaps between settings are of the size of their paired SD, so the ranking within the grid is not stable; the choice is mainly one of iteration cost. The old engine's runs are not reproducible (its parallel workers are not seeded), so its numbers are means over replicates and would move somewhat on a rerun.", length(REP_SEEDS)), "",
    "## Definitions (choices for Carl to confirm)", "",
    "- **Scaled RMSE** (selection metric): per replicate, the root mean squared error over all imputed cells of numeric columns, each error divided by the SD of that column's observed values (the scale the stopping rule uses), pooled over columns; then the mean over replicates. Constant columns are left out (filled exactly). **Raw RMSE** is the same without scaling. **Factor accuracy** (scenario iii) is the share of imputed factor-like cells (factor, ordered, character, logical) with the true label; it is reported, not used to select.",
    "- **Total iterations**: stage 1 + stage 2 iterations. For the old engine, counted from the lines starting `Iteration` in its printed log.",
    "- **Selection rule** (spec section 5.4): per scenario, the setting with the smallest mean scaled RMSE whose mean total iterations do not exceed the old engine's mean.",
    "- **Diff. to selected**: per replicate, the setting's scaled RMSE minus the selected setting's on the same data set; mean (SD) over replicates. Positive means worse than the selected setting.",
    "- **Builder's aggregations** (added for this report; not part of the spec rule, and with margins far below the replicate noise): the *common setting* is the setting whose worst-case (over scenarios) mean paired scaled-RMSE excess over that scenario's pick is smallest; the *efficient setting* of a scenario is, among settings whose mean paired difference to the pick is at most 2 standard errors (SD / sqrt(5)), the one with the fewest mean total iterations (ties: the larger tol).",
    "- **Stage 2 at cap**: replicates whose stage 2 ran to its iteration cap (60) with cells still not converged. **Cells left**: mean number of imputed cells still not converged at the end of stage 2 (old engine: from its log; in scenario iii these are cells of its indicator matrix, see the note under that table).", "",
    "## Setup", "",
    "- Scenarios: (i) `sim_numeric(200, 500)` (log-normal, rank-3 signal plus noise) with 5 % of cells missing completely at random (MCAR); (ii) the same data with the lowest 5 % of each column missing (left-censored, a detection-limit pattern); (iii) `sim_mixed(900, extra = 10)` (20 columns: numeric, integer, factor, ordered factor, character, logical, a constant and 10 noise columns) with 5 % MCAR over all columns.",
    sprintf("- Replicate r (r = %s) uses data seed r, hole seed %d + r, and `seed = r` for `mvImpute()`. The old engine is called after a scoped `set.seed(r)`, as the brief asks; that has no effect on its parallel workers, which seed themselves.",
            paste(REP_SEEDS, collapse = ", "), HOLE_SEED_OFFSET),
    "- `mvImpute()`: defaults (ranger, 100 trees, mtry 5, guess \"minVar\", maxIter 15 and 60) except `tol`; one 2-worker multisession pool per process.",
    sprintf("- Old engine: StatTools %s, commit `%s`; `mvImpWrap(method = \"RF\", rfMeth = \"ranger\", tol1 = 0.05, n1 = 15, tol2 = 0.025, n2 = 60, nCore = 2)`; first guess column minimum (its default). Its tolerance is relative to the current value, not to the column SD, so its tol values (marked *) are not on the same scale as the new ones.",
            desc$Version, if (is.null(desc$RemoteSha)) "not recorded" else desc$RemoteSha),
    "- Old engine, scenario iii: character and logical columns made factors; each factor (including the ordered one) became K indicator columns (`model.matrix(~ f - 1)`, missing rows kept); imputed indicators decoded to the level with the largest value (ties to the first level).",
    sprintf("- %s; mvImputer %s; ranger %s; future %s. Sum of per-run times %s minutes. Runs were spread over several concurrent processes, so the seconds columns include contention and are indicative only; the timing comparison belongs to the benchmark.",
            R.version.string, utils::packageVersion("mvImputer"), utils::packageVersion("ranger"),
            utils::packageVersion("future"), fmt(sum(res$seconds) / 60, 0)), "")
  for (sc in names(SCENARIOS))
    body <- c(body, sprintf("## Scenario %s: %s", sc, SCENARIOS[[sc]]), "",
              scenario_table(s, res, sc, picks[[sc]]), "",
              if (sc == "iii") c(OLD_ENGINE_III_NOTE, ""))
  writeLines(body, path, useBytes = TRUE)
  invisible(list(summary = s, picks = picks, efficient = eff, common = common))
}

# --- main ------------------------------------------------------------------------------------
# Modes: "probe" (timing); "run <scenario> <rep>" (one slice into the cache, no report, so
# slices can run in parallel processes); "report" (reads the cache; every run must be
# there); no argument: every run in this process, then the report.
main <- function(args) {
  mode <- if (length(args)) args[1] else "all"
  stopifnot(mode %in% c("probe", "run", "report", "all"))           # BLOCKING
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  design <- build_design(mode == "probe")
  if (mode == "run")
    design <- design[design$scenario == args[2] & design$rep == as.integer(args[3]), ]
  if (mode == "report") {
    keys <- sprintf("%s_%s_%s_%s_r%d", design$engine, design$scenario, design$tol1,
                    design$tol2, design$rep)
    absent <- !file.exists(file.path(CACHE_DIR, paste0(keys, ".rds")))
    if (any(absent)) stop(sum(absent), " runs missing from the cache, e.g. ", keys[absent][1])  # BLOCKING
  }
  stopifnot(nrow(design) > 0)                                       # BLOCKING: a valid slice
  # One pool for all mvImpute() calls (workers = NULL uses it); results do not depend on
  # the pool, only start-up time does.
  old_plan <- future::plan(future::multisession, workers = WORKERS)
  on.exit(future::plan(old_plan), add = TRUE)
  res <- run_design(design)
  if (mode == "probe") print(res[c("engine", "scenario", "iter_total", "rmse_scaled", "seconds")])
  if (mode %in% c("report", "all")) {
    utils::write.csv(res, file.path(OUT_DIR, "calibration-runs.csv"), row.names = FALSE)
    write_report(res, file.path(OUT_DIR, "calibration.md"))
  }
  invisible(res)
}

main(commandArgs(trailingOnly = TRUE))
