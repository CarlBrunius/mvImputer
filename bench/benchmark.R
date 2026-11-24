# benchmark.R -- speed and accuracy of mvImpute() against the old engine (Task 10)
#
# Purpose: time mvImpute() (final defaults: tol = c(0.05, 0.025), common random numbers on)
#   and StatTools::mvImpWrap() at the same worker count (2 and 4) on three simulated
#   scenarios, score both against the simulated truth, and write bench/output/benchmark.md.
# Estimand: per scenario x worker count, ONE run of each engine on one simulated data set
#   (data seed 1): the elapsed wall seconds of the engine call, the ratio new/old of those
#   seconds, the RMSE of the imputed numeric cells against the truth (raw, and scaled per
#   column by the SD of its observed values) and, for the mixed table, the share of imputed
#   factor-like cells with the true label. A single run has no uncertainty estimate; the
#   replicate-level accuracy comparison is in bench/output/calibration.md.
# Usage (from the package root, after R CMD INSTALL --no-docs .):
#   Rscript bench/benchmark.R smoke                      # small tables, no cache, no report
#   Rscript bench/benchmark.R run <engine> <scenario> <workers>   # one timed run into the cache
#   Rscript bench/benchmark.R all                        # every missing run, one at a time, then the report
#   Rscript bench/benchmark.R report                     # the report from a complete cache
# Dependencies: mvImputer (installed: multisession workers load the installed package);
#   StatTools 0.0.917 (the old engine, reference only); withr (scoped seed, via bench-common.R).
#
# Timing protocol (the machine is shared with other jobs):
# - One R process per timed run, never two at once; the driver ("all") only spawns them.
# - A run starts its clock only when total CPU load is below LOAD_MAX percent; the load just
#   before and just after the run is recorded, so a contaminated run can be seen as such.
# - Each run process lowers its OWN priority to BelowNormal before building anything; the
#   PSOCK / multisession workers it starts inherit it, so both engines run at the same
#   priority and other sessions are not disturbed.
# - Clock: proc.time() elapsed around the engine call alone. Both engines start their own
#   worker pool inside the call, so pool start-up is inside both timings. Loading the
#   packages, building the data, decoding and scoring are outside.
# - The engine that goes first alternates between scenarios (ENGINE_ORDER), so slow drift
#   in the machine's load does not favour one engine throughout.
#
# Choices recorded here:
# - mvImpute() runs with its defaults (nothing but seed and workers is passed); the old
#   engine with mvImpWrap defaults written out (old_args() in bench-common.R).
# - Old engine on the mixed table: one-hot route, as in calibrate_tol.R (bench-common.R).
# - missRanger is not run: it is not installed here and packages are never installed.

stopifnot(file.exists("DESCRIPTION"),                                   # BLOCKING
          read.dcf("DESCRIPTION", fields = "Package")[1, 1] == "mvImputer")
source("tests/testthat/helper-sim.R")
source("bench/bench-common.R")

DATA_SEED   <- 1L                       # data seed; holes use HOLE_SEED_OFFSET + 1; engines seed 1
P_NUMERIC   <- 1000L                    # columns of the numeric scenarios (200 rows)
WORKER_SET  <- c(2L, 4L)
SCENARIOS   <- c(i = "numeric 200 x 1000, 5 % MCAR", ii = "numeric 200 x 1000, 5 % left-censored",
                 iii = "mixed 900 x 20, 5 % MCAR")
ENGINE_ORDER <- list(i = c("new", "old"), ii = c("old", "new"), iii = c("new", "old"))
LOAD_MAX    <- 50                       # percent total CPU load allowed at the start of a run
LOAD_POLL_SECS <- 30
LOAD_WAIT_MAX_SECS <- 7200              # give up (error) rather than time on a busy machine
FLAG_RATIO  <- 1.25
OUT_DIR     <- "bench/output"
CACHE_DIR   <- file.path(OUT_DIR, "cache-benchmark")

# --- machine state (PowerShell; Windows only, as the benchmark's machine is) -----------------
ps <- function(cmd) system2("powershell", c("-NoProfile", "-NonInteractive", "-Command",
                                           shQuote(cmd, type = "cmd")), stdout = TRUE)

# Mean total CPU load over 5 one-second samples. CIM, not performance counters: the counter
# names are translated on non-English Windows.
cpu_load <- function() {
  v <- as.numeric(ps("$s=@(); 1..5 | ForEach-Object { $s += (Get-CimInstance Win32_Processor | Measure-Object LoadPercentage -Average).Average; Start-Sleep 1 }; ($s | Measure-Object -Average).Average"))
  stopifnot(length(v) == 1L, is.finite(v))                        # BLOCKING: a reading was obtained
  v
}

priority_of <- function(pid) ps(sprintf("(Get-Process -Id %d).PriorityClass", pid))

# Lowers this process only (children inherit it); read back so a silent failure is an error.
set_below_normal <- function() {
  ps(sprintf("(Get-Process -Id %d).PriorityClass = 'BelowNormal'", Sys.getpid()))
  stopifnot(identical(priority_of(Sys.getpid()), "BelowNormal"))  # BLOCKING
}

# Failure state: other jobs keep the machine above LOAD_MAX. Waiting is bounded and ends in
# an error, never in a timing taken on a busy machine.
wait_for_quiet <- function() {
  waited <- 0
  repeat {
    load <- cpu_load()
    if (load < LOAD_MAX) return(load)
    if (waited >= LOAD_WAIT_MAX_SECS) stop("CPU load stayed above ", LOAD_MAX, " % for ",
                                           LOAD_WAIT_MAX_SECS, " s; not timing a busy machine")
    message(format(Sys.time(), "%H:%M:%S"), " load ", round(load), " % >= ", LOAD_MAX, ", waiting")
    Sys.sleep(LOAD_POLL_SECS); waited <- waited + LOAD_POLL_SECS
  }
}

# --- one timed run ----------------------------------------------------------------------------
# The final-stage "not converged" warning is a measured outcome (cells_left), so it is
# muffled; any other warning surfaces.
time_new <- function(sc, workers) {
  t0 <- proc.time()[["elapsed"]]
  out <- withCallingHandlers(
    suppressMessages(mvImputer::mvImpute(sc$holes, seed = DATA_SEED, workers = workers)),
    warning = function(w) if (grepl("not converged after maxIter", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  secs <- proc.time()[["elapsed"]] - t0
  d <- attr(out, "mvImpute")
  stopifnot(d$workers == workers)                                 # BLOCKING: ran on the pool asked for
  list(res = out, it1 = d$iterations[1], it2 = d$iterations[2],
       left = nrow(d$not_converged[[2]]), seconds = secs)
}

# Engine call, then scoring. The clock is inside time_new() / old_engine(); the load is read
# immediately before and after those calls.
run_one <- function(engine, scenario, workers, p_numeric = P_NUMERIC) {
  stopifnot(engine %in% c("new", "old"), scenario %in% names(SCENARIOS),
            workers %in% WORKER_SET)                              # BLOCKING
  # BLOCKING: the old engine's ranger threads are its default only if nothing sets them;
  # the thread note in the report depends on it.
  stopifnot(Sys.getenv("R_RANGER_NUM_THREADS") == "", is.null(getOption("ranger.num.threads")),
            is.null(getOption("Ncpus")))
  set_below_normal()
  assert_old_engine()
  sc <- make_scenario(scenario, DATA_SEED, p_numeric)
  load_before <- wait_for_quiet()
  o <- if (engine == "new") time_new(sc, workers) else old_engine(sc$old_input, workers, DATA_SEED)
  load_after <- cpu_load()
  imp <- if (engine == "new") o$res else sc$decode_old(o$res)
  cbind(data.frame(engine = engine, scenario = scenario, workers = workers,
                   seconds = o$seconds, load_before = load_before, load_after = load_after,
                   iter1 = o$it1, iter2 = o$it2, left = o$left,
                   priority = priority_of(Sys.getpid()), finished = format(Sys.time())),
        as.data.frame(as.list(score_run(imp, sc))))
}

cache_path <- function(engine, scenario, workers)
  file.path(CACHE_DIR, sprintf("%s_%s_w%d.rds", engine, scenario, workers))

# The run order: scenarios in order; within a scenario the first engine of ENGINE_ORDER,
# then the other, for each worker count.
build_design <- function() {
  rows <- list()
  for (sc in names(SCENARIOS)) for (w in WORKER_SET) for (e in ENGINE_ORDER[[sc]])
    rows[[length(rows) + 1L]] <- data.frame(engine = e, scenario = sc, workers = w,
                                            stringsAsFactors = FALSE)
  do.call(rbind, rows)
}

# --- smoke run ----------------------------------------------------------------------------------
# Small tables, no cache, no report: checks that every code path runs and that the workers
# inherit the lowered priority (DIAGNOSTIC, printed).
smoke <- function() {
  set_below_normal()
  cl <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  wpid <- parallel::clusterEvalQ(cl, Sys.getpid())[[1]]
  message("worker priority: ", priority_of(wpid))
  rows <- list()
  for (sc in c("i", "ii", "iii")) for (e in c("new", "old")) {
    rows[[length(rows) + 1L]] <- run_one(e, sc, 2L, p_numeric = 30L)
  }
  print(do.call(rbind, rows)[c("engine", "scenario", "workers", "seconds", "load_before",
                               "iter1", "iter2", "left", "rmse_scaled", "factor_acc", "priority")])
}

# --- report ----------------------------------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

read_cache <- function(design) {
  paths <- cache_path(design$engine, design$scenario, design$workers)
  absent <- !file.exists(paths)
  if (any(absent)) stop(sum(absent), " runs missing from the cache, e.g. ", paths[absent][1])  # BLOCKING
  do.call(rbind, lapply(paths, readRDS))
}

# Ratio new / old of wall seconds, per scenario x worker count.
ratios <- function(res) {
  new <- res[res$engine == "new", ]; old <- res[res$engine == "old", ]
  m <- merge(new, old, by = c("scenario", "workers"), suffixes = c("_new", "_old"))
  data.frame(scenario = m$scenario, workers = m$workers, new = m$seconds_new, old = m$seconds_old,
             ratio = m$seconds_new / m$seconds_old)
}

ratio_cell <- function(r) if (r > FLAG_RATIO) sprintf("**FLAG %s**", fmt(r, 2)) else fmt(r, 2)

scenario_table <- function(res, rt, scen) {
  r <- res[res$scenario == scen, ]
  r <- r[order(r$workers, r$engine != "new"), ]
  has_acc <- scen == "iii"
  head <- c("Workers", "Engine", "Goes first", "Load before / after (%)", "Wall seconds",
            "Ratio new/old", "Iter. stage 1 + 2", "Cells left", "Scaled RMSE", "Raw RMSE",
            if (has_acc) "Factor accuracy")
  lines <- c(paste0("| ", paste(head, collapse = " | "), " |"),
             paste0("|", paste(rep("---", length(head)), collapse = "|"), "|"))
  for (k in seq_len(nrow(r))) {
    x <- r[k, ]
    cells <- c(x$workers, x$engine, if (x$engine == ENGINE_ORDER[[scen]][1]) "yes" else "",
               sprintf("%s / %s", fmt(x$load_before, 0), fmt(x$load_after, 0)), fmt(x$seconds, 0),
               if (x$engine == "new") ratio_cell(rt$ratio[rt$scenario == scen & rt$workers == x$workers]) else "",
               sprintf("%d + %d", x$iter1, x$iter2), fmt(x$left, 0), fmt(x$rmse_scaled),
               fmt(x$rmse_raw), if (has_acc) fmt(x$factor_acc))
    lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  lines
}

write_report <- function(res, path) {
  rt <- ratios(res)
  flagged <- rt[rt$ratio > FLAG_RATIO, ]
  desc <- utils::packageDescription("StatTools")
  flag_txt <- if (nrow(flagged)) paste0("**FLAG: ", paste(sprintf("scenario %s at %d workers, ratio %s",
      flagged$scenario, flagged$workers, fmt(flagged$ratio, 2)), collapse = "; "),
      " exceed(s) ", FLAG_RATIO, ".**") else sprintf("No ratio exceeds %s.", FLAG_RATIO)
  sum_rows <- sprintf("(%s, %d workers) %s, new %s s, old %s s", rt$scenario, rt$workers,
                      vapply(rt$ratio, ratio_cell, ""), fmt(rt$new, 0), fmt(rt$old, 0))
  old <- res[res$engine == "old", ]; new <- res[res$engine == "new", ]
  body <- c(
    "# Benchmark: mvImpute() against the old engine", "",
    "## Summary", "",
    sprintf("- Ratio of wall time new/old, one run each: %s.", paste(sum_rows, collapse = "; ")),
    sprintf("- %s A ratio compares whole runs, not cost per iteration: the old engine ran %s iterations (caps 15 + 60 = 75), the new engine %s (it stops a stage when its cells have converged).",
            flag_txt, paste(sprintf("(%s, %dw) %d", old$scenario, old$workers, old$iter1 + old$iter2), collapse = ", "),
            paste(sprintf("(%s, %dw) %d", new$scenario, new$workers, new$iter1 + new$iter2), collapse = ", ")),
    "- This is a single run per cell (data seed 1). The old engine's runs are not reproducible (its parallel workers are not seeded) and timings carry machine noise, so differences of a few percent in time or in RMSE mean nothing. Replicate-level accuracy (5 data sets per scenario, paired) is in `bench/output/calibration.md`.", "",
    "## Definitions", "",
    "- **Wall seconds**: `proc.time()` elapsed around the engine call. Both engines start their own worker pool inside the call, so pool start-up is included for both; data building, decoding and scoring are not.",
    "- **Ratio new/old**: wall seconds of `mvImpute()` divided by wall seconds of `StatTools::mvImpWrap()` at the same scenario and worker count. Above 1.25 it is flagged.",
    "- **Scaled RMSE**: RMSE over all imputed numeric cells against the truth, each error divided by the SD of that column's observed values, pooled over columns (constant columns left out). **Raw RMSE** is the same without scaling. **Factor accuracy** (mixed table) is the share of imputed factor-like cells (factor, ordered, character, logical) with the true label.",
    "- **Load before / after**: total CPU load (%, mean of 5 one-second samples) just before and just after the engine call. A run starts only when the load before is below 50 %.",
    "- **Cells left**: imputed cells not converged at the end of the last stage (old engine: from its log; for the mixed table these are cells of its one-hot indicator matrix, not table cells, so not comparable with the new engine's).", "",
    "## Setup", "",
    sprintf("- Scenarios: (i) `sim_numeric(200, %d)` with 5 %% of cells missing completely at random; (ii) the same data with the lowest 5 %% of each column missing (left-censored); (iii) `sim_mixed(900, extra = 10)` (20 columns of numeric, integer, factor, ordered, character, logical, a constant and 10 noise columns) with 5 %% MCAR over all columns. Data seed %d, hole seed %d, `seed = %d` for `mvImpute()`; the old engine after a scoped `set.seed(%d)` (no effect on its workers).",
            P_NUMERIC, DATA_SEED, HOLE_SEED_OFFSET + DATA_SEED, DATA_SEED, DATA_SEED),
    "- `mvImpute()`: its final defaults (ranger, 100 trees, mtry 5, guess \"minVar\", `tol = c(0.05, 0.025)`, `maxIter = c(15, 60)`, common random numbers on); only `seed` and `workers` are passed.",
    sprintf("- Old engine: StatTools %s, commit `%s`; `mvImpWrap(method = \"RF\", rfMeth = \"ranger\", tol1 = 0.05, n1 = 15, tol2 = 0.025, n2 = 60, nCore = <workers>)`. Mixed table through the one-hot route: character and logical columns made factors, each factor (ordered included) turned into K indicator columns, decoded by the largest indicator (ties to the first level).",
            desc$Version, if (is.null(desc$RemoteSha)) "not recorded" else desc$RemoteSha),
    "- **Thread note**: the old engine calls ranger without `num.threads`, and ranger's default is 2 threads (no `R_RANGER_NUM_THREADS`, `ranger.num.threads` or `Ncpus` was set; checked at each run), so it runs 2 ranger threads in each of its workers: 2 workers use up to 4 threads, 4 workers up to 8. `mvImpute()` fits each forest with 1 thread, so the same worker count gives it half the threads. The comparison is at equal worker count, not equal thread count.",
    "- **Iterations**: the old engine ran to both caps (75 iterations) in the runs in `calibration.md`, as it does here wherever its count is 75; in the mixed table that follows from its stopping rule on indicator columns and is not a meaningful reference cost (see the note in `calibration.md`). The new engine stops when converged.",
    "- **Timing conditions**: one R process per run, runs one at a time; each run waited for total CPU load below 50 % before starting and ran at BelowNormal priority (own process tree only, children inherit), as did the other engine; the engine that goes first alternates between scenarios (the 'Goes first' column). The machine is shared, so the load columns are part of the result.",
    "- Outside reference: `missRanger` is not installed here, so there is no missRanger row (packages are never installed by this benchmark).",
    sprintf("- %s; mvImputer %s; ranger %s; future %s; %d logical cores.", R.version.string,
            utils::packageVersion("mvImputer"), utils::packageVersion("ranger"),
            utils::packageVersion("future"), parallel::detectCores()), "")
  for (sc in names(SCENARIOS))
    body <- c(body, sprintf("## Scenario %s: %s", sc, SCENARIOS[[sc]]), "",
              scenario_table(res, rt, sc), "")
  writeLines(body, path, useBytes = TRUE)
  invisible(rt)
}

# --- main --------------------------------------------------------------------------------------------
# "all" is a driver: it spawns one Rscript per missing run, in design order, waiting for each
# (so never two timed processes at once), and stops at the first failure.
main <- function(args) {
  mode <- if (length(args)) args[1] else "all"
  stopifnot(mode %in% c("smoke", "run", "all", "report"))             # BLOCKING
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  if (mode == "smoke") return(smoke())
  if (mode == "run") {
    stopifnot(length(args) == 4L)                                     # BLOCKING
    r <- run_one(args[2], args[3], as.integer(args[4]))
    saveRDS(r, cache_path(args[2], args[3], as.integer(args[4])))
    message(format(Sys.time(), "%H:%M:%S"), " ", args[2], " ", args[3], " ", args[4], " workers: ",
            round(r$seconds, 1), " s, load before ", round(r$load_before), " %")
    return(invisible(r))
  }
  design <- build_design()
  if (mode == "all") for (k in seq_len(nrow(design))) {
    d <- design[k, ]
    if (file.exists(cache_path(d$engine, d$scenario, d$workers))) next
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("bench/benchmark.R", "run", d$engine, d$scenario, d$workers))
    if (status != 0L) stop("run failed: ", d$engine, " ", d$scenario, " ", d$workers)  # BLOCKING
  }
  res <- read_cache(design)
  utils::write.csv(res, file.path(OUT_DIR, "benchmark-runs.csv"), row.names = FALSE)
  write_report(res, file.path(OUT_DIR, "benchmark.md"))
}

main(commandArgs(trailingOnly = TRUE))
