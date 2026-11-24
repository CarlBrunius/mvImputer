# experiments.R -- common random numbers and the co-missing factor pair (spec section 5.5)
#
# Purpose: (a) compare the stage engine with and without common random numbers (crn: every
#   iteration of a stage reuses a column's random stream) on omics-like numeric data;
#   (b) check that two strongly associated factors missing in the same rows do not fall
#   into a 2-cycle (each flipping the other's hard label). Writes bench/output/experiments.md.
# Estimand: (a) per crn setting, iterations per stage and the RMSE of the imputed cells
#   against the simulated truth (scaled: each error divided by the SD of the column's
#   observed values; raw beside it), per seed and as a mean over 5 seeds. (b) per seed and
#   crn setting, whether the final stage reaches maxIter, and per iteration the number of
#   imputed factor labels that changed and that flipped back to their value two iterations
#   earlier (the 2-cycle signature).
# Usage (from the package root, after R CMD INSTALL --no-docs .):
#   Rscript bench/experiments.R          # all runs, then the report
#   Rscript bench/experiments.R report   # the report again from the saved runs
# Dependencies: mvImputer internals via ::: (run_stage() takes crn; mvImpute() does not
#   expose it); future (one 2-worker pool); withr (scoped seeds). Test helpers source()d.
#
# run_engine() below repeats mvImpute()'s wiring (encode, first guess, scales, streams, one
# run_stage() per stage, decode) so crn can be switched. A BLOCKING check shows it gives
# mvImpute()'s exact output with crn = TRUE, so (a) measures the package, not a copy.
#
# Per-iteration labels in (b): run_stage() returns only its final table. Its first k
# iterations do not depend on maxIter (the loop only stops earlier), so the table after
# iteration k is run_stage(..., maxIter = k)$W. That is exact and costs k reruns per stage;
# the data are small (300 rows), so it is affordable.

stopifnot(file.exists("DESCRIPTION"),                                   # BLOCKING
          read.dcf("DESCRIPTION", fields = "Package")[1, 1] == "mvImputer")
source("tests/testthat/helper-sim.R")

SEEDS      <- 1:5
# (b) uses more seeds: a scan of 20 seeds before this script was fixed showed the final stage
# at its cap with crn in about 1 data set in 20, which 5 seeds would likely miss.
SEEDS_PAIR <- 1:20
HOLE_SEED_OFFSET <- 1000L          # holes and the factor pair use seeds apart from the data's
TOL        <- c(0.05, 0.025)       # mvImpute() defaults at the time of the experiment
MAX_ITER   <- c(15, 60)
WORKERS    <- 2L
N_COMISS   <- 40                   # (b) rows where both factors are missing
PERMUTE_FRAC <- 0.10               # (b) share of fac3 labels permuted to make fac3b
OUT_PATH   <- "bench/output/experiments.md"

# --- engine ----------------------------------------------------------------------------
# Returns the decoded table, iterations per stage, final-stage unconverged cells and, when
# `trace_cols` names factor columns, per stage the list of tables after each iteration.
run_engine <- function(X, seed, crn, trace_cols = NULL) {
  # run_stage() advances the caller's .Random.seed and, in a session without one, future
  # can leave RNGkind() at L'Ecuyer-CMRG; the package's guard restores both.
  mvImputer:::local_rng_guard()
  enc     <- mvImputer:::encode_input(X)
  W       <- mvImputer:::first_guess(enc$W, enc, "minVar")
  scales  <- mvImputer:::col_scales(enc$W, enc$miss, enc$type)
  streams <- mvImputer:::make_streams(seed, length(TOL), ncol(W))
  stage <- function(W, s, max_iter) mvImputer:::run_stage(W, enc, scales, streams[[s]], TOL[s],
    max_iter, "ranger", ntree = 100, mtry = 5, nComp = 2, forceZero = FALSE, crn = crn)
  iters <- integer(length(TOL)); traces <- list(); left <- NA_integer_
  for (s in seq_along(TOL)) {
    res <- stage(W, s, MAX_ITER[s])
    iters[s] <- res$iterations
    if (length(trace_cols)) {
      steps <- lapply(seq_len(res$iterations - 1L), function(k) stage(W, s, k)$W)
      traces[[s]] <- c(list(W), steps, list(res$W))   # tables after iteration 0..K
    }
    left <- nrow(res$not_converged)
    W <- res$W
  }
  list(out = mvImputer:::decode_output(W, enc), iterations = iters, left = left,
       traces = traces)
}

rmse_vs_truth <- function(imp, truth, holes) {
  miss <- is.na(holes)
  s <- apply(holes, 2, stats::sd, na.rm = TRUE)
  err <- (imp - truth)[miss]
  c(scaled = sqrt(mean((err / s[col(holes)[miss]])^2)), raw = sqrt(mean(err^2)))
}

# --- (a) common random numbers ----------------------------------------------------------------
experiment_a <- function() {
  rows <- list()
  for (seed in SEEDS) {
    truth <- sim_numeric(200, 300, seed = seed)
    holes <- punch_holes(truth, 0.05, "mcar", seed = HOLE_SEED_OFFSET + seed)
    for (crn in c(TRUE, FALSE)) {
      t0 <- proc.time()[["elapsed"]]
      r <- run_engine(holes, seed, crn)
      secs <- proc.time()[["elapsed"]] - t0
      rows[[length(rows) + 1L]] <- data.frame(seed, crn, iter1 = r$iterations[1],
        iter2 = r$iterations[2], left = r$left, t(rmse_vs_truth(r$out, truth, holes)),
        seconds = secs)
      message(format(Sys.time(), "%H:%M:%S"), " (a) seed ", seed, " crn ", crn, ": ",
              round(secs), " s")
    }
  }
  do.call(rbind, rows)
}

# BLOCKING: run_engine(crn = TRUE) must reproduce mvImpute() exactly, or (a) and (b)
# describe a copy of the engine rather than the package.
check_engine_matches_package <- function() {
  holes <- punch_holes(sim_numeric(60, 15, seed = 99), 0.05, "mcar", seed = 99)
  pkg <- mvImputer::mvImpute(holes, seed = 7, workers = NULL, tol = TOL, maxIter = MAX_ITER)
  attr(pkg, "mvImpute") <- NULL
  stopifnot(identical(pkg, run_engine(holes, 7, crn = TRUE)$out))
}

# --- (b) two associated factors, co-missing ------------------------------------------------
# fac3b copies fac3, then the labels of a random 10 % of rows are permuted among those rows
# (so fac3b is strongly but not perfectly tied to fac3); both go missing in the same 40 rows.
make_factor_pair <- function(seed) {
  d <- sim_mixed(300, seed = seed)
  # sim_seeded() (helper-sim.R) pins the RNG kinds: with_seed() alone would draw under the
  # session's current kind, so the pair would change with whatever ran before.
  sim_seeded(HOLE_SEED_OFFSET + seed, {
    idx <- sample.int(nrow(d), round(PERMUTE_FRAC * nrow(d)))
    d$fac3b <- d$fac3
    d$fac3b[idx] <- d$fac3[idx][sample.int(length(idx))]
    gap <- sample.int(nrow(d), N_COMISS)
  })
  truth <- d
  d$fac3[gap] <- NA; d$fac3b[gap] <- NA
  list(truth = truth, holes = d, agree = mean(truth$fac3 == truth$fac3b))
}

# Per iteration k of a stage: labels that changed from k-1, and labels that changed and are
# back at their value from k-2 (a 2-cycle in progress).
label_changes <- function(trace, cols) {
  K <- length(trace) - 1L
  changed <- vapply(seq_len(K), function(k)
    sum(trace[[k + 1L]][, cols] != trace[[k]][, cols]), numeric(1))
  flipback <- vapply(seq_len(K), function(k) if (k < 2L) 0 else
    sum(trace[[k + 1L]][, cols] != trace[[k]][, cols] &
        trace[[k + 1L]][, cols] == trace[[k - 1L]][, cols]), numeric(1))
  list(changed = changed, flipback = flipback)
}

# "10 8 2 2 2" -> "10 8 2x3": run-length text, so a 60-iteration sequence fits a table cell.
run_text <- function(x) {
  r <- rle(x)
  paste(ifelse(r$lengths > 1, paste0(r$values, "x", r$lengths), r$values), collapse = " ")
}

experiment_b <- function() {
  rows <- list()
  for (seed in SEEDS_PAIR) {
    fp <- make_factor_pair(seed)
    for (crn in c(TRUE, FALSE)) {
      r <- run_engine(fp$holes, seed, crn, trace_cols = c("fac3", "fac3b"))
      lc <- lapply(r$traces, label_changes, cols = c("fac3", "fac3b"))
      gap <- is.na(fp$holes$fac3)
      rows[[length(rows) + 1L]] <- data.frame(seed, crn, agree = fp$agree,
        iter1 = r$iterations[1], iter2 = r$iterations[2],
        at_cap = r$iterations[2] == MAX_ITER[2], left = r$left,
        changed1 = run_text(lc[[1]]$changed), changed2 = run_text(lc[[2]]$changed),
        flipback = sum(lc[[1]]$flipback) + sum(lc[[2]]$flipback),
        flip_last = utils::tail(lc[[2]]$flipback, 1),
        acc_fac3 = mean(r$out$fac3[gap] == fp$truth$fac3[gap]),
        acc_fac3b = mean(r$out$fac3b[gap] == fp$truth$fac3b[gap]))
      message(format(Sys.time(), "%H:%M:%S"), " (b) seed ", seed, " crn ", crn)
    }
  }
  do.call(rbind, rows)
}

# --- report ------------------------------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
md_table <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"),
                           paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
                           apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))

# Spec 5.5: crn is kept unless it does harm. Harm, as read here (fixed before the full run):
# more total iterations on average, a scaled RMSE worse by more than 2 standard errors of the
# paired difference, or a final stage at its cap in (b) more often than without crn. The
# verdict is reported with the size of each difference; the decision is Carl's.
crn_verdict <- function(a, b) {
  it  <- tapply(a$iter1 + a$iter2, a$crn, mean)
  # Paired by seed explicitly, so the pairing does not depend on the loop's row order.
  m   <- merge(a[a$crn, c("seed", "scaled")], a[!a$crn, c("seed", "scaled")], by = "seed",
               suffixes = c("_crn", "_free"))
  stopifnot(nrow(m) == length(SEEDS))                  # BLOCKING: one pair per seed
  d   <- m$scaled_crn - m$scaled_free
  se  <- stats::sd(d) / sqrt(length(d))
  harm <- c(iterations = it[["TRUE"]] > it[["FALSE"]], rmse = mean(d) > 2 * se,
            cap_b = sum(b$at_cap[b$crn]) > sum(b$at_cap[!b$crn]))
  rel <- mean(d) / mean(a$scaled[!a$crn])
  list(it = it, d = d, se = se, rel = rel, n_worse = sum(d > 0), harm = harm)
}

write_report <- function(a, b, minutes) {
  v <- crn_verdict(a, b)
  ta <- data.frame(seed = a$seed, crn = a$crn, `iter. stage 1` = a$iter1,
                   `iter. stage 2` = a$iter2, `final-stage cells not converged` = a$left,
                   `scaled RMSE` = fmt(a$scaled, 4), `raw RMSE` = fmt(a$raw, 4),
                   seconds = fmt(a$seconds, 0), check.names = FALSE)
  ma <- stats::aggregate(cbind(iter1, iter2, left, scaled, raw, seconds) ~ crn, a, mean)
  tm <- data.frame(crn = ma$crn, `iter. stage 1` = fmt(ma$iter1, 1), `iter. stage 2` = fmt(ma$iter2, 1),
                   `final-stage cells not converged` = fmt(ma$left, 1), `scaled RMSE` = fmt(ma$scaled, 4),
                   `raw RMSE` = fmt(ma$raw, 4), seconds = fmt(ma$seconds, 0), check.names = FALSE)
  tb <- data.frame(seed = b$seed, crn = b$crn, `fac3 = fac3b` = fmt(b$agree, 2),
                   `iter. 1 + 2` = paste(b$iter1, "+", b$iter2), `stage 2 at cap` = b$at_cap,
                   `labels changed per iteration, stage 1` = b$changed1,
                   `stage 2` = b$changed2, `flip-backs` = b$flipback,
                   `flip-backs, last iteration` = b$flip_last,
                   `accuracy fac3 / fac3b` = paste(fmt(b$acc_fac3, 2), "/", fmt(b$acc_fac3b, 2)),
                   check.names = FALSE)
  body <- c(
    "# Experiments: common random numbers and a co-missing factor pair", "",
    "## Summary", "",
    sprintf("- Compared the stage engine with common random numbers (crn, the current default: a column reuses one random stream in every iteration of a stage) and without (a fresh sub-stream per iteration), on numeric data (a) and on two associated factors missing in the same rows (b); %d seeds in (a), %d in (b), tol = (%s), maxIter = (%s), ranger, 2 workers.",
            length(SEEDS), length(SEEDS_PAIR), paste(TOL, collapse = ", "),
            paste(MAX_ITER, collapse = ", ")),
    sprintf("- (a) Mean total iterations: %s with crn, %s without. Paired difference in scaled RMSE (crn minus no crn): mean %s, SD over seeds %s.",
            fmt(v$it[["TRUE"]], 1), fmt(v$it[["FALSE"]], 1), fmt(mean(v$d), 4), fmt(stats::sd(v$d), 4)),
    sprintf("- (b) Final stage at its cap: %d of %d runs with crn, %d of %d without; flip-backs (2-cycle signature) summed over runs: %d with crn, %d without.",
            sum(b$at_cap[b$crn]), sum(b$crn), sum(b$at_cap[!b$crn]), sum(!b$crn),
            sum(b$flipback[b$crn]), sum(b$flipback[!b$crn])),
    sprintf("- Rule verdict (spec 5.5, \"kept unless it does harm\"; harm = more total iterations on average, a scaled RMSE worse by more than 2 standard errors of the paired difference, or the final stage at its cap in (b) more often than without crn): harm found in %s. Size: with crn the scaled RMSE is %s %% higher (worse in %d of %d seeds) at %s times fewer total iterations. The RMSE gap is confounded with iteration count (%s against %s iterations): the extra iterations, rather than the random streams, may explain it (not tested).",
            if (any(v$harm)) paste(names(v$harm)[v$harm], collapse = ", ") else "none",
            fmt(100 * v$rel, 2), v$n_worse, length(v$d), fmt(v$it[["FALSE"]] / v$it[["TRUE"]], 1),
            fmt(v$it[["TRUE"]], 1), fmt(v$it[["FALSE"]], 1)),
    "- **Carl's decision (2026-11-24): keep crn = TRUE.** Uncertainty: 5 seeds in (a), 20 in (b); seconds were measured while other jobs ran on the machine and are inflated.", "",
    "## Definitions", "",
    "- **Scaled RMSE**: root mean squared error over all imputed cells against the simulated truth, each error divided by the SD of its column's observed values (the stopping rule's scale), pooled over columns. **Raw RMSE**: the same without scaling.",
    "- **Labels changed**: imputed labels of the two factors (both columns together) that differ from the previous iteration. **Flip-backs**: labels that changed and now equal their value two iterations earlier; a 2-cycle shows as a steady count of flip-backs and a final stage that runs to its cap. **Flip-backs, last iteration**: flip-backs in the last iteration of stage 2 (non-zero means a cycle still running when the stage stopped). Sequences are run-length coded: `2x13` is 2 in 13 consecutive iterations.", "",
    "## (a) Common random numbers", "",
    "Data: `sim_numeric(200, 300)` (log-normal, rank-3 signal plus noise), 5 % of cells missing completely at random; seed s uses data seed s, hole seed 1000 + s and engine seed s.", "",
    md_table(ta), "", "Means over seeds:", "", md_table(tm), "",
    "Without crn the engine rebuilds each column's sub-stream from the start of the stage in every iteration (cost grows with the square of the iteration count); the seconds columns show what that costs at these sizes.", "",
    "## (b) Two associated factors missing in the same rows", "",
    sprintf("Data: `sim_mixed(300)` plus `fac3b`, a copy of `fac3` whose labels are permuted within a random %d %% of rows; both missing in the same %d random rows (no other gaps). The table after iteration k is recomputed as the engine run with maxIter = k, which is exact (the first k iterations do not depend on maxIter).",
            round(100 * PERMUTE_FRAC), N_COMISS), "",
    "Reading the crn rows: the two factors are each other's only imputed predictors and their training rows hold no imputed cells, so under crn each factor's model is the same in every iteration and an iteration is a fixed, deterministic map. A change therefore shows twice (`Nx2`: a label changes, then its partner follows) and then stops (`0x2`: two unchanged predictions freeze a cell). Seed 6's 2-cycle with crn is a genuine cycle of that map, not Monte Carlo noise.", "",
    md_table(tb), "",
    "## Setup", "",
    sprintf("- %s; mvImputer %s; ranger %s; future %s. One 2-worker multisession pool for all runs. Total run time %s minutes.",
            R.version.string, utils::packageVersion("mvImputer"), utils::packageVersion("ranger"),
            utils::packageVersion("future"), fmt(minutes, 0)),
    "- The engine is driven through its internal stage function so crn can be switched; a check in the script confirms that with crn it returns exactly what `mvImpute()` returns.", "")
  writeLines(body, OUT_PATH, useBytes = TRUE)
}

# --- main ------------------------------------------------------------------------------------
# Modes: no argument runs everything and writes the results (RDS, gitignored) and the
# report; "report" rewrites the report from the saved results without rerunning.
RUNS_PATH <- "bench/output/experiments-runs.rds"
main <- function(args) {
  if (identical(args[1], "report")) {
    r <- readRDS(RUNS_PATH)
    return(write_report(r$a, r$b, r$minutes))
  }
  old_plan <- future::plan(future::multisession, workers = WORKERS)
  on.exit(future::plan(old_plan), add = TRUE)
  t0 <- proc.time()[["elapsed"]]
  check_engine_matches_package()
  a <- experiment_a()
  b <- experiment_b()
  minutes <- (proc.time()[["elapsed"]] - t0) / 60
  saveRDS(list(a = a, b = b, minutes = minutes), RUNS_PATH)
  write_report(a, b, minutes)
}

main(commandArgs(trailingOnly = TRUE))
