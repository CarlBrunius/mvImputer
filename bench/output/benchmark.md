# Benchmark: mvImpute() against the old engine

## Summary

- Ratio of wall time new/old, one run each: (i, 2 workers) 0.51, new 167 s, old 328 s; (i, 4 workers) 0.54, new 104 s, old 190 s; (ii, 2 workers) 0.36, new 134 s, old 377 s; (ii, 4 workers) 0.39, new 82 s, old 210 s; (iii, 2 workers) **FLAG 1.53**, new 36 s, old 23 s; (iii, 4 workers) **FLAG 1.92**, new 33 s, old 17 s.
- **FLAG: scenario iii at 2 workers, ratio 1.53; scenario iii at 4 workers, ratio 1.92 exceed(s) 1.25.** A ratio compares whole runs, not cost per iteration: the old engine ran (i, 2w) 75, (i, 4w) 75, (ii, 2w) 75, (ii, 4w) 75, (iii, 2w) 75, (iii, 4w) 75 iterations (caps 15 + 60 = 75), the new engine (i, 2w) 14, (i, 4w) 14, (ii, 2w) 11, (ii, 4w) 11, (iii, 2w) 75, (iii, 4w) 75 (it stops a stage when its cells have converged).
- Time per iteration (wall seconds divided by iterations, averaged over the run), new/old, in scenarios i and ii: 2.4 to 2.9. The new engine is faster there only because it stops early.
- This is a single run per cell (data seed 1). The old engine's runs are not reproducible (its parallel workers are not seeded) and timings carry machine noise, so differences of a few percent in time mean nothing.
- Raw RMSE: in scenarios i and ii the new engine's raw RMSE is 7 % to 15 % above the old engine's, while the scaled RMSE differs by at most 0.5 %. Raw RMSE weights the highest-SD, heavy-tailed columns heavily, and a few imputed cells there carried most of the scenario ii gap in a decomposition of that data set. Single-run raw RMSE is also unstable: one repeat run of the old engine in scenario i (2 workers) moved its raw RMSE from 1.98 to 2.18. Over 5 paired replicates (`bench/output/calibration.md`) the raw RMSE was: random gaps 2.868 (new) against 2.880 (old), left-censored 0.593 against 0.578; the cause of the left-censored gap is not established.
- Replicate-level accuracy (5 data sets per scenario, paired) is in `bench/output/calibration.md`.

## Definitions

- **Wall seconds**: `proc.time()` elapsed around the engine call. Both engines start their own worker pool inside the call, so pool start-up is included for both; data building, decoding and scoring are not.
- **Ratio new/old**: wall seconds of `mvImpute()` divided by wall seconds of `StatTools::mvImpWrap()` at the same scenario and worker count. Above 1.25 it is flagged.
- **Scaled RMSE**: RMSE over all imputed numeric cells against the truth, each error divided by the SD of that column's observed values, pooled over columns (constant columns left out). **Raw RMSE** is the same without scaling. **Factor accuracy** (mixed table) is the share of imputed factor-like cells (factor, ordered, character, logical) with the true label.
- **Load before / after**: total CPU load (%, mean of 5 one-second samples) just before and just after the engine call. A run starts only when the load before is below 50 %.
- **Cells left**: imputed cells not converged at the end of the last stage (old engine: from its log; for the mixed table these are cells of its one-hot indicator matrix, not table cells, so not comparable with the new engine's).

## Setup

- Scenarios: (i) `sim_numeric(200, 1000)` with 5 % of cells missing completely at random; (ii) the same data with the lowest 5 % of each column missing (left-censored); (iii) `sim_mixed(900, extra = 10)` (20 columns of numeric, integer, factor, ordered, character, logical, a constant and 10 noise columns) with 5 % MCAR over all columns. Data seed 1, hole seed 1001, `seed = 1` for `mvImpute()`; the old engine after a scoped `set.seed(1)` (no effect on its workers).
- `mvImpute()`: its final defaults (ranger, 100 trees, mtry 5, guess "minVar", `tol = c(0.05, 0.025)`, `maxIter = c(15, 60)`, common random numbers on); only `seed` and `workers` are passed.
- Old engine: StatTools 0.0.917, commit `7ee1dc09bda153de7b62a344484db53b09613516`; `mvImpWrap(method = "RF", rfMeth = "ranger", tol1 = 0.05, n1 = 15, tol2 = 0.025, n2 = 60, nCore = <workers>)`. Mixed table through the one-hot route: character and logical columns made factors, each factor (ordered included) turned into K indicator columns, decoded by the largest indicator (ties to the first level).
- **Thread note**: the old engine calls ranger without `num.threads`, and ranger's default is 2 threads (no `R_RANGER_NUM_THREADS`, `ranger.num.threads` or `Ncpus` was set; checked at each run), so it runs 2 ranger threads in each of its workers: 2 workers use up to 4 threads, 4 workers up to 8. `mvImpute()` fits each forest with 1 thread, so the same worker count gives it half the threads. The comparison is at equal worker count, not equal thread count.
- **Iterations**: the old engine ran to both caps (75 iterations) in the runs in `calibration.md`, as it does here wherever its count is 75; in the mixed table that follows from its stopping rule on indicator columns and is not a meaningful reference cost (see the note in `calibration.md`). The new engine stops when converged.
- **Timing conditions**: one R process per run, runs one at a time; each run waited for total CPU load below 50 % before starting and ran at BelowNormal priority (own process tree only, children inherit), as did the other engine; the engine that goes first alternates between scenarios (the 'Goes first' column). The machine is shared, so the load columns are part of the result.
- Outside reference: `missRanger` is not installed here, so there is no missRanger row (packages are never installed by this benchmark).
- R version 4.6.1 (2026-06-24 ucrt); mvImputer 0.1.0; ranger 0.18.0; future 1.70.0; 24 logical cores.

## Scenario i: numeric 200 x 1000, 5 % MCAR

| Workers | Engine | Goes first | Load before / after (%) | Wall seconds | Ratio new/old | Iter. stage 1 + 2 | Cells left | Scaled RMSE | Raw RMSE |
|---|---|---|---|---|---|---|---|---|---|
| 2 | new | yes | 1 / 6 | 167 | 0.51 | 6 + 8 | 0 | 0.774 | 2.280 |
| 2 | old |  | 6 / 5 | 328 |  | 15 + 60 | 8 | 0.771 | 1.980 |
| 4 | new | yes | 7 / 6 | 104 | 0.54 | 6 + 8 | 0 | 0.774 | 2.280 |
| 4 | old |  | 4 / 4 | 190 |  | 15 + 60 | 9 | 0.770 | 2.014 |

## Scenario ii: numeric 200 x 1000, 5 % left-censored

| Workers | Engine | Goes first | Load before / after (%) | Wall seconds | Ratio new/old | Iter. stage 1 + 2 | Cells left | Scaled RMSE | Raw RMSE |
|---|---|---|---|---|---|---|---|---|---|
| 2 | new |  | 3 / 6 | 134 | 0.36 | 5 + 6 | 0 | 0.472 | 0.632 |
| 2 | old | yes | 4 / 7 | 377 |  | 15 + 60 | 37 | 0.471 | 0.592 |
| 4 | new |  | 4 / 5 | 82 | 0.39 | 5 + 6 | 0 | 0.472 | 0.632 |
| 4 | old | yes | 3 / 3 | 210 |  | 15 + 60 | 20 | 0.471 | 0.580 |

## Scenario iii: mixed 900 x 20, 5 % MCAR

| Workers | Engine | Goes first | Load before / after (%) | Wall seconds | Ratio new/old | Iter. stage 1 + 2 | Cells left | Scaled RMSE | Raw RMSE | Factor accuracy |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | new | yes | 3 / 7 | 36 | **FLAG 1.53** | 15 + 60 | 11 | 0.778 | 0.891 | 0.877 |
| 2 | old |  | 4 / 5 | 23 |  | 15 + 60 | 256 | 0.778 | 0.891 | 0.882 |
| 4 | new | yes | 4 / 4 | 33 | **FLAG 1.92** | 15 + 60 | 11 | 0.778 | 0.891 | 0.877 |
| 4 | old |  | 5 / 7 | 17 |  | 15 + 60 | 255 | 0.778 | 0.892 | 0.845 |

