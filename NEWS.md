# mvImputer 0.1.0

First release. `mvImpute()` continues the `mvImp()` / `mvImpWrap()` imputer of StatTools
0.0.917 as a standalone package. Everything below is a change against that version.

## Fixes of known problems in StatTools 0.0.917

* **Reproducible seeding (F1).** `mvImpWrap()` started an unseeded parallel cluster, so
  `set.seed()` did not reproduce a run. `mvImpute(seed = )` gives identical imputed values
  for 1, 2 or more workers (same machine, same R and package versions). Each (stage, column)
  has its own random-number stream; the caller's `.Random.seed` and `RNGkind()` are restored
  on exit, except that `seed = NULL` draws exactly one integer from the caller's generator.
* **Stopping rule (F2).** The old rule divided the change by the current value,
  `abs(old - new) / old < tolerance`. A cell moving from -1 to -2 passed as converged (negative
  denominator); a cell at 0 that moved could never converge (infinite ratio), so 0/1
  indicator columns ran to the iteration cap; a ratio exactly equal to the tolerance was
  neither converged nor left; and the rule was not shift-invariant, so centred or log-scale
  data near 0 behaved erratically. The new rule has no division by the current value: a numeric
  cell has converged when its change is below `tol` times the standard deviation of that
  column's observed values; a factor cell when its predicted level has been unchanged for two
  consecutive iterations. "Still iterating" is the exact complement of "converged". `tol` is
  therefore on a different scale from the old `tol1`/`tol2`.
* **Oscillation with one-hot columns (F3).** With full K-column one-hot coding each indicator
  is an exact function of its siblings, and under this update scheme the imputed indicators of
  a missing row swing between about 1 and about 0 on alternate iterations and never settle. In
  a simulation (n = 300, 3-level factor, 45 rows missing) decoded accuracy was 0.58 with
  siblings as predictors against 0.82 without. Factors are now handled directly: each factor
  is one imputation target and its own indicator columns are never its predictors. Analyses
  that one-hot encoded factors and imputed them with the old engine probably got near-random
  factor imputations and should be re-run. Exact linear dependencies among numeric columns that
  the user supplies still can cause the oscillation; this is documented under "Assumptions".
* **Housekeeping (F4).** Progress goes through `message()` (silenceable), not `cat()`;
  `options(warn = -1)` is no longer set globally; there is no `library()` call inside
  functions; all dependencies are declared in `DESCRIPTION`; ranger runs with one thread per
  worker instead of its default two.

## Training rows

A column's model is now trained only on the rows where that column was originally observed.
StatTools 0.0.917 trained on every row not still iterating, so converged imputed cells became
training data for their own column. Imputed values therefore differ from the old engine's even
where nothing else changed. (Decided 2026-11-23: pooled imputation, no `train` argument.)

## New features

* Data frames with numeric, integer, factor, ordered factor, character and logical columns;
  the same class and column types come back, including unused factor levels. Matrices are
  still accepted.
* `method = "ranger"` (default), `"randomForest"` or `"pls"`, replacing `method` plus `rfMeth`;
  the Rborist backend is dropped.
* `workers` replaces `nCore`: parallelism through `future` / `future.apply`; an active
  `future::plan()` is used and left in place.
* `exclude`: columns passed through untouched (identifiers, the outcome); never predictors,
  never imputed, gaps kept.
* `attr(out, "mvImpute")`: seed, settings, iterations and time per stage, the cells not
  converged per stage, an out-of-bag error per imputed column, and rows with no observed value.
* Input checks with explicit messages (infinite values, identifier-like columns, columns with
  fewer than two observed values, user-supplied one-hot blocks); constant columns are filled
  with their value and take no part in the iteration.
* A warning when the final stage ends with cells not converged.

## Speed against StatTools 0.0.917

Simulated data, same number of workers, one run each (`bench/output/benchmark.md` in the
repository). The old engine ran two ranger threads in each worker and the new one runs one.

* Numeric tables (200 x 1000, 5 % gaps): wall time new/old 0.36 to 0.54, because the new
  engine stopped after 11 to 14 iterations while the old one ran to its cap of 75.
* Small mixed table (900 x 20, one-hot route for the old engine): 1.53 to 1.92, both engines at
  75 iterations. Higher cost per iteration, possibly from single-thread forests and dispatch
  overhead (a hypothesis).

## Default settings (decided 2026-11-24, from simulation)

* `tol = c(0.05, 0.025)` stays the default. Over a grid of 12 tolerance pairs on simulated
  omics-like data (random and left-censored gaps) and a mixed covariate table, imputation
  accuracy was flat (largest mean difference 0.003 in SD-scaled RMSE), so no other setting
  was better supported. StatTools 0.0.917 ran to its iteration caps in every run, so
  its iteration count did not narrow the choice.
* Common random numbers stay on: within a stage, each column reuses one random stream in
  every iteration. Compared with fresh random numbers per iteration, this gave a 0.4 %
  higher SD-scaled RMSE with about 6 times fewer iterations, and fewer label cycles in a test
  with two associated factors missing in the same rows (final stage at its cap in 1 of 20
  runs, against 6 of 20).

## Known limitations

* On tables with factor-like columns the final stage often reaches `maxIter`, because
  predicted factor levels can cycle between iterations, so the "not converged" warning is
  common there (54 of 60 simulated runs on a mixed table). The cells concerned are listed in
  `attr(out, "mvImpute")$not_converged`. To be revisited after first real use.
* On simulated left-censored data (200 x 500, 5 % of the lowest values missing) the raw RMSE
  of `mvImpute()` with the default `tol` was about 3 % higher than that of StatTools 0.0.917
  (mean 0.593 against 0.578, which is 2.6 %; higher in all 5 replicates, by 0.9 % to 4.2 %),
  while the SD-scaled RMSE was equal (0.494 against 0.492). In a decomposition on one
  200 x 1000 data set the gap came from a few cells in the highest-SD, heavy-tailed columns,
  and a smaller `tol` did not close it (0.595 at `tol = c(0.025, 0.005)`). Its cause is not
  established.
