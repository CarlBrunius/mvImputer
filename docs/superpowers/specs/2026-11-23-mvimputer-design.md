# mvImputer: design

Date: 2026-11-23. Status: **approved by Carl, 2026-11-23** (after two independent reviews, both
folded in). Requirements decided by Carl on 2026-11-18; design decisions from the
brainstorming session on 2026-11-23.

## 1. Intent and why

`mvImp()` and `mvImpWrap()` (StatTools, GitLab `CarlBrunius/StatTools`, version 0.0.917, commit
`7ee1dc09`) are Carl's two-stage iterative imputer for omics and covariate tables: each variable with
gaps is predicted from all the others by a random forest or PLS, iterating until the imputed values
stop changing. The functions are in regular use in the author's research projects. They move into
a new standalone, public R package, **mvImputer** (this repository), with:

1. reproducible seeding: the same seed gives identical imputations whatever the number of workers;
2. parallelism through `future` / `future.apply`, replacing foreach/doParallel;
3. automatic handling of factor columns (one-hot encoding as predictors);
4. factors returned on their original scale as the most likely category;
5. a stopping rule that works for zeros, negative values, centred data and indicator columns.

**Why it matters (pre-mortem, Carl):** the function will keep being used only if the new version is
as easy to call and about as fast as the old one. Seeding, better parallelism and factor support are
the gains; clunkiness or a large slowdown are the ways it fails.

## 2. What is wrong today (evidence)

All demonstrations used simulated data; results are quoted here.

- **F1, not reproducible.** `mvImpWrap` registers a plain doParallel cluster without doRNG, so
  `set.seed()` does not reproduce it. A project-level workaround (doRNG) reproduces only for a
  fixed worker count, or serially.
- **F2, stopping rule** (`mvImp`, the two lines computing `whichConverged` and `stillLeft`):
  `abs(old - new) / old < tolerance`.
  - Negative denominator: a cell moving from −1 to −2 (100 % change) gives ratio −1 and passes as
    converged. Confirmed.
  - Zero denominator: a cell at 0 that moves gets ratio `Inf` and can never converge; 0/1 indicator
    columns therefore run to the iteration cap (seen in practice).
  - The two conditions are not complements: a ratio exactly equal to `tolerance` is in neither set.
  - Relative change is not shift-invariant: centred or log-scale data near 0 behave erratically.
- **F3, sibling oscillation with one-hot columns.** With full K-column coding (as
  `stats::model.matrix(~ f - 1)` produces), each indicator is an exact function of its siblings in observed rows (`a = 1 − b − c`).
  Under mvImp's update scheme (every feature predicted from the previous iteration's values) the
  imputed indicators of a missing row swing between about 1 and about 0 on alternate iterations and
  never settle. Simulation (n = 300, 3-level factor, 45 rows missing, ranger as in mvImp): decoded
  accuracy 0.58 with siblings as predictors versus 0.82 without them. **Consequence outside this
  package:** analyses that one-hot encoded factors and imputed them with the old engine probably
  got near-random factor imputations and should be re-run once mvImputer exists.
- **F4, minor:** messages via `cat()` cannot be silenced; `options(warn = -1)` is set globally during
  the run; `library(doParallel)` is called inside functions; doParallel and foreach are not declared
  in StatTools' DESCRIPTION; ranger runs with its default 2 threads inside every worker.

## 3. Goals, non-goals, anti-goals

**Goals**
- One exported imputation function with sensible defaults: `mvImpute(X)` must do the right thing.
- Bit-identical imputed data for a given seed across 1, 2 or more workers and across `future`
  backends, on the same machine with the same R and package versions (BLAS threading may differ
  between machines; not promised there).
- Mixed data frames: numeric, integer, factor, ordered factor, logical, character columns.
- A stopping rule without division by the current value.
- CRAN-ready from the start (`R CMD check --as-cran` clean; examples and tests use at most 2 cores);
  submission is decided later.
- Run time comparable to `StatTools::mvImpWrap` at the same worker count (graded, see section 9).

**Non-goals (decided)**
- Imputation anchored to training rows (a `train` argument that keeps test rows out of the models).
  Carl, 2026-11-23: pooled imputation only.

**Assumptions the help page must state** (estimand and missingness)
- Pooled: every row, including rows a later cross-validation fold holds out, informs every imputed
  value. Downstream CV or test performance is therefore optimistic, more so if the outcome is a
  column of `X`; the help page advises leaving the outcome out of `X`.
- Missingness: the models assume values are missing at random given the other columns. Values
  missing because they are below a detection limit (missing not at random) are imputed too high:
  forests cannot predict below the observed range.
- Single imputation: one deterministic completed dataset; downstream standard errors that treat
  imputed values as observed are too small.
- Exact linear dependencies among numeric columns the user supplies (their own one-hot columns, parts
  and a total) can make the iterations oscillate (finding F3); factor columns are handled safely.
- Storing fitted models to impute new data later.
- The Rborist backend (dropped).
- A user-extensible backend framework.

**Anti-goals** (a failure even if it works)
- Noticeably clunkier to call than `mvImpWrap(MAT)`.
- Output that silently differs from input outside the missing cells (observed values must come back
  bit-identical).
- Seeding that touches the caller's global RNG state beyond one documented draw.

**Wish list** (not in this round)
- Evaluate the `mtry` default: balance speed and accuracy; 5 (the old value) may not be optimal
  (Carl, 2026-11-23). Evaluate further forest settings for speed (node size, depth, row sampling
  fraction) in the same study.
- Predictive mean matching (as missRanger's `pmm.k`): imputations become observed donor values.
- Multiple imputation (m seeds); needs predictive mean matching first, or the between-imputation
  variance is only forest noise.
- An exported evaluation helper: mask observed cells, impute, report error (as missForest's
  `prodNA`/`mixError`, missRanger's `generateNA`); built internally for tests and benchmark now.

Prior art reviewed 2026-11-23: missForest 1.6.1 and missRanger 2.6.1 (CRAN sources). Rejected from
them: out-of-sample prediction with stored forests (non-goal above), case weights (randomForest and
PLS lack them), class weights and cutoffs, sequential variable updates (would give up parallelism).

## 4. Approaches considered

- **A. Port and restructure mvImp (chosen).** Keep the algorithm (two stages, all features updated
  from the previous iteration so they can run in parallel, cells freeze once converged); split it
  into small internal steps (prepare, iterate, fit-and-predict per backend, decode); backends in an
  internal lookup list.
- **B. Build on missRanger (CRAN).** Mature, handles factors and seeds, but a different algorithm
  (sequential updates, no parallelism across features, no PLS, global stopping rule). Rejected; may
  serve as an accuracy reference in the benchmark.
- **C. User-pluggable backend framework (exported S3 methods).** Speculative; A keeps the option open.

On StatTools (separate repository): Carl chose to **freeze the old functions as deprecated**,
following the precedent of `rfImp`/`rfImpWrap`: code unchanged, `.Deprecated()` pointing to
mvImputer, README naming the known stopping-rule bugs, removal in a later StatTools release. Tag the
last StatTools commit before the change. Alternatives rejected: deleting outright (breaks callers),
forwarding stubs (old scripts would silently get new numbers).

## 5. Design

### 5.1 Interface (clean break, new names)

```r
mvImpute(
  X,                          # data.frame or matrix; the same class is returned
  method    = c("ranger", "randomForest", "pls"),
  exclude   = NULL,           # names of columns passed through untouched (IDs, outcome)
  seed      = NULL,
  workers   = NULL,
  tol       = c(0.05, 0.025), # one per stage; fraction of the feature's SD; defaults to be recalibrated (5.4)
  maxIter   = c(15, 60),      # one per stage; length(tol) == length(maxIter) = number of stages
  guess     = "minVar",       # "minVar" | "median" | "minTotal" | matrix/data.frame of X's shape
  forceZero = FALSE,          # floor numeric imputations at 0
  ntree     = 100,
  mtry      = 5,              # kept from mvImp for speed; capped at (number of predictors); see wish list
  nComp     = 2,              # PLS components
  verbose   = FALSE
)
```

- `method`: replaces the old `method` + `rfMeth` pair.
- `exclude` (Carl, 2026-11-23, after the prior-art review): a character vector of column names.
  These columns are returned unchanged (gaps stay gaps), are not imputed and are never predictors.
  All input checks except "unknown column name" skip them (so an ID column can simply be excluded).
  Fewer than 2 remaining columns is an error. The help page recommends excluding the outcome.
- `seed`: an integer. `NULL` draws one integer from the caller's RNG, so `set.seed(1); mvImpute(X)`
  reproduces; this is the only touch of the caller's RNG state. The seed used is recorded in the
  output.
- `workers`:
  - `NULL`: if the active `future::plan()` is not sequential, use it; otherwise run a multisession
    plan with `min(8, future::availableCores() - 1)` workers (at least 1) for this call and restore
    the previous plan on exit.
  - `1`: serial (`sequential` plan for the call).
  - `k > 1`: a multisession plan with `k` workers for the call.
  - `future` cannot tell "never set" from "explicitly sequential"; users who want serial pass
    `workers = 1`. Multisession on every OS (forking is unavailable on Windows and unsafe in RStudio).
- `guess`: read only at missing cells (observed values always come from `X`); a user matrix or
  data frame must be non-missing there. Factor columns always start at their most frequent observed
  level (ties: first level). `"minVar"` and `"median"` are per numeric column; `"minTotal"` is the
  minimum over all numeric columns' observed values. `"minVar"` stays the default (Carl,
  2026-11-23: omics is the main use); the help page advises `guess = "median"` for covariate tables.
- `verbose = TRUE` reports progress per stage and iteration through `message()`. No `cat()`.
- The **final** stage reaching `maxIter` with unconverged cells gives a `warning()` naming the count;
  earlier stages (the "sloppy" stage 1 by design) only record it in the diagnostics. *(Builder's
  call, from review; Carl may change.)*

**Output:** the same class as `X`. Column types are restored: factors with their original levels
(including unused levels), ordered factors stay ordered, logical comes back logical, character comes
back character. An integer column without gaps stays integer; an integer column with gaps comes back
as double (unrounded), as documented. Row and column names are kept. Diagnostics
go in `attr(out, "mvImpute")`, a list with: `seed`, `rng_kind`, `method`, the settings used,
`iterations` per stage, `not_converged` (row, column) per stage, `time` per stage, and `oob_error`.

`oob_error` (Carl, 2026-11-23, after the prior-art review): one value per imputed variable, from
the last fit of that variable in the final stage, estimating how well it can be predicted from the
others. Numeric: out-of-bag RMSE divided by `s_j` (the stopping-rule scale). Factor: out-of-bag
misclassification rate. `NA` for PLS. Computed from the out-of-bag predictions the fitted forest
stores, not from ranger's `prediction.error` (hypothesis to verify in the plan: for probability
forests that is a Brier score, not a misclassification rate). The help page notes that it is
estimated on observed rows, so under missing-not-at-random gaps it overstates imputation quality.

**Input checks (errors, with the reason):** not a data.frame or matrix; unsupported column type
(e.g. Date, list); `Inf`/`-Inf` present ("set them to NA if they mean missing"; the old code silently
treated them as missing; `NaN` counts as missing, as `is.na()` does); fewer than 2 columns; a column
with gaps and fewer than 2 observed values; `length(tol) != length(maxIter)`; any `tol <= 0` or
`maxIter < 1`; a character or factor column whose observed values are all distinct (an identifier,
which would become one class per row); `method = "pls"` without
mixOmics installed. A column whose observed values are constant (a numeric with one distinct value, or a
factor with one observed level) is filled with that constant, with a `message()`, and takes no part
in iteration as a target. Rows with no observed value give a `warning()` (they are imputed from
first guesses and other imputations only) and are listed in the diagnostics. Numeric 0/1 columns
that sum to 1 in every complete row (a user-supplied one-hot block) give a `warning()` pointing to
F3 and suggesting a factor column instead. A matrix input must be numeric. No missing values: `X`
is returned unchanged (with the diagnostics attribute) and a message.

### 5.2 Internal representation

- **Numeric** (double, integer): one column, as is.
- **Nominal factor** (factor, character, logical): one-hot indicator columns, used **only as
  predictors** for other variables. The factor itself is one imputation target.
- **Ordered factor** (hybrid, Carl 2026-11-23): one integer column 1..K, used only as a predictor;
  the factor itself is one imputation target, imputed like a nominal factor.
- During iteration a missing factor cell enters other variables' models through its current most
  likely level (hard 0/1 indicators, or the integer code for an ordered factor), never as
  probabilities.
- A factor's own indicator (or integer) columns are excluded from the predictors of that factor.
  This removes the F3 oscillation by construction.

### 5.3 Models per target type

| Target | ranger | randomForest | pls |
|---|---|---|---|
| numeric | regression forest | regression forest | PLS regression, `nComp` components |
| factor (nominal or ordered) | probability forest; most likely class | classification forest; class with the most votes (the probability analogue for randomForest) | PLS2 on the factor's indicator block (PLS-DA); level with the largest predicted score |

- ranger runs with `num.threads = 1` inside workers (parallelism comes from `future`). Verified:
  ranger 0.18.0 with a fixed seed gives identical predictions for 1, 2, 4 and 8 threads.
- **Per fit:** predictors with zero variance in the training rows are dropped (verified: mixOmics
  6.36 PLS warns "the standard deviation is zero" on a constant column; randomForest stops on empty
  classes); factor targets are fitted on `droplevels()` of the training labels and mapped back to
  the full level set; `nComp` is capped at the number of remaining predictors; any non-finite
  prediction is an error, never silently "converged" or "left".
- PLS: predictors and the indicator response block are centred and scaled (mixOmics default
  `scale = TRUE`).
- For PLS, "most likely" means the largest predicted indicator score, not a probability; the help
  page says so.
- `forceZero` applies to numeric targets only.

### 5.4 Stopping rule (requirement 5)

Per cell, per stage, over the cells still iterating:

- **Numeric cell:** converged when `abs(new - old) / s_j < tol[stage]`, where `s_j` is the SD of the
  **observed** values of column j, computed once before stage 1. `s_j > 0` is guaranteed because
  constant columns are removed beforehand (5.1). No division by the current value, so zeros,
  negative values, centred data and log-scale data are handled alike.
- **Factor cell:** converged when its predicted level is unchanged for **two consecutive**
  iterations, so a cell whose first prediction happens to equal the modal first guess does not
  freeze on a model fitted to first-guess predictors. *(Builder's call, from review.)*
- **Still left** is the exact complement of **converged**. The rule is one internal pure function
  (old values, new values, scales, types, tol → logical), tested directly, including the
  exactly-at-`tol` case that forest output cannot reach on purpose.
- Converged cells freeze (as in mvImp). A feature is refitted in an iteration while any of its cells
  is still iterating; its predictions are written only to its unfrozen cells.
- **Training rows** (Carl, 2026-11-23; a deliberate change from mvImp): a feature's model is trained
  only on rows where that feature was **originally observed**. Frozen imputed cells are never
  training targets; they serve only as predictors for other features. (Old mvImp trained on every
  row not still iterating, so converged imputations became training data for their own feature.)
- A stage ends when no cells are left or at `maxIter[stage]`. Stage 2 starts from stage 1's result
  with all missing cells unfrozen.
- **Recalibration:** `tol` now means "fraction of the feature's SD". The defaults `c(0.05, 0.025)`
  are placeholders. The plan includes a simulation that picks defaults: the smallest-RMSE setting
  whose total iteration count does not exceed the old engine's, on (i) omics-like positive data with
  random (MCAR) gaps, (ii) the same with left-censored gaps, (iii) a mixed covariate table; reported
  per scenario. Carl signs off the chosen values before first real use. *(Builder's call on the
  rule; Carl may change.)*
- **Known trade-off (accepted):** for right-skewed intensities the SD is driven by high values, so
  low imputed values (near the detection limit, typically missing-not-at-random) are judged more
  loosely than under the old relative rule.

### 5.5 Seeding

- From `seed`, generate one L'Ecuyer-CMRG stream per **(stage, original column)** with
  `parallel::nextRNGStream`. The caller's RNG kind and `.Random.seed` are saved on entry and
  restored on exit (`on.exit`), including when `.Random.seed` did not exist. This is required, not
  optional: verified that `future_lapply` with a list of seeds still advances the caller's
  `.Random.seed` (future.apply 1.20.2). With `seed = NULL` the state on exit equals the state after
  exactly one draw.
- Each iteration passes the streams of the features being fitted to `future.apply::future_lapply`
  as `future.seed = <list of seeds>`. Verified: a list of pre-generated seeds gives identical results
  sequentially and with 3 multisession workers (future 1.70.0, future.apply 1.20.2).
- **Common random numbers:** a feature reuses the same stream in every iteration of a stage, so a
  change between iterations reflects the data, not new random draws. This keeps forest Monte Carlo
  noise from blocking convergence at small `tol`. The plan's first experiment measures iteration
  counts with and without it; it is kept unless it does harm (it costs nothing). The same experiment
  includes two strongly associated factors missing in the same rows, to check for a 2-cycle under
  hard-level predictors.
- Streams are keyed by the column's position among the non-excluded columns, not by position in
  the list of features with gaps, so the set of features still iterating cannot shift another
  feature's random numbers, and adding or removing excluded columns changes nothing. *(Amended
  at plan review, 2026-11-23; was "original column index".)*
- ranger gets `seed = NULL` (it then draws its seed from R's RNG in the worker, per its help page);
  randomForest uses R's RNG directly; PLS is deterministic.

### 5.6 Package

- R package `mvImputer`, roxygen2, testthat 3rd edition. Exports: `mvImpute()` only.
- Imports: `future`, `future.apply`, `ranger`, `randomForest`, `stats`, `parallel`, `withr`
  (Carl, 2026-11-23: randomForest required). Suggests: `mixOmics` (Bioconductor; PLS only),
  `testthat`. The benchmark against StatTools lives in a `bench/` folder listed in `.Rbuildignore`,
  so StatTools (GitLab only) never enters DESCRIPTION.
- Namespace-qualified calls everywhere; no `library()` in package code; no global options changed.
- Public GitHub repo under `CarlBrunius`; commits as `carl.brunius@gmail.com`. Nothing committed or
  pushed until Carl says so.
- Licence: **MIT** (`License: MIT + file LICENSE`), Carl 2026-11-23. Maintainer (`cre`) e-mail:
  `carl.brunius@gmail.com` (CRAN allows one); Carl's NBIS address appears in the README as an
  additional contact. `BugReports:` points to the GitHub issues page.

## 6. Tests (testthat; simulated data only)

- **Seed reproducibility** (identity on the imputed data with the diagnostics attribute stripped,
  plus the attribute's `seed` and `iterations`; `time` and `workers` differ legitimately): workers 1
  and 2 (CRAN limit), plus 3+ in a test skipped on CRAN; same seed twice; different seeds differ;
  `seed = NULL` after `set.seed()` reproduces; caller's `.Random.seed` and RNG kind unchanged with a
  seed given, advanced by exactly one draw with `seed = NULL`, and still absent if it was absent.
- **Observed values untouched:** bit-identical values for every column type (integer columns with
  gaps compared by value, since they return as double).
- **Factor round trips:** levels, unused levels, ordering, character and logical types restored; a
  missing factor gets a valid level; no indicator oscillation (iterations below the cap on a
  simulated factor with informative predictors).
- **Stopping-rule edge cases:** zeros, negative values, centred data, a constant column, a ratio
  exactly at `tol`, a column that is all zeros except one value.
- **Input checks:** each error in 5.1 fires with its message.
- **`exclude`:** excluded columns come back identical (including their gaps); changing an excluded
  column's values does not change any imputation; an excluded ID column avoids the identifier error.
- **`oob_error`:** present for every imputed variable with ranger and randomForest, `NA` with pls;
  near 0 for a variable that is an exact function of another, near 1 (numeric) for pure noise.
- **Backends:** each of ranger, randomForest, pls on numeric and factor targets (pls skipped without
  mixOmics).
- **Accuracy sanity:** on simulated data with known truth, imputed values beat column-mean
  imputation (RMSE for numeric, accuracy for factors).

## 7. Decisions from the interview (2026-11-23)

| Question | Answer | Class |
|---|---|---|
| API | Clean break; new function names | build now |
| Parallel control | Hybrid: use an active plan, else multisession with `min(8, cores - 1)` | build now |
| Backends | ranger + randomForest + pls; Rborist dropped; randomForest required | build now |
| Factor imputation | One model per factor (option b); one-hot only as predictors | build now |
| Ordered factors | Hybrid: integer code as predictor, imputed as a class | build now |
| Stopping rule | SD-scaled for numeric, class stability for factors, common random numbers | build now |
| Training-row anchoring | Not wanted; pooled only | non-goal |
| CRAN | CRAN-ready; submit later | build now |
| First real use | Within 2 weeks | review tiers |
| StatTools | Freeze as deprecated, README pointer, tag | separate repo task |
| `mtry` | Keep 5 now; evaluate later | wish list |
| Training rows after freezing | Originally observed rows only (not old mvImp behaviour) | build now |
| Default first guess | `"minVar"`; help advises `"median"` for covariates | build now |
| From missForest/missRanger | `oob_error` diagnostic and `exclude` now; PMM, multiple imputation, forest speed settings, evaluation helper later | build now / wish list |

**Review tiers** (first use within 2 weeks): *full* review for seeding, the stopping rule, factor
encoding and decoding, and observed-value invariance; *light* for docs, the benchmark script and
the StatTools deprecation.

## 8. StatTools follow-up (separate repository, GitLab `CarlBrunius/StatTools`)

Add `.Deprecated("mvImputer::mvImpute")` to `mvImp` and `mvImpWrap`; README rows marked deprecated
with a pointer to mvImputer and one line on the stopping-rule bugs; NEWS/version bump. Tag the last
pre-change commit. Not committed or pushed until Carl says so; it ships after mvImputer is public,
so the pointer resolves.

Other follow-up (outside this work): analyses that used the one-hot route with the old engine
should be re-run with mvImputer (finding F3).

## 9. Acceptance

- All tests in section 6 pass; `R CMD check --as-cran` with no errors or warnings.
- **Benchmark** (simulated, omics-sized, e.g. 200 × 1000 with 5 % gaps, under MCAR and left-censored
  gaps reported separately, plus a 900 × 20 mixed covariate table): run time of `mvImpute` relative
  to `StatTools::mvImpWrap` at the same worker count, and imputation accuracy against the simulated
  truth. For the mixed table the comparator is one-hot encoding plus the old engine (the route
  earlier analyses used). The old engine ran ranger with 2 threads per worker, the new one with 1; recorded next to
  the ratio. The time ratio is reported, not pass/fail: Carl weighs it against the added
  functionality; anything above about 1.25× is flagged.
- `tol` defaults chosen from the recalibration simulation (5.4), with the choice written up.

## 10. Left to the builder

- Internal function boundaries and file layout.
- The exact RNG isolation mechanism (5.5), subject to its test.
- How the benchmark compares against StatTools (installed package or the vendored copy).
- Message wording.
