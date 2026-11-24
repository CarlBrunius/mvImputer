# mvImputer

Multivariate imputation of missing values in numeric and mixed-type tables. Each column with
gaps is predicted from all the other columns by a random forest (or PLS), the predictions
replace the gaps, and this is repeated until the imputed values stop moving. Numeric, integer,
factor, ordered factor, character and logical columns are supported; factors are imputed as
one target each and come back as the most likely level.

The package exports one function, `mvImpute()`. It continues the `mvImp()` / `mvImpWrap()`
imputer of the StatTools package (version 0.0.917) with a reproducible seed, parallelism
through `future`, factor support and a corrected stopping rule (see `NEWS.md`).

## Install

```r
remotes::install_github("CarlBrunius/mvImputer")
```

## Example

```r
library(mvImputer)

set.seed(1)
d <- data.frame(a = rnorm(60), b = rnorm(60))
d$c <- d$a + d$b + rnorm(60, sd = 0.1)
d$g <- factor(ifelse(d$a > 0, "up", "down"))
d$a[c(3, 9)] <- NA
d$g[c(5, 20)] <- NA

out <- mvImpute(d, seed = 1, workers = 1)
attr(out, "mvImpute")$iterations   # iterations run in each stage
```

* The same seed gives identical imputed values whatever `workers` is (same machine, same R and
  package versions). The caller's random-number state is restored on exit, except that
  `seed = NULL` (the default) takes one draw from the caller's random-number stream.
* `workers = NULL` (default) uses an active non-sequential `future::plan()` if there is one,
  otherwise starts a multisession plan for the call. `workers = 1` runs serially.
* Leave the outcome out of `X`, or pass it in `exclude = "outcome"`: excluded columns are
  returned unchanged and never used as predictors.
* `guess` sets only the starting values (default `"minVar"`, the column minimum; `"median"`
  is the natural start for a covariate table). The final values hardly depend on it.
* Diagnostics (seed, iterations per stage, cells that did not converge, out-of-bag error per
  column) are in `attr(out, "mvImpute")`. See `?mvImpute`.

## Choosing settings for your data

Each recommendation is marked **Simulation** when simulations in the repository support it
(`bench/output/detection-limit.md`, `bench/output/calibration.md`) and **Reasoning** when it has
not been evaluated here. The simulations use simulated, low-rank, log-normal data, so they say
what happened there, not what will happen on your data.

**General rules**

* **Impute on the scale you will analyse.** For right-skewed intensities (metabolomics, MS
  proteomics, SomaScan, microarray), log-transform first and back-transform afterwards if you
  need original units. *Simulation:* on simulated data with the lowest values missing, imputing
  the log values gave a lower original-scale error than imputing the raw values, for ranger
  (0.360 against 0.469) and for PLS (0.318 against 0.698 with `forceZero = TRUE` and 1.534
  without).
* **A 0 that means "not detected" must become `NA` before you log.** *Reasoning:* `log(0)` is
  `-Inf`, which `mvImpute()` refuses, and a 0 left as a value is treated as observed.
* **Usual order:** sample-level normalisation (for example PQN for dilution), then log
  transform, then imputation. *Reasoning:* the package does none of these, and the choice is
  yours.
* **`forceZero`** only matters for `method = "pls"` on non-negative, untransformed data, and PLS
  on untransformed intensities is not recommended. *Simulation:* it set 54 % of the imputed
  cells to zero (59 % were zero or negative without it). Keep `forceZero = FALSE` on log-scale
  data, where negative values are legitimate. Forests never need it, as they cannot predict
  outside the observed range.
* **`guess`** hardly changes the final values. *Simulation (PLS only):* with 3 components on log
  data, `"minVar"` against `"median"` gave 1.154 against 1.155, 0.953 against 0.954 and 0.944
  against 0.945 in the three patterns below. For the forests this is expected, because the
  iterations replace the starting values, but it was not evaluated (*reasoning*). Use
  `"minVar"` for omics data where low values go missing and `"median"` for covariates
  (*reasoning*).
* **Filter features with too many gaps before imputing; choose and report the threshold**
  (*reasoning*). Run time grows with the number of columns: the engine fits one model per column
  with gaps in every iteration.
* **Check the diagnostics:** `attr(out, "mvImpute")$not_converged` (cells still moving at the
  iteration cap) and `$oob_error` (out-of-bag error per column; it is estimated on observed rows,
  so when gaps are not missing at random it overstates how well the imputed values can be
  predicted).

**Gaps below a detection limit** (*Simulation*, log scale, 5 simulated data sets, SD-scaled
RMSE against the truth; lower is better; `bench/output/detection-limit.md`)

| Missing values are | Closest | Others |
|---|---|---|
| strictly the lowest 5 % of each column | column-minimum fill, 0.62 | PLS 1.13, ranger 1.54 |
| low values more likely but not only the lowest (a soft limit) | PLS, 0.92 (better than the minimum fill in 5 of 5 data sets) | minimum fill 1.12, ranger 1.17 |
| half lowest values, half at random | PLS, 0.92 | ranger 1.19, minimum fill 1.68 |

* Forests impute low values too high and never go below the observed minimum; PLS on log data
  can (26 % of its imputed cells when the lowest values were missing).
* PLS ran 4.5 to 5.6 times faster than ranger here (one worker).
* **Recommendation for metabolomics-type data where low values go missing:** log-transform, then
  `method = "pls"` (*simulation*; two or three components gave similar results). If you know the
  gaps are strictly below a detection limit, a minimum fill (each gap set to the lowest observed
  value of its column) was closer than PLS in this simulation, in all 5 simulated data sets
  (0.62 against 1.13; *simulation*, first row), so it is worth considering.
* Caveats: the simulated data are low-rank and log-normal, which favours PLS; "strictly the
  lowest" is the case most favourable to a minimum fill; only 2 and 3 components and only ranger
  among the forests were tried.

**By data type.** In the metabolomics row only the log scale and the choice of PLS rather than
ranger or a minimum fill are backed by simulation; the rest of that row and every other row is
*reasoning, not evaluated*. Vendor and platform details vary: check them against your own data.

| Data | Typical missingness | Impute on | `method` | `guess` | Notes |
|---|---|---|---|---|---|
| LC/GC-MS metabolomics | mostly low values below detection, some random (*reasoning*) | log (*simulation*) | `"pls"` (*simulation*) | `"minVar"` (*reasoning*) | see the table above; zeros meaning "not detected" to `NA` and normalisation before the log, and filtering features with many gaps first (*reasoning*) |
| MS proteomics, label-free DDA | a mix of low-abundance and random gaps | log2 intensities | `"pls"` or `"ranger"` | `"minVar"` | filter features first; report the threshold |
| MS proteomics, DIA | often reported to have fewer gaps than DDA; check whether those left concentrate at low intensity | log2 intensities | `"ranger"` or `"pls"` | `"median"` or `"minVar"` | |
| MS proteomics, TMT | often block-wise: a protein missing in a whole plex or batch | log2 intensities | only for scattered gaps | `"minVar"` or `"median"` | do not impute whole-batch gaps: it mixes the batch effect into the values |
| Olink NPX | often few; some values below the limit of detection (LOD) | already log2: no further transform | `"pls"` or `"ranger"` | `"minVar"` if the gaps are below the LOD, else `"median"` | NPX can be negative, so `forceZero = FALSE`; check whether your data already report values below the LOD, in which case little may be missing |
| SomaScan and other affinity panels | usually few | log (the signal is right-skewed) | `"ranger"` or `"pls"` | `"median"` | |
| Gene-expression microarray | often failed spots; some pipelines also set values missing by a detection p-value filter, which is tied to low intensity (check yours) | log2 intensities | `"ranger"` or `"pls"` | `"median"` (`"minVar"` if the gaps follow a detection filter) | |
| RNA-seq | zeros are real counts, not missing | a log-CPM-type scale | `"ranger"` or `"pls"` | `"median"` | do not impute zeros; impute only true `NA`; about 20,000 genes make runs slow (filter first) |
| DNA methylation | failed detection, linked to low total signal rather than low methylation | M-values, not beta values | `"ranger"` or `"pls"` | `"median"` | PLS on beta values can leave [0, 1] and `forceZero` gives no upper bound; 450K/EPIC arrays have hundreds of thousands of CpGs and one model is fitted per column per iteration, so the package is not designed for array-wide methylation imputation |
| Clinical covariates | any | original scale, factors as factors | `"ranger"` | `"median"` | |

## Assumptions

* **Pooled imputation.** Every row, including rows that a later cross-validation fold holds out,
  informs every imputed value. Downstream cross-validated or test-set performance is therefore
  optimistic, more so if the outcome is a column of `X`. Leave the outcome out of `X`, or pass it
  in `exclude`: excluded columns are returned unchanged (gaps stay gaps), are never predictors and
  are never imputed. `exclude` covers only the columns it names: a copy or a derived version of
  the outcome left in `X` still informs every model. `mvImpute()` cannot fit on training rows and
  apply the fit to test rows; to keep a test set untouched, impute it in a separate call (the
  models then learn from fewer rows, so its imputations are likely less accurate; not measured).
* **Missing at random.** The models assume that values are missing at random given the other
  columns. Values missing because they are below a detection limit (missing not at random) are
  imputed too high by the forest backends, which cannot predict below the observed range. In a
  simulation where the lowest 5 % of each column were missing, all forest imputations lay above
  the true values, and leaving the column minimum in place was much closer to the truth (SD-scaled
  RMSE 0.065 against 0.469 on the original scale; `bench/output/detection-limit.md`). That is one
  simulated setting, the one most favourable to a minimum fill. `method = "pls"` is linear and can
  extrapolate below the observed range, including to negative values on untransformed intensities
  (59 % of its imputations in that simulation). On log-scale data PLS was closer than ranger in
  the detection-limit simulations (`randomForest` was not run).
* **Single imputation.** The result is one deterministic completed table. Downstream standard
  errors that treat imputed values as observed are too small.
* **Linear dependencies.** Exact linear dependencies among numeric columns that you supply (your
  own one-hot columns, or parts and their total) can make the iterations oscillate. Recode a
  one-hot block as one factor column; factor columns are handled safely.

## Moving from StatTools::mvImpWrap

| `mvImpWrap()` (StatTools 0.0.917) | `mvImpute()` |
|---|---|
| `MAT` (matrix) | `X` (data frame or numeric matrix; the same class is returned) |
| `method = "RF", rfMeth = "ranger"` (what a default call ran) | `method = "ranger"` (the default) |
| `method = "RF", rfMeth = "rf"` | `method = "randomForest"` |
| `method = "PLS"` | `method = "pls"` (needs `mixOmics`) |
| `rfMeth = "Rborist"` | dropped |
| `guess = NULL` (column minimum) | `guess = "minVar"` (the default, as before); `"median"`, `"minTotal"` and a matrix of the shape of `X` are accepted as before |
| `forceZero`, `nComp` | `forceZero`, `nComp`, unchanged |
| `nCore` | `workers` (`NULL`: active plan, else up to 8 workers; `1`: serial) |
| `tol1`, `tol2` | `tol = c(tol1, tol2)`: now a fraction of the column's standard deviation, not relative to the current value, so old values do not carry over (defaults `c(0.05, 0.025)`) |
| `n1`, `n2` | `maxIter = c(n1, n2)` (defaults `c(15, 60)`) |
| `mvImp()` arguments `ntree`, `mtry`, `verbose` | `ntree = 100`, `mtry = 5`, `verbose = FALSE` |
| (none) | `seed`, `exclude` |
| messages by `cat()`, always on | messages by `message()`, silenceable; `verbose = TRUE` for per-iteration progress |

A model is now trained only on the rows where its column was originally observed (the old
engine also trained on converged imputed cells). Tables that were one-hot encoded for the old
engine should be passed as factor columns instead; see `NEWS.md`.

## Speed

Benchmarks on simulated data against StatTools 0.0.917 at the same number of workers, one run
each (`bench/output/benchmark.md` in the repository):

* Numeric tables (200 x 1000, 5 % gaps): new/old wall time 0.36 to 0.54, because the new engine
  stopped after 11 to 14 iterations where the old one ran to its cap of 75.
* The old engine ran two ranger threads in each worker and the new one runs one, so equal
  worker counts give the new engine half the threads.
* A small mixed table (900 x 20): new/old wall time 1.53 to 1.92, both engines at 75
  iterations. The higher cost per iteration may come from single-thread forests and
  dispatch overhead; this is a hypothesis, not tested.

## Related work

The iterative random-forest imputation in mvImputer is related to missForest:
Stekhoven DJ, Bühlmann P (2012). MissForest: non-parametric missing value imputation for
mixed-type data. *Bioinformatics* 28(1):112-118. <https://doi.org/10.1093/bioinformatics/btr597>.
mvImputer differs in its stopping rule (scaled by each column's SD), its PLS backend,
training only on originally observed rows, two stages, and seeding that is reproducible
across worker counts.

## Contact

Carl Brunius, carl.brunius@gmail.com (maintainer) or carl.brunius@nbis.se. Bug reports:
<https://github.com/CarlBrunius/mvImputer/issues>.
