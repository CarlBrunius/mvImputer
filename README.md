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

## Assumptions

* **Pooled imputation.** Every row, including rows that a later cross-validation fold holds out,
  informs every imputed value. Downstream cross-validated or test-set performance is therefore
  optimistic, more so if the outcome is a column of `X`. Leave the outcome out of `X`, or pass
  it in `exclude`: excluded columns are returned unchanged (gaps stay gaps), are never
  predictors and are never imputed. `exclude` covers only the columns it names: a copy or a
  derived version of the outcome left in `X` still informs every model. `mvImpute()` cannot
  fit on training rows and apply the fit to test rows; to keep a test set untouched, impute it
  in a separate call.
* **Missing at random.** The models assume that values are missing at random given the other
  columns. Values missing because they are below a detection limit (missing not at random) are
  imputed too high by the forest backends, which cannot predict below the observed range. In a
  simulation where the lowest 5 % of each column were missing, all forest imputations lay above
  the true values, and leaving the column minimum in place was much closer to the truth
  (SD-scaled RMSE about 0.06 against 0.44). That is one simulated setting, the one most
  favourable to a minimum fill. `method = "pls"` is linear and can extrapolate below the
  observed range, including to negative values on untransformed intensities (57 % of its
  imputations in one such simulation); no report evaluates PLS.
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

## Contact

Carl Brunius, carl.brunius@gmail.com (maintainer) or carl.brunius@nbis.se. Bug reports:
<https://github.com/CarlBrunius/mvImputer/issues>.
