# mvImputer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the R package mvImputer: one exported function, `mvImpute()`, a two-stage iterative
RF/PLS imputer with worker-count-invariant seeding, `future` parallelism, native factor handling
and an SD-scaled stopping rule.

**Architecture:** Port of StatTools `mvImp`/`mvImpWrap`, restructured into small internal
functions: validate → encode (working matrix `W` + design matrix) → per-stage iteration (all
features fitted from the previous iteration's values via `future.apply::future_lapply`, one
pre-generated L'Ecuyer stream per (stage, column)) → decode. Backends (ranger, randomForest, pls)
sit behind one internal `fit_predict()`.

**Tech Stack:** R (>= 4.1.0), future, future.apply, ranger, randomForest, withr, stats, parallel;
mixOmics (Suggests); testthat 3e; roxygen2.

**Spec:** `docs/superpowers/specs/2026-11-23-mvimputer-design.md` (approved 2026-11-23). Executors
read both. Section references below (§5.4 etc.) point into the spec.

## Global Constraints

- Package root: the repository root. Package name `mvImputer`, version `0.1.0`.
- Exports: `mvImpute()` only. Every other function is internal (`@noRd` or `@keywords internal`).
- Imports: `future`, `future.apply`, `ranger`, `randomForest`, `stats`, `parallel`, `withr`.
  Suggests: `mixOmics`, `testthat (>= 3.0.0)`. StatTools never appears in DESCRIPTION.
- License: `MIT + file LICENSE`; `Authors@R`: Carl Brunius, `carl.brunius@gmail.com`, roles
  `aut`, `cre`. Carl's NBIS address appears only in README.md.
- Namespace-qualified calls everywhere (`stats::sd`, `future::plan`); no `library()`/`require()`
  in `R/`; no global options changed; no `cat()`/`print()` for progress (use `message()`).
- RNG: the caller's `.Random.seed` and `RNGkind()` are restored on exit (§5.5); `set.seed()` is
  never called outside a guarded scope.
- Simulated data only, in tests, examples and bench scripts. Never read any folder outside the
  package root and a temporary directory.
- CRAN-ready: examples and tests use at most 2 workers; tests needing more use
  `skip_on_cran()`; `R CMD check --as-cran` with 0 errors, 0 warnings.
- `bench/`, `dev/`, `docs/` are listed in `.Rbuildignore`; logs go to `dev/logs/` (gitignored).
- **Commits:** Carl's rule is no commit or push until he says so. Every "Commit" step below runs
  only if Carl has allowed local commits for this execution; otherwise skip it. Never push.
  Commit author e-mail `carl.brunius@gmail.com`.
- Test command (keeps logs out of context; `<name>` is one plain word, one test file per run):
  - Tasks 1–6: `Rscript -e "testthat::test_local(filter = '<name>', reporter = 'summary', stop_on_failure = TRUE)" > dev/logs/test-<name>.log 2>&1; tail -n 40 dev/logs/test-<name>.log`
  - Tasks 7 onward (multisession workers must load an **installed** mvImputer; `load_all` code
    is invisible to them): `R CMD INSTALL --no-docs . > dev/logs/install.log 2>&1 && Rscript -e "testthat::test_local(filter = '<name>', load_package = 'installed', reporter = 'summary', stop_on_failure = TRUE)" > dev/logs/test-<name>.log 2>&1; tail -n 40 dev/logs/test-<name>.log`
  - Full suite: the Task 7 form without `filter`, log name `all`.
- R code that produces a result follows the `writing-auditable-r` contract: file header with
  purpose and usage, comments that say why, assertions marked BLOCKING or DIAGNOSTIC.

## Review Focus

1. **A data frame with an excluded column holding gaps** (for example an outcome with NA): the
   excluded column must come back with its gaps, and no imputation may change if its values
   change. Test in Task 8 (`exclude is invisible to the models`).
2. **A factor level that never occurs in the observed rows but is declared in `levels()`**: it
   must survive decoding and never be imputed. Test in Task 5 (`unused levels survive`) and
   Task 6 (`droplevels mapping`).
3. **A user who already set `future::plan(multisession, workers = 2)`**: `mvImpute()` must use
   it, not replace it, and leave it in place afterwards. Test in Task 3.
4. **Matrix input with column names missing**: output keeps the same (absent) dimnames, and
   messages refer to columns as `V1`, `V2`, ... Test in Task 5.
5. **A tiny table (n = 20) with one factor**: completes without error and without a warning
   from a backend (zero-variance predictors in small training sets). Test in Task 8.

## File Structure

| File | Responsibility |
|---|---|
| `DESCRIPTION`, `LICENSE`, `NAMESPACE` (roxygen), `.Rbuildignore`, `.gitignore` | Package metadata |
| `R/mvImputer-package.R` | Package doc, `@importFrom` nothing (qualified calls only) |
| `R/validate.R` | `check_input()`: errors and warnings of §5.1 |
| `R/encode.R` | `encode_input()`, `design_matrix()`, `first_guess()`, `decode_output()` |
| `R/converge.R` | `col_scales()`, `converge_cells()` (pure, §5.4) |
| `R/rng.R` | `resolve_seed()`, `local_rng_guard()`, `make_streams()` (§5.5) |
| `R/workers.R` | `local_workers()` (§5.1 `workers`) |
| `R/backends.R` | `fit_predict()` for ranger / randomForest / pls, numeric / factor (§5.3) |
| `R/stage.R` | `stage_jobs()`, `run_stage()` (§5.4 freezing, training rows) |
| `R/mvImpute.R` | Exported `mvImpute()`: orchestration, diagnostics, warnings |
| `tests/testthat/helper-sim.R` | `sim_numeric()`, `sim_mixed()`, `punch_holes()` |
| `tests/testthat/test-*.R` | One file per `R/` file, plus `test-mvImpute.R` |
| `bench/experiments.R`, `bench/calibrate_tol.R`, `bench/benchmark.R` | Studies, not shipped |
| `README.md`, `NEWS.md` | User docs |

## Internal data model (shared by Tasks 5–8)

- `enc` (from `encode_input()`), a list:
  - `W`: numeric matrix n × p over the **non-excluded** columns, in original order. Numeric and
    integer columns hold their values; factor-like columns hold integer level codes 1..K. `NA`
    where missing.
  - `miss`: logical n × p, `TRUE` where the original value was missing.
  - `type`: character p, each `"numeric"`, `"nominal"` or `"ordered"`.
  - `class`: character p, original class: `"numeric"`, `"integer"`, `"factor"`, `"ordered"`,
    `"character"`, `"logical"`.
  - `levels`: named list (factor-like columns only) of the full level set, as character.
  - `constant`: character, names of columns whose observed values are constant.
  - `exclude`: character, excluded column names.
  - `X`: the original input, kept for decoding.
  - `colnames`: names used internally (`V1..Vp` for unnamed matrices); `had_names`: logical.
- Design matrix `dm` (from `design_matrix(W, enc)`): list with `D` (numeric n × q: numeric
  columns as is, nominal one-hot columns named `<col>_<level>` for every declared level,
  ordered as their integer code) and `map` (named list: for each column of `W`, the integer
  indices of its columns in `D`).
- A **job** (from `stage_jobs()`): list with `col` (index into `W`), `train` (integer rows where
  that column is originally observed), `pred` (integer rows of its unfrozen missing cells).

---

### Task 1: Package scaffold and simulation helpers

Review tier: **light**.

**Files:**
- Create: `DESCRIPTION`, `LICENSE`, `.Rbuildignore`, `.gitignore`, `R/mvImputer-package.R`,
  `tests/testthat.R`, `tests/testthat/helper-sim.R`, `tests/testthat/test-helpers.R`,
  `dev/logs/.gitkeep`

**Interfaces:**
- Produces:
  - `sim_numeric(n, p, seed)` → numeric matrix n × p, colnames `F1..Fp`. Pinned process:
    `Z` n × 3 and loadings `L` 3 × p standard normal, `E` n × p normal with SD 0.5,
    `X = exp(0.5 * Z %*% L + E)` (positive, log-normal, omics-like).
  - `sim_mixed(n, seed, extra = 0)` → data.frame. Pinned process (`e(s)` = fresh normal noise
    with SD `s`; `z1`, `z2` standard normal latents): `num1 = z1 + e(0.3)`,
    `num2 = z1 + z2 + e(0.3)`, `num3 = z2 + e(0.3)`, `num4 = z1 - z2 + e(0.3)`,
    `int1 = as.integer(round(5 + 2 * z1 + e(0.5)))`,
    `fac3 = factor(c("a","b","c")[max.col(cbind(num1, num2, num3) + e(0.3))], levels = c("a","b","c"))`,
    `ord4` = `num4 + e(0.3)` cut at its quartiles into ordered `L1`<`L2`<`L3`<`L4`,
    `chr2 = ifelse(num2 + e(0.3) > 0, "x", "y")`, `lgl = num3 + e(0.3) > 0`, `const = 5`, then
    `extra` columns `noise1..` = `0.5 * z1 + e(1)` (benchmarks use `extra = 10` for 20 columns).
    Threshold tests in Tasks 6–8 were set against this process. If one fails, report it
    (DONE_WITH_CONCERNS); never loosen a threshold silently.
  - `punch_holes(X, frac, mechanism = c("mcar", "left"), seed, cols = NULL)` → X with `NA`s.
    Accepts a matrix or a data.frame. `"mcar"`: exactly `round(frac * <cells in cols>)` cells,
    drawn uniformly in `cols` (all if `NULL`), any column type. `"left"`: in each numeric column of `cols`, the lowest `frac` of values
    become NA (detection-limit pattern).
  All helpers scope their randomness with `withr::with_seed()`.

- [ ] **Step 1: Create the scaffold.** The git repository already exists (keep its `.gitignore`
  entries, including `NOTES.md`, and append the ones below). DESCRIPTION per Global
  Constraints (`Depends: R (>= 4.1.0)`, `Config/testthat/edition: 3`, `Encoding: UTF-8`,
  `BugReports: https://github.com/CarlBrunius/mvImputer/issues`). LICENSE in the CRAN MIT
  template: `YEAR: 2026` / `COPYRIGHT HOLDER: Carl Brunius`. `.Rbuildignore`: `^bench$`,
  `^dev$`, `^docs$`, `^\.git`, `^.*\.Rproj$`, `^NOTES.*\.md$`. `.gitignore` (append):
  `dev/logs/*.log`, `.Rhistory`.

- [ ] **Step 2: Write the failing tests** in `tests/testthat/test-helpers.R`:

```r
test_that("punch_holes mcar hits the requested fraction and leaves other cells alone", {
  X <- sim_numeric(100, 20, seed = 1)
  Y <- punch_holes(X, 0.1, "mcar", seed = 2)
  expect_equal(mean(is.na(Y)), 0.1, tolerance = 0.01)
  expect_identical(Y[!is.na(Y)], X[!is.na(Y)])
})
test_that("punch_holes left removes the lowest values of each column", {
  X <- sim_numeric(100, 5, seed = 1)
  Y <- punch_holes(X, 0.2, "left", seed = 2)
  for (j in 1:5) expect_lt(max(X[is.na(Y[, j]), j]), min(X[!is.na(Y[, j]), j]))
})
test_that("sim_mixed has the documented columns and classes", {
  d <- sim_mixed(300, seed = 1)
  expect_identical(vapply(d, function(x) class(x)[1], ""),
    c(num1 = "numeric", num2 = "numeric", num3 = "numeric", num4 = "numeric", int1 = "integer",
      fac3 = "factor", ord4 = "ordered", chr2 = "character", lgl = "logical", const = "numeric"))
  expect_false(anyNA(d))
})
test_that("helpers are reproducible and do not touch the caller's RNG", {
  set.seed(9); s <- .Random.seed
  expect_identical(sim_mixed(50, seed = 3), sim_mixed(50, seed = 3))
  expect_identical(s, .Random.seed)
})
```

- [ ] **Step 3: Run** the test command with `filter = 'helpers'`. Expected: FAIL, `could not find function "sim_numeric"`.
- [ ] **Step 4: Implement the three helpers** in `tests/testthat/helper-sim.R`. The factor
  columns must be predictable from their drivers (for example `fac3` = argmax of
  `num1..num3` plus noise) so accuracy tests later have signal.
- [ ] **Step 5: Run** again. Expected: PASS (4 tests). Also run
  `Rscript -e "devtools::document()"` and check it produces a NAMESPACE without errors.
- [ ] **Step 6: Commit** (if allowed): `git add DESCRIPTION LICENSE .Rbuildignore .gitignore R tests dev/logs/.gitkeep NAMESPACE` / `git commit -m "Scaffold package and simulation helpers"`

---

### Task 2: Stopping rule

Review tier: **full**.

**Files:**
- Create: `R/converge.R`, `tests/testthat/test-converge.R`

**Interfaces:**
- Produces:
  - `col_scales(W, miss, type)` → **unnamed** numeric vector length `ncol(W)`: SD of observed values for
    `"numeric"` columns (0 for a constant column), `NA_real_` for factor-like columns.
  - `converge_cells(old, new, prev_same, scale, is_factor, tol)` → list with `converged`
    (logical, same length as `new`) and `same` (logical, `new == old`). All arguments are
    vectors of equal length except `tol` (scalar). Numeric cell: converged when
    `old == new` or `abs(new - old) / scale < tol`. Factor cell: converged when `same & prev_same`
    (two consecutive unchanged predictions). Stops with an error containing
    `"non-finite prediction"` if any `new` is not finite. "Still left" is `!converged`.

- [ ] **Step 1: Write the failing tests:**

```r
test_that("a negative cell that doubles is not converged (old rule said it was)", {
  expect_false(converge_cells(-1, -2, FALSE, 1, FALSE, 0.05)$converged)
})
test_that("a cell at zero can converge", {
  expect_true(converge_cells(0, 0.001, FALSE, 1, FALSE, 0.05)$converged)
})
test_that("a change exactly at tol is not converged", {
  expect_false(converge_cells(0, 0.5, FALSE, 10, FALSE, 0.05)$converged)
})
test_that("the rule is invariant to shifting the data", {
  old <- c(-3, 0, 2, 7); new <- old + c(0.01, -0.2, 0.03, 0.5); s <- rep(2, 4)
  expect_identical(converge_cells(old, new, rep(FALSE, 4), s, rep(FALSE, 4), 0.05)$converged,
                   converge_cells(old + 1000, new + 1000, rep(FALSE, 4), s, rep(FALSE, 4), 0.05)$converged)
})
test_that("a factor cell needs two consecutive unchanged predictions", {
  expect_false(converge_cells(2, 2, FALSE, NA, TRUE, 0.05)$converged)
  expect_true(converge_cells(2, 2, TRUE, NA, TRUE, 0.05)$converged)
  expect_false(converge_cells(2, 3, TRUE, NA, TRUE, 0.05)$converged)
})
test_that("non-finite predictions are an error, never silently converged", {
  expect_error(converge_cells(1, NaN, FALSE, 1, FALSE, 0.05), "non-finite prediction")
})
test_that("col_scales uses observed values only", {
  W <- cbind(a = c(1, 2, 3, 100), b = c(5, 5, 5, NA), f = c(1, 2, 1, 2))
  miss <- cbind(a = c(FALSE, FALSE, FALSE, TRUE), b = c(FALSE, FALSE, FALSE, TRUE), f = FALSE)
  expect_equal(col_scales(W, miss, c("numeric", "numeric", "nominal")), c(1, 0, NA))
})
```

- [ ] **Step 2: Run** with `filter = 'converge'`. Expected: FAIL, function not found.
- [ ] **Step 3: Implement** both functions in `R/converge.R`, with a header comment recording why
  the old relative rule was replaced (spec F2) and why factors need two unchanged predictions.
- [ ] **Step 4: Run.** Expected: PASS (7 tests).
- [ ] **Step 5: Commit** (if allowed): `git commit -m "Add SD-scaled stopping rule"` (add the two files by name).

---

### Task 3: Seeds, RNG isolation and workers

Review tier: **full**.

**Files:**
- Create: `R/rng.R`, `R/workers.R`, `tests/testthat/test-rng.R`, `tests/testthat/test-workers.R`

**Interfaces:**
- Produces:
  - `resolve_seed(seed)` → integer scalar. `NULL`: one draw
    `sample.int(.Machine$integer.max, 1L)` from the caller's RNG. Otherwise validates a single
    finite whole number and returns it as integer.
  - `local_rng_guard(envir = parent.frame())`: restores `.Random.seed` (including its absence)
    with `withr::local_preserve_seed(envir)` and the three `RNGkind()` values with a
    `withr::defer(RNGkind(...), envir)` registered after it. No hand-written `assign()` into
    `globalenv()`.
  - `make_streams(seed, n_stage, n_col)` → list of `n_stage` lists of `n_col` integer vectors of
    length 7 (L'Ecuyer-CMRG `.Random.seed` values), all distinct, a pure function of the
    arguments: sets all three kinds explicitly, `RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")`,
    then `set.seed(seed)`, inside `local_rng_guard()`; chains `parallel::nextRNGStream()`.
    `n_col` counts the columns of `W` (non-excluded), so streams are keyed by position in `W`.
    *(Plan-review decision: adding or removing excluded columns then cannot change any
    imputation. Spec §5.5 said "original column index"; the spec is amended to "column of W".)*
  - `stream_for_iteration(stream, iter, crn)` → the stream for an iteration: `stream` itself when
    `crn = TRUE`; otherwise `parallel::nextRNGSubStream()` applied `iter - 1` times.
  - `local_workers(workers, envir = parent.frame())` → integer, the number of workers in effect.
    Rules (§5.1): `NULL` and current plan not sequential → leave the plan; `NULL` and sequential →
    `k <- min(8, future::availableCores() - 1)`, multisession with `k` workers if `k >= 2`, else
    stay sequential (this also stops nested worker pools when `mvImpute()` is itself called inside
    a future, where `availableCores()` is 1); `1` → `future::plan(future::sequential)`; `k > 1` →
    multisession with `k`. Any plan it sets is restored when `envir` exits.

- [ ] **Step 1: Write the failing tests** (`test-rng.R`):

```r
test_that("make_streams is deterministic, distinct and L'Ecuyer", {
  s <- make_streams(42L, 2L, 3L)
  expect_identical(s, make_streams(42L, 2L, 3L))
  expect_length(s, 2); expect_length(s[[1]], 3)
  flat <- unlist(lapply(s, function(x) vapply(x, paste, "", collapse = ",")))
  expect_false(anyDuplicated(flat) > 0)
  expect_identical(s[[1]][[1]][1], 10407L)
})
test_that("make_streams leaves the caller's RNG exactly as it was", {
  set.seed(1); s0 <- .Random.seed
  make_streams(42L, 2L, 3L)
  expect_identical(.Random.seed, s0); expect_identical(RNGkind()[1], "Mersenne-Twister")
})
test_that("an absent .Random.seed stays absent", {
  withr::local_preserve_seed()
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  make_streams(1L, 1L, 1L)
  expect_false(exists(".Random.seed", globalenv()))
})
test_that("resolve_seed(NULL) follows set.seed and costs exactly one draw", {
  set.seed(7); a <- resolve_seed(NULL); s1 <- .Random.seed
  set.seed(7); b <- resolve_seed(NULL)
  expect_identical(a, b)
  set.seed(7); sample.int(.Machine$integer.max, 1L)
  expect_identical(s1, .Random.seed)
})
test_that("crn = FALSE gives a different stream each iteration, crn = TRUE the same", {
  s <- make_streams(1L, 1L, 1L)[[1]][[1]]
  expect_identical(stream_for_iteration(s, 3L, TRUE), s)
  expect_false(identical(stream_for_iteration(s, 3L, FALSE), stream_for_iteration(s, 2L, FALSE)))
})
```

  And `test-workers.R`:

```r
test_that("workers = 2 sets multisession for the call and restores the old plan", {
  future::plan(future::sequential); withr::defer(future::plan(future::sequential))
  f <- function() { n <- local_workers(2L); c(n, future::nbrOfWorkers()) }
  expect_equal(f(), c(2, 2))
  expect_true(inherits(future::plan(), "sequential"))
})
test_that("an active non-sequential plan is used and left in place", {
  future::plan(future::multisession, workers = 2); withr::defer(future::plan(future::sequential))
  f <- function() local_workers(NULL)
  expect_equal(f(), 2)
  expect_true(inherits(future::plan(), "multisession"))
})
test_that("workers = 1 runs sequentially", {
  f <- function() { local_workers(1L); future::nbrOfWorkers() }
  expect_equal(f(), 1)
})
```

- [ ] **Step 2: Run** with `filter = 'rng'`, then `filter = 'workers'`. Expected: FAIL, functions not found.
- [ ] **Step 3: Implement** `R/rng.R` and `R/workers.R`. Header of `rng.R` records why the guard
  is required (spec §5.5: `future_lapply` advances the caller's `.Random.seed` even with
  pre-generated seeds; verified 2026-11-23 with future.apply 1.20.2).
- [ ] **Step 4: Run.** Expected: PASS (8 tests).
- [ ] **Step 5: Commit** (if allowed): `git commit -m "Add scoped seeding and worker plan handling"`.

**Checkpoint:** seeds and the worker contract are what Tasks 7–8 build on; review before the data model.

---

### Task 4: Input checks

Review tier: **full** (guards against silent failure).

**Files:**
- Create: `R/validate.R`, `tests/testthat/test-validate.R`

**Interfaces:**
- Produces: `check_input(X, exclude, tol, maxIter, method, nComp)` → invisible `NULL`;
  `stop()` with the messages below, `warning()` for the two warnings. Runs before encoding.
  Checks apply to non-excluded columns only, except "unknown column".

| Condition | Kind | Message contains |
|---|---|---|
| `X` neither data.frame nor matrix | error | `"data.frame or matrix"` |
| matrix not numeric | error | `"numeric matrix"` |
| `exclude` names not in `colnames(X)` | error | `"unknown column"` |
| fewer than 2 non-excluded columns | error | `"at least 2 columns"` |
| column type not numeric/integer/factor/ordered/character/logical | error | `"unsupported column type"` |
| any `Inf`/`-Inf` | error | `"Inf"` and `"set them to NA"` |
| a column with gaps and < 2 observed values | error | `"fewer than 2 observed"` |
| character/factor column with ≥ 3 observed values, all distinct | error | `"identifier"` and `"exclude"` |
| `length(tol) != length(maxIter)` | error | `"one value per stage"` |
| any `tol <= 0` or `maxIter < 1` | error | `"tol"` / `"maxIter"` |
| `method == "pls"` and mixOmics not installed | error | `"mixOmics"` |
| a row with no observed value | warning | `"no observed value"` |
| numeric 0/1 columns forming a block that sums to 1 in every complete row | warning | `"one-hot"` and `"factor"` |

  The mixOmics check goes through an internal `has_pkg(pkg)` wrapper around
  `requireNamespace(pkg, quietly = TRUE)`, so the test can use
  `testthat::local_mocked_bindings(has_pkg = function(pkg) FALSE)`.
  Approach for the one-hot warning: among numeric columns whose observed values are all in
  {0, 1}, two columns are linked when never both 1 in rows where both are observed; within each
  connected group of ≥ 2 columns, warn if the row sum is 1 in every row complete for the group.

- [ ] **Step 1: Write the failing tests:** one `test_that()` per row of the table, each building a
  minimal input from `sim_mixed(40, seed = 1)` or a small literal matrix and asserting
  `expect_error(..., "<message>")` or `expect_warning(..., "<message>")`; plus
  `test_that("a valid mixed table passes silently", expect_silent(check_input(punch_holes(sim_mixed(100, 1), 0.1, "mcar", 2, cols = 1:9), NULL, c(0.05, 0.025), c(15L, 60L), "ranger", 2L)))`
  and `test_that("an excluded identifier column passes", ...)` with `id = as.character(1:40)`
  and `exclude = "id"`.
- [ ] **Step 2: Run** with `filter = 'validate'`. Expected: FAIL.
- [ ] **Step 3: Implement** `check_input()` in `R/validate.R`.
- [ ] **Step 4: Run.** Expected: PASS (15 tests).
- [ ] **Step 5: Commit** (if allowed): `git commit -m "Add input checks"`.

---

### Task 5: Encoding, design matrix, first guess, decoding

Review tier: **full** (factor round trips, observed-value invariance).

**Files:**
- Create: `R/encode.R`, `tests/testthat/test-encode.R`

**Interfaces:**
- Consumes: nothing from earlier tasks except `check_input()` having run.
- Produces (data model above):
  - `encode_input(X, exclude)` → `enc`.
  - `design_matrix(W, enc)` → `dm` (`D`, `map`). Missing cells of factor-like columns must
    already be filled (codes) before this is called; numeric NA stays NA only before the first
    guess.
  - `first_guess(W, enc, guess)` → `W` with every missing cell filled. Constant columns: always
    their constant, whatever `guess` is. Numeric: `"minVar"` column minimum of observed values,
    `"median"` column median, `"minTotal"` minimum over all numeric observed values;
    factor-like: always the most frequent observed code (ties: lowest code), per spec §5.1. A
    matrix/data.frame `guess` of `X`'s dimensions is read at missing cells of non-constant
    numeric columns only and must be non-NA there, else error `"guess"`.
  - `decode_output(W, enc)` → object of `class(enc$X)` with dimnames of `X`; types restored per
    §5.1 (integer with gaps → double); excluded columns copied from `enc$X` untouched; constant
    columns filled with their constant.

- [ ] **Step 1: Write the failing tests:**

```r
test_that("encode then decode returns the input unchanged, gaps included", {
  d <- sim_mixed(200, seed = 1)
  d$int1 <- as.numeric(d$int1)              # integer-with-gaps returns double by design
  dh <- punch_holes(d, 0.1, "mcar", seed = 2, cols = 1:9)
  enc <- encode_input(dh, exclude = NULL)
  expect_identical(decode_output(enc$W, enc), dh)
})
test_that("an integer column without gaps stays integer", {
  d <- sim_mixed(50, seed = 1); d$num1[3] <- NA
  expect_type(decode_output(encode_input(d, NULL)$W, encode_input(d, NULL))$int1, "integer")
})
test_that("unused levels survive and are never imputed", {
  d <- sim_mixed(100, seed = 1); d$fac3 <- factor(d$fac3, levels = c("a", "b", "c", "zz"))
  d$fac3[1:5] <- NA
  enc <- encode_input(d, NULL); W <- first_guess(enc$W, enc, "minVar")
  out <- decode_output(W, enc)
  expect_identical(levels(out$fac3), c("a", "b", "c", "zz")); expect_false(any(out$fac3 == "zz"))
})
test_that("design matrix: one-hot for nominal, integer code for ordered, own columns mapped", {
  d <- sim_mixed(30, seed = 1); enc <- encode_input(d, NULL); dm <- design_matrix(enc$W, enc)
  expect_true(all(c("fac3_a", "fac3_b", "fac3_c", "chr2_x", "chr2_y", "lgl_FALSE", "lgl_TRUE") %in% colnames(dm$D)))
  expect_identical(colnames(dm$D)[dm$map[["fac3"]]], c("fac3_a", "fac3_b", "fac3_c"))
  expect_identical(colnames(dm$D)[dm$map[["ord4"]]], "ord4")
  expect_setequal(unique(dm$D[, "ord4"]), 1:4)
})
test_that("excluded columns are absent from W and D but returned untouched", {
  d <- sim_mixed(30, seed = 1); d$y <- c(NA, rnorm(29)); enc <- encode_input(d, exclude = "y")
  expect_false("y" %in% colnames(enc$W))
  expect_identical(decode_output(first_guess(enc$W, enc, "minVar"), enc)$y, d$y)
})
test_that("first guesses: minVar, median, minTotal, factor mode, user matrix at missing cells only", {
  m <- cbind(a = c(1, NA, 3, 10), b = c(-5, 2, NA, 4)); enc <- encode_input(m, NULL)
  expect_equal(first_guess(enc$W, enc, "minVar")[2, "a"], 1)
  expect_equal(first_guess(enc$W, enc, "median")[3, "b"], 2)
  expect_equal(first_guess(enc$W, enc, "minTotal")[2, "a"], -5)
  g <- m; g[] <- 99; expect_equal(first_guess(enc$W, enc, g)[1, "a"], 1)
  expect_error(first_guess(enc$W, enc, m), "guess")
})
test_that("unnamed matrix: dimnames stay absent, internal names are V1..Vp", {
  m <- matrix(c(1, NA, 3, 4, 5, 6), 3); enc <- encode_input(m, NULL)
  expect_identical(enc$colnames, c("V1", "V2"))
  expect_null(dimnames(decode_output(first_guess(enc$W, enc, "minVar"), enc)))
})
test_that("constant columns are detected", {
  expect_identical(encode_input(sim_mixed(30, seed = 1), NULL)$constant, "const")
})
```

- [ ] **Step 2: Run** with `filter = 'encode'`. Expected: FAIL.
- [ ] **Step 3: Implement** `R/encode.R` (four functions; ≤150 lines per function).
- [ ] **Step 4: Run.** Expected: PASS (8 tests).
- [ ] **Step 5: Commit** (if allowed): `git commit -m "Add encoding and decoding of mixed tables"`.

**Checkpoint:** the data model (`enc`, `dm`) is what backends and the stage engine consume.

---

### Task 6: Backends

Review tier: **full**.

**Files:**
- Create: `R/backends.R`, `tests/testthat/test-backends.R`

**Interfaces:**
- Consumes: `dm$D` columns as predictors (Task 5).
- Produces: `fit_predict(method, y, x_train, x_pred, is_factor, n_levels, scale, ntree, mtry, nComp)`
  → list with `pred` (numeric for numeric targets; integer codes in `1..n_levels` for factors)
  and `oob` (numeric scalar; `NA_real_` for pls). `y` is numeric (numeric target) or integer
  codes (factor target). Draws all randomness from the current R RNG (the caller sets it).
  Per fit, in this order: drop `x_train` columns with zero variance (and the same columns of
  `x_pred`); factor targets fitted on `droplevels(factor(y))` and predicted codes mapped back to
  `1..n_levels`; `mtry <- min(mtry, ncol(x_train))`; `nComp <- min(nComp, ncol(x_train))`;
  `stop("non-finite prediction")` on any non-finite numeric prediction. `x_train`/`x_pred` are
  always matrices (callers subset with `drop = FALSE`); one prediction row and one predictor
  column must work. If no predictor survives the zero-variance drop: predict the training mean
  (numeric) or most frequent code (factor), `oob = NA_real_`. randomForest's regression warning
  for a response with five or fewer unique values is muffled for that one message only
  (`withCallingHandlers`), because 0/1 and small-integer numeric columns are legitimate targets.
  - ranger: `ranger::ranger(x =, y =, num.trees = ntree, mtry =, num.threads = 1, seed = NULL)`;
    factor targets with `probability = TRUE`, class = argmax. `oob`: numeric
    `sqrt(fit$prediction.error) / scale`; factor = OOB misclassification computed from
    `fit$predictions` (OOB class probabilities) versus `y`.
  - randomForest: `randomForest::randomForest(x =, y =, ntree =, mtry =)`; `oob`: numeric
    `sqrt(fit$mse[ntree]) / scale`; factor `fit$err.rate[ntree, "OOB"]`.
  - pls: `mixOmics::pls(X =, Y =, ncomp = nComp, mode = "regression", scale = TRUE)`; factor `Y` is
    the indicator matrix of the observed levels; class = largest predicted score. `oob = NA_real_`.

- [ ] **Step 1: Write the failing tests:**

```r
for (m in c("ranger", "randomForest", "pls")) {
  test_that(paste(m, "predicts a numeric target"), {
    if (m == "pls") skip_if_not_installed("mixOmics")
    x <- withr::with_seed(1, matrix(rnorm(200 * 4), 200)); colnames(x) <- paste0("x", 1:4)
    y <- x[, 1] + withr::with_seed(2, rnorm(200, sd = 0.2))
    r <- withr::with_seed(3, fit_predict(m, y[1:150], x[1:150, ], x[151:200, ], FALSE, NA, sd(y), 100, 2, 2))
    expect_gt(cor(r$pred, y[151:200]), 0.9)
  })
  test_that(paste(m, "predicts a factor target as a valid code"), {
    if (m == "pls") skip_if_not_installed("mixOmics")
    x <- withr::with_seed(1, matrix(rnorm(300 * 3), 300)); colnames(x) <- paste0("x", 1:3)
    y <- max.col(x)                          # class = which column is largest
    r <- withr::with_seed(3, fit_predict(m, y[1:200], x[1:200, ], x[201:300, ], TRUE, 4L, NA, 100, 2, 2))
    expect_true(all(r$pred %in% 1:4)); expect_gt(mean(r$pred == y[201:300]), 0.8)
  })
}
test_that("a constant predictor is dropped without backend warnings", {
  skip_if_not_installed("mixOmics")
  x <- cbind(withr::with_seed(1, matrix(rnorm(60 * 3), 60)), k = 1); colnames(x)[1:3] <- paste0("x", 1:3)
  for (m in c("ranger", "randomForest", "pls"))
    expect_no_warning(withr::with_seed(2, fit_predict(m, x[1:50, 1], x[1:50, -1], x[51:60, -1], FALSE, NA, 1, 50, 2, 2)))
})
test_that("droplevels mapping: a level absent from training is never predicted", {
  x <- withr::with_seed(1, matrix(rnorm(100 * 2), 100)); colnames(x) <- c("x1", "x2")
  y <- ifelse(x[, 1] > 0, 1L, 3L)                    # level 2 of 3 never observed
  r <- withr::with_seed(2, fit_predict("randomForest", y[1:80], x[1:80, ], x[81:100, ], TRUE, 3L, NA, 50, 1, 2))
  expect_true(all(r$pred %in% c(1L, 3L)))
})
test_that("oob: near 0 for an exact function, near 1 for noise, NA for pls", {
  x <- withr::with_seed(1, matrix(rnorm(300 * 3), 300)); colnames(x) <- paste0("x", 1:3)
  exact <- withr::with_seed(2, fit_predict("ranger", 2 * x[, 1], x, x[1:2, ], FALSE, NA, sd(2 * x[, 1]), 200, 3, 2))
  noise <- withr::with_seed(2, fit_predict("ranger", rnorm(300), x, x[1:2, ], FALSE, NA, 1, 200, 3, 2))
  expect_lt(exact$oob, 0.3); expect_gt(noise$oob, 0.8)
  skip_if_not_installed("mixOmics")
  expect_true(is.na(fit_predict("pls", x[, 1], x[, -1], x[1:2, -1], FALSE, NA, 1, 10, 2, 2)$oob))
})
test_that("one prediction row, one predictor, and no surviving predictor all work", {
  x <- withr::with_seed(1, matrix(rnorm(50 * 2), 50)); colnames(x) <- c("x1", "x2"); y <- x[, 1]
  for (m in c("ranger", "randomForest"))
    expect_length(withr::with_seed(2, fit_predict(m, y[1:49], x[1:49, , drop = FALSE], x[50, , drop = FALSE], FALSE, NA, 1, 50, 2, 2))$pred, 1)
  expect_length(withr::with_seed(2, fit_predict("ranger", y, x[, 1, drop = FALSE], x[1:3, 1, drop = FALSE], FALSE, NA, 1, 50, 2, 2))$pred, 3)
  k <- cbind(k = rep(1, 50))
  r <- fit_predict("ranger", y, k, k[1:2, , drop = FALSE], FALSE, NA, 1, 50, 2, 2)
  expect_equal(r$pred, rep(mean(y), 2)); expect_true(is.na(r$oob))
})
test_that("randomForest regression on a 0/1 numeric target gives no warning", {
  x <- withr::with_seed(1, matrix(rnorm(80 * 2), 80)); colnames(x) <- c("x1", "x2")
  expect_no_warning(withr::with_seed(2, fit_predict("randomForest", as.numeric(x[, 1] > 0), x, x[1:3, ], FALSE, NA, 0.5, 50, 1, 2)))
})
test_that("same RNG state gives identical ranger and randomForest output", {
  x <- withr::with_seed(1, matrix(rnorm(100 * 3), 100)); colnames(x) <- paste0("x", 1:3); y <- x[, 1]
  for (m in c("ranger", "randomForest"))
    expect_identical(withr::with_seed(5, fit_predict(m, y, x, x[1:5, ], FALSE, NA, 1, 50, 2, 2)),
                     withr::with_seed(5, fit_predict(m, y, x, x[1:5, ], FALSE, NA, 1, 50, 2, 2)))
})
```

- [ ] **Step 2: Run** with `filter = 'backends'`. Expected: FAIL.
- [ ] **Step 3: Verify the spec's open hypothesis first** (§5.1 `oob_error`): fit a ranger
  probability forest on a toy factor and compare `fit$prediction.error` with the OOB
  misclassification from `fit$predictions`. Record the finding (one line, with the ranger
  version) in the header of `R/backends.R`. The implementation uses the misclassification rate
  either way.
- [ ] **Step 4: Implement** `fit_predict()` in `R/backends.R` (one internal helper per backend if
  the function would exceed 150 lines).
- [ ] **Step 5: Run.** Expected: all PASS (pls tests skip without mixOmics).
- [ ] **Step 6: Commit** (if allowed): `git commit -m "Add ranger, randomForest and PLS backends"`.

---

### Task 7: Stage engine

Review tier: **full** (training rows, freezing, seeding, the no-oscillation property).

**Files:**
- Create: `R/stage.R`, `tests/testthat/test-stage.R`

**Interfaces:**
- Consumes: `encode_input`, `design_matrix`, `first_guess` (Task 5); `col_scales`,
  `converge_cells` (Task 2); `make_streams`, `stream_for_iteration` (Task 3); `fit_predict` (Task 6).
- Produces:
  - `stage_jobs(miss, frozen, targets)` → list of jobs (data model above), one per target
    column with at least one unfrozen missing cell. `targets`: integer indices of columns that
    are imputed (not constant). `train` is always `which(!miss[, col])`, whatever is frozen.
  - `run_stage(W, enc, scales, streams, tol, maxIter, method, ntree, mtry, nComp, forceZero, crn = TRUE, verbose = FALSE)`
    → list with `W` (updated), `iterations` (integer), `not_converged` (data.frame with columns
    `row` (integer) and `column` (character, original column name)), `oob` (named numeric, last
    fit per target column). Targets are the columns with missing cells that are not in
    `enc$constant`; constant columns' cells start frozen (already filled by `first_guess()`), so
    they never keep a stage iterating. `verbose = TRUE`: one `message()` per iteration
    (`"Iteration <i>: <n> cells still iterating"`). All subsetting of `dm$D` uses
    `drop = FALSE`. One iteration: build
    `dm <- design_matrix(W, enc)` once; `jobs <- stage_jobs(...)`; run
    `future.apply::future_lapply(jobs, <fit job>, future.seed = <list: stream_for_iteration(streams[[job$col]], iter, crn) per job>)`;
    each job calls `fit_predict()` with `y = W[job$train, job$col]`,
    `x_train = dm$D[job$train, -dm$map[[job$col]]]`, `x_pred = dm$D[job$pred, -dm$map[[job$col]]]`;
    after all jobs return, write predictions into unfrozen cells only (`forceZero` floors
    numeric predictions at 0), apply `converge_cells()` per cell, freeze converged cells. Stop
    when no unfrozen cells remain or at `maxIter`.

- [ ] **Step 1: Write the failing tests:**

```r
test_that("training rows are the originally observed rows, even after freezing", {
  miss <- cbind(a = c(TRUE, TRUE, FALSE, FALSE), b = c(FALSE, TRUE, FALSE, FALSE))
  frozen <- cbind(a = c(TRUE, FALSE, FALSE, FALSE), b = FALSE)
  jobs <- stage_jobs(miss, frozen, targets = 1:2)
  expect_identical(jobs[[1]]$train, 3:4); expect_identical(jobs[[1]]$pred, 2L)
  expect_identical(jobs[[2]]$train, c(1L, 3L, 4L))
})
test_that("observed cells are never changed by a stage", {
  X <- punch_holes(sim_numeric(80, 10, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(1L, 1L, ncol(W0))[[1]]
  out <- run_stage(W0, enc, col_scales(enc$W, enc$miss, enc$type), st, 0.05, 5L, "ranger", 50, 3, 2, FALSE)
  expect_identical(out$W[!enc$miss], W0[!enc$miss])
})
test_that("a stage gives identical results with 1 and 2 workers", {
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(11L, 1L, ncol(W0))[[1]]; s <- col_scales(enc$W, enc$miss, enc$type)
  run <- function(k) { future::plan(future::multisession, workers = k); on.exit(future::plan(future::sequential))
    run_stage(W0, enc, s, st, 0.05, 4L, "ranger", 50, 3, 2, FALSE) }
  expect_identical(run(1), run(2))
})
test_that("factors with informative predictors converge without oscillating", {
  d <- sim_mixed(300, seed = 4); truth <- d$fac3
  d$fac3[withr::with_seed(5, sample(300, 45))] <- NA
  enc <- encode_input(d, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(2L, 1L, ncol(W0))[[1]]
  out <- run_stage(W0, enc, col_scales(enc$W, enc$miss, enc$type), st, 0.025, 20L, "ranger", 100, 5, 2, FALSE)
  # Accuracy is the real check: with sibling indicators as predictors it fell to ~0.58 (spec F3).
  j <- which(colnames(enc$W) == "fac3"); m <- enc$miss[, j]
  expect_gt(mean(enc$levels$fac3[out$W[m, j]] == truth[m]), 0.75)
})
test_that("oob is reported per imputed column", {
  X <- punch_holes(sim_numeric(80, 6, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  out <- run_stage(W0, enc, col_scales(enc$W, enc$miss, enc$type), make_streams(1L, 1L, 6L)[[1]], 0.05, 3L, "ranger", 50, 3, 2, FALSE)
  expect_setequal(names(out$oob), colnames(X)[colSums(is.na(X)) > 0])
})
```

- [ ] **Step 2: Run** with `filter = 'stage'`. Expected: FAIL.
- [ ] **Step 3: Implement** `stage_jobs()` and `run_stage()` in `R/stage.R`. Header records the
  two deliberate departures from mvImp: training on originally observed rows only (Carl,
  2026-11-23) and factors fitted as one target (spec F3).
- [ ] **Step 4: Run.** Expected: PASS (5 tests).
- [ ] **Step 5: Commit** (if allowed): `git commit -m "Add stage engine with frozen cells and per-column streams"`.

**Checkpoint:** the engine is complete; review the seeding path end to end before the public function wraps it.

---

### Task 8: `mvImpute()`

Review tier: **full**.

**Files:**
- Create: `R/mvImpute.R`, `tests/testthat/test-mvImpute.R`

**Interfaces:**
- Consumes: everything above.
- Produces: the exported `mvImpute()` with the exact signature of spec §5.1:
  `mvImpute(X, method = c("ranger", "randomForest", "pls"), exclude = NULL, seed = NULL, workers = NULL, tol = c(0.05, 0.025), maxIter = c(15, 60), guess = "minVar", forceZero = FALSE, ntree = 100, mtry = 5, nComp = 2, verbose = FALSE)`.
  Order: `method <- match.arg(method)`; `check_input()`; `seed <- resolve_seed(seed)` (before the
  guard, so `seed = NULL` costs exactly one draw); `local_rng_guard()`; `encode_input()`; return
  early with a `message()` if nothing is missing; `local_workers(workers)` (after the early
  return, so no worker pool starts for nothing); `first_guess()`;
  `make_streams(seed, length(tol), ncol(W))`; one `run_stage()` per stage, all missing cells
  unfrozen at the start of each stage; `decode_output()`; attach diagnostics.
  Diagnostics `attr(out, "mvImpute")`: `seed`, `rng_kind = "L'Ecuyer-CMRG"`, `method`,
  `settings` (list of the remaining arguments), `workers` (integer), `iterations` (integer per
  stage), `not_converged` (list per stage of data.frames `row`, `column`), `time` (numeric
  seconds per stage), `oob_error` (final stage), `rows_no_observed` (integer).
  `warning()` only when the final stage ends with unconverged cells, naming the count.
  `message()` for each constant column that had gaps, exactly
  `"Column '<name>' is constant; its missing values were filled with that value."`
  `verbose = TRUE`: one `message()` per stage start and per iteration (iteration number, cells
  still iterating); `verbose = FALSE` emits no progress messages.
  Roxygen help includes a section "Assumptions" with the four statements of spec §3
  ("Assumptions the help page must state"), the `exclude` advice for outcomes, and the
  `guess = "median"` advice for covariate tables.

- [ ] **Step 1: Write the failing tests** (one helper strips the attribute:
  `strip <- function(x) { attr(x, "mvImpute") <- NULL; x }`):

```r
d0 <- sim_mixed(150, seed = 1)
dh <- punch_holes(d0, 0.08, "mcar", seed = 2, cols = c(1:4, 6:9))

test_that("same seed, same data: identical with 1 and 2 workers", {
  a <- mvImpute(dh, seed = 1, workers = 1, maxIter = c(3, 5))
  b <- mvImpute(dh, seed = 1, workers = 2, maxIter = c(3, 5))
  expect_identical(strip(a), strip(b))
  expect_identical(attr(a, "mvImpute")$iterations, attr(b, "mvImpute")$iterations)
})
test_that("identical with 4 workers too", {
  skip_on_cran()
  expect_identical(strip(mvImpute(dh, seed = 1, workers = 1, maxIter = c(3, 5))),
                   strip(mvImpute(dh, seed = 1, workers = 4, maxIter = c(3, 5))))
})
test_that("different seeds differ; set.seed reproduces seed = NULL", {
  expect_false(identical(strip(mvImpute(dh, seed = 1, workers = 1, maxIter = c(2, 2))),
                         strip(mvImpute(dh, seed = 2, workers = 1, maxIter = c(2, 2)))))
  set.seed(5); a <- mvImpute(dh, workers = 1, maxIter = c(2, 2))
  set.seed(5); b <- mvImpute(dh, workers = 1, maxIter = c(2, 2))
  expect_identical(strip(a), strip(b))
  expect_identical(attr(a, "mvImpute")$seed, attr(b, "mvImpute")$seed)
})
test_that("caller RNG untouched with a seed, one draw without", {
  set.seed(3); s0 <- .Random.seed; mvImpute(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_identical(.Random.seed, s0)
  set.seed(3); mvImpute(dh, workers = 1, maxIter = c(2, 2)); s1 <- .Random.seed
  set.seed(3); sample.int(.Machine$integer.max, 1L); expect_identical(s1, .Random.seed)
})
test_that("observed values and types come back unchanged; no NA left", {
  out <- mvImpute(dh, seed = 1, workers = 1, maxIter = c(3, 5))
  for (nm in names(dh)) { o <- !is.na(dh[[nm]])
    expect_equal(out[[nm]][o], dh[[nm]][o]) }
  expect_false(anyNA(out)); expect_identical(class(out$ord4), c("ordered", "factor"))
  expect_type(out$chr2, "character"); expect_type(out$lgl, "logical")
})
test_that("exclude is invisible to the models", {
  # y is built from num1 so it WOULD inform the models if it leaked in; d2 has unrelated values.
  d1 <- dh; d1$y <- d0$num1 + withr::with_seed(9, rnorm(150, sd = 0.1)); d1$y[1:10] <- NA
  d2 <- d1; d2$y <- withr::with_seed(10, rnorm(150))
  a <- mvImpute(d1, exclude = "y", seed = 1, workers = 1, maxIter = c(3, 5))
  b <- mvImpute(d2, exclude = "y", seed = 1, workers = 1, maxIter = c(3, 5))
  none <- mvImpute(dh, seed = 1, workers = 1, maxIter = c(3, 5))
  expect_identical(a$y, d1$y)
  expect_identical(strip(a)[names(dh)], strip(b)[names(dh)])
  expect_identical(strip(a)[names(dh)], strip(none))
})
test_that("numeric matrix input: class, dimnames and diagnostics", {
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.05, "mcar", seed = 2)
  out <- mvImpute(X, seed = 1, workers = 1, maxIter = c(3, 5))
  expect_true(is.matrix(out)); expect_identical(dimnames(out), dimnames(X)); expect_false(anyNA(out))
  expect_named(attr(out, "mvImpute"), c("seed", "rng_kind", "method", "settings", "workers",
    "iterations", "not_converged", "time", "oob_error", "rows_no_observed"), ignore.order = TRUE)
})
test_that("forceZero floors numeric imputations at zero", {
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.05, "mcar", seed = 2); X <- log(X)  # half negative
  out <- mvImpute(X, seed = 1, workers = 1, maxIter = c(2, 2), forceZero = TRUE)
  expect_true(all(out[is.na(X)] >= 0))
})
test_that("every guess option runs and observed values stay put", {
  X <- punch_holes(sim_numeric(60, 6, seed = 1), 0.05, "mcar", seed = 2)
  for (g in list("minVar", "median", "minTotal", ifelse(is.na(X), 1, X)))
    expect_equal(mvImpute(X, seed = 1, workers = 1, maxIter = c(2, 2), guess = g)[!is.na(X)], X[!is.na(X)])
})
test_that("a row with no observed value warns and is listed", {
  X <- punch_holes(sim_numeric(40, 5, seed = 1), 0.05, "mcar", seed = 2); X[7, ] <- NA
  expect_warning(out <- mvImpute(X, seed = 1, workers = 1, maxIter = c(2, 2)), "no observed value")
  expect_identical(attr(out, "mvImpute")$rows_no_observed, 7L)
})
test_that("an absent .Random.seed stays absent", {
  withr::local_preserve_seed()
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  mvImpute(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_false(exists(".Random.seed", globalenv()))
})
test_that("NaN counts as missing; an all-zero-but-one column is fine", {
  X <- sim_numeric(40, 5, seed = 1); X[3, 2] <- NaN
  X[, 5] <- 0; X[1, 5] <- 4; X[10, 5] <- NA
  out <- mvImpute(X, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_false(anyNA(out))
})
test_that("each backend runs on mixed data and reports oob", {
  for (m in c("ranger", "randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    out <- mvImpute(dh, method = m, seed = 1, workers = 1, maxIter = c(2, 3))
    expect_false(anyNA(out))
    oob <- attr(out, "mvImpute")$oob_error
    if (m == "pls") expect_true(all(is.na(oob))) else expect_true(all(is.finite(oob)))
  }
})
test_that("imputation beats column means and modes on simulated truth", {
  out <- mvImpute(dh, seed = 1, workers = 1)
  m <- is.na(dh$num1); mean_err <- mean((mean(dh$num1, na.rm = TRUE) - d0$num1[m])^2)
  expect_lt(mean((out$num1[m] - d0$num1[m])^2), mean_err)
  f <- is.na(dh$fac3); mode_acc <- mean(names(which.max(table(dh$fac3))) == d0$fac3[f])
  expect_gt(mean(out$fac3[f] == d0$fac3[f]), mode_acc)
})
test_that("a tiny table with one factor completes without backend warnings", {
  d <- sim_mixed(20, seed = 3)[, c("num1", "num2", "fac3")]; d$num1[2] <- NA; d$fac3[5] <- NA
  # Only backend warnings are disallowed; the package's own cap warning is not what this tests.
  expect_no_warning(mvImpute(d, seed = 1, workers = 1, maxIter = c(3, 5)),
                    message = "standard deviation is zero|empty classes|five or fewer unique")
})
test_that("only the final stage warns at the cap", {
  expect_warning(mvImpute(dh, seed = 1, workers = 1, tol = c(1e-9, 1e-9), maxIter = c(1, 1)), "not converged")
  expect_length(attr(suppressWarnings(mvImpute(dh, seed = 1, workers = 1, tol = c(1e-9, 1e-9), maxIter = c(1, 1))), "mvImpute")$not_converged, 2)
})
test_that("a constant column with gaps is filled with a message and never stalls a stage", {
  d <- dh; d$const[1:3] <- NA
  expect_message(out <- mvImpute(d, seed = 1, workers = 1), "Column 'const' is constant")
  expect_true(all(out$const == 5))
  expect_equal(nrow(attr(out, "mvImpute")$not_converged[[2]][attr(out, "mvImpute")$not_converged[[2]]$column == "const", ]), 0)
})
test_that("verbose controls progress messages", {
  expect_no_message(mvImpute(dh, seed = 1, workers = 1, maxIter = c(2, 2)))
  expect_message(mvImpute(dh, seed = 1, workers = 1, maxIter = c(2, 2), verbose = TRUE), "Iteration")
})
test_that("no missing values: input returned with a message", {
  expect_message(out <- mvImpute(d0, seed = 1, workers = 1), "no missing")
  expect_identical(strip(out), d0)
})
```

- [ ] **Step 2: Run** with `filter = 'mvImpute'`. Expected: FAIL.
- [ ] **Step 3: Implement** `R/mvImpute.R` with full roxygen documentation (all arguments,
  return value, diagnostics fields, Assumptions section, one example on
  `sim`-style inline data using `workers = 1`). Run `devtools::document()`.
- [ ] **Step 4: Run** the full suite (log name `all`). Expected: all PASS.
- [ ] **Step 5: Commit** (if allowed): `git commit -m "Add mvImpute()"`.

---

### Task 9: Experiments and `tol` recalibration

Review tier: **full** (produces the numbers behind the package defaults).

**Files:**
- Create: `bench/experiments.R`, `bench/calibrate_tol.R`, `bench/output/` (gitignored except
  `*.md`), `bench/output/experiments.md`, `bench/output/calibration.md`

**Interfaces:**
- Consumes: the installed package (`devtools::install()` first) and the helpers from
  `tests/testthat/helper-sim.R` (`source()`d by path).
- Produces: two short reports for Carl. They decide the `crn` default and the `tol` defaults;
  code defaults change only after Carl signs off.

- [ ] **Step 1: `bench/experiments.R`.** (a) Common random numbers: run `run_stage()` (via
  `mvImputer:::`) with `crn = TRUE` and `crn = FALSE` on `sim_numeric(200, 300)` with 5 % MCAR
  holes, 5 seeds; report iterations per stage and RMSE against truth. (b) Two strongly
  associated factors missing in the same rows (add to `sim_mixed(300)` a factor `fac3b` equal
  to `fac3` with 10 % of labels permuted; co-missing in 40 rows): report whether the final
  stage reaches `maxIter` and the per-iteration count of changed class labels (2-cycle check).
  Write `bench/output/experiments.md`.
- [ ] **Step 2: `bench/calibrate_tol.R`.** Grid `tol1 ∈ {0.1, 0.05, 0.025}` × `tol2 ∈ {0.05, 0.025, 0.01, 0.005}`,
  scenarios (i) `sim_numeric(200, 500)` with 5 % MCAR, (ii) same with 5 % `"left"`, (iii)
  `sim_mixed(900, extra = 10)` (20 columns) with 5 % MCAR; 5 replicate seeds each; `workers = 2`.
  Old-engine reference: installed StatTools 0.0.917 (commit `7ee1dc09`; record
  `packageDescription("StatTools")$RemoteSha` in the report), `StatTools::mvImpWrap()` with its
  defaults (`tol1 = 0.05, n1 = 15, tol2 = 0.025, n2 = 60`, `method = "RF"`, `rfMeth = "ranger"`),
  `nCore = 2`, after `set.seed(<replicate seed>)` (not reproducible across runs: report means
  over replicates, not single runs). Its iterations are counted from `utils::capture.output()`
  as lines starting `Iteration`. For scenario iii it gets the full-K one-hot route of spec §9
  (`stats::model.matrix(~ f - 1)` per factor, decoded by the largest indicator). Selection rule (spec §5.4): per
  scenario, the smallest mean RMSE (numeric) whose mean total iterations do not exceed the old
  engine's. Write `bench/output/calibration.md` with one table per scenario and the proposed
  defaults. First time one replicate of each scenario and extrapolate; if the full grid exceeds
  30 minutes, run it in the background (`run_in_background`) with its log in `dev/logs/`.
- [ ] **Step 3: Stop for Carl's sign-off.** Present both reports. If Carl changes the `tol`
  defaults, update `mvImpute()`'s formals and docs; if he drops common random numbers, change
  `run_stage()`'s `crn` default to `FALSE` (it is not a public argument). Rerun the full test
  suite and record the decision in `NEWS.md`.
- [ ] **Step 4: Commit** (if allowed): `git commit -m "Add experiments and tol calibration"`.

**Checkpoint:** Carl signs off the defaults before benchmarking against the old version.

---

### Task 10: Benchmark against StatTools

Review tier: **light**.

**Files:**
- Create: `bench/benchmark.R`, `bench/output/benchmark.md`

- [ ] **Step 1:** Time `mvImpute()` and `StatTools::mvImpWrap()` at the same worker count (2 and
  4) on `sim_numeric(200, 1000)` with 5 % MCAR and 5 % left-censored holes, and on
  `sim_mixed(900)` (old engine via one-hot, spec §9). Report: time ratio new/old, RMSE (numeric)
  and accuracy (factors) against truth, and the thread note (old engine: ranger 2 threads per
  worker; new: 1). Mixed table: `sim_mixed(900, extra = 10)`. Flag a ratio above 1.25. Optional row: `missRanger::missRanger()` accuracy
  as an outside reference, only if it is installed.
- [ ] **Step 2:** Run it with output to `dev/logs/benchmark.log`; check `bench/output/benchmark.md` exists and has both scenarios.
- [ ] **Step 3: Commit** (if allowed): `git commit -m "Add benchmark against StatTools"`.

---

### Task 11: Documentation and CRAN check

Review tier: **light**.

**Files:**
- Create: `README.md`, `NEWS.md`; Modify: `DESCRIPTION` (Description field)

- [ ] **Step 1:** README: what it does, install line
  (`remotes::install_github("CarlBrunius/mvImputer")`), a minimal example, a "Moving from
  StatTools::mvImpWrap" table (old argument → new), the Assumptions section, contact (gmail and
  Carl's NBIS address). NEWS: 0.1.0 changes versus StatTools 0.0.917 (findings F1–F4, the
  training-row change, new features).
- [ ] **Step 2:** `R CMD build .` then `R CMD check --as-cran mvImputer_0.1.0.tar.gz`, output to
  `dev/logs/check.log`; `tail -n 40`. If `Sys.which("pdflatex")` is empty, add `--no-manual` and
  say so in the report. Expected: `Status: OK` or NOTEs only ("New submission" and an
  unreachable `BugReports` URL before the repo is public are fine). Fix any error or warning.
- [ ] **Step 3: Commit** (if allowed): `git commit -m "Add README, NEWS; pass R CMD check"`.

---

### Task 12: StatTools deprecation (separate repository)

Review tier: **light**. Repository: a local clone of GitLab `CarlBrunius/StatTools`. Runs only after Carl says
mvImputer is public, so the pointer resolves.

- [ ] **Step 1:** `git pull`, confirm `git log -1` is still the last commit before this change
  (it was `7ee1dc09` on 2026-11-23), then `git tag pre-mvImputer-split <that hash>` (local only).
- [ ] **Step 2:** Add `.Deprecated("mvImputer::mvImpute")` as the first line of `mvImp()` and
  `mvImpWrap()`; README rows marked deprecated with the pointer and one line on the stopping-rule
  bugs; version 0.0.918 in DESCRIPTION; README version table row.
- [ ] **Step 3:** `Rscript -e "devtools::install('<path to the StatTools clone>')"` then
  `Rscript -e "expect <- tryCatch(StatTools::mvImp(matrix(c(1,NA,3,4,5,6,7,8),4)), warning = function(w) conditionMessage(w)); cat(expect)"`.
  Expected: a message naming `mvImputer::mvImpute`.
- [ ] **Step 4:** Commit and push only on Carl's word.
