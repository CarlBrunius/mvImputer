# Detection-limit simulation: which settings impute low values best

## Summary

- Simulated data (200 rows x 100 columns, log-normal, 5 data sets), 5 % of cells missing, imputed on the log scale. **Strictly the lowest values missing (`left`):** the column-minimum fill was closest (0.621), then PLS (1.131), then ranger (1.544); PLS (2 components) beat the minimum fill in 0 of 5 seeds.
- **Low values more likely, not only the lowest, missing (`graded`, a soft limit):** PLS was closest (0.924, better than the minimum fill in 5 of 5 seeds); minimum fill 1.120; ranger 1.171.
- **Half lowest-in-column, half at random (`mixed`):** PLS 0.918, ranger 1.190, minimum fill 1.675. The minimum fill is good on the lowest cells (0.554) and poor on the random ones (2.303); PLS is the reverse (1.140 and 0.621).
- **Scale and forceZero (pattern `left`):** on the original scale ranger gave original-scale RMSE 0.469 against 0.360 on the log scale; PLS on original-scale intensities gave 0.698 with `forceZero = TRUE` (54 % of imputed cells set to exactly zero) and 1.534 with `forceZero = FALSE` (59 % zero or negative). Log-scale PLS was 0.318. Log-transforming first lowered the original-scale error of both methods, and PLS on untransformed intensities gave many zero or negative values.
- **Forests never go below the observed minimum** (share of imputed cells below it in `left`: 0 % with ranger on the original scale, 0 % on the log scale): a forest predicts averages of observed values. PLS on log data did go below it (26 % of its imputed cells in `left`).
- **Speed (one worker):** ranger took 4.5 to 5.6 times as long as PLS (15.3 s against 3.4 s in `left`).
- **PLS components and `guess`:** 2 and 3 components differ little (log RMSE 0.924 against 0.953 in `graded`); `guess = "median"` against `"minVar"` (default) with 3 components: 0.954 against 0.953 in `graded`.
- Only simulated, low-rank data, 5 data sets and the settings listed here: see Caveats. These numbers say what happened in this simulation, not what will happen on a given real data set.

## Definitions

- **Log-scale RMSE**: root mean squared error over the imputed cells of the log of the imputed value minus the log of the true value, each error divided by the SD of that column's observed log values, pooled over columns. Imputations made on the log scale are exponentiated and logged again, so the same cells are compared. Not defined (n/a) when an imputed value is zero or negative.
- **Original-scale RMSE**: the same with the values themselves, each error divided by the SD of that column's observed values. It is dominated by the largest values, which are rarely the missing ones.
- **Above truth**: share of imputed cells whose imputed value is larger than the true value. **Below observed minimum**: share of imputed cells below the smallest observed value of their column. **Zero or negative**: share of imputed cells at or below zero after back-transformation.
- **Mean error, original scale**: mean of the same scaled errors, with sign (positive: imputed too high).
- **Lowest cells / random cells** (`mixed` only): the log-scale RMSE over the cells that are lowest-in-column by construction, and over the cells removed at random.
- **Cells not converged**: imputed cells still changing after the iteration cap in the final stage (mean over seeds). **Seconds**: elapsed time of the `mvImpute()` call, one worker, no worker pool; mean over seeds.
- **Minimum fill**: each gap set to the smallest observed value of its column.

## Setup

- Data: `sim_numeric(200, 100, seed = s)`, s = 1 to 5 (log-normal, rank-3 signal plus noise). Gaps: 5 % of cells. `left`: the lowest 5 % of every column (deterministic). `graded`: in each column 5 % of rows drawn with weight exp(-2 z), z the standardised log value (hole seed 3000 + s). `mixed`: the lowest 2.5 % of each column plus 2.5 % of the remaining cells at random (hole seeds 1000 + s and 2000 + s).
- `mvImpute()` with its defaults (100 trees, mtry 5, `tol = c(0.05, 0.025)`, `maxIter = c(15, 60)`, `guess = "minVar"`) except what the row says; `seed = s`, `workers = 1`. Log-scale rows impute `log(table)` and exponentiate.
- R version 4.6.1 (2026-06-24 ucrt); mvImputer 0.1.0; ranger 0.18.0; mixOmics 6.36.0; 24 logical cores. One machine; other jobs may have run at the same time, so seconds are indicative.
- Script: `bench/detection_limit.R`. It replaces two exploratory scripts (not committed) whose designs it keeps; their log-scale `left` fits are not run twice.

## Experiment 1: scale, engine and forceZero (pattern `left`)

Log-scale RMSE is n/a where an imputed value was zero or negative (log undefined); the original-scale RMSE then still applies.

| Method (pattern left) | Log-scale RMSE (SD over seeds) | Original-scale RMSE | Mean error, original scale | Above truth | Below observed min. | Zero or negative | Cells not converged | Seconds |
|---|---|---|---|---|---|---|---|---|
| min fill | 0.621 (0.031) | 0.065 | 0.042 | 100 % | 0 % | 0 % | 0.0 | - |
| ranger, original | 1.878 (0.060) | 0.469 | 0.358 | 100 % | 0 % | 0 % | 0.0 | 12.7 |
| pls 2 comp, original, forceZero = TRUE | n/a | 0.698 | 0.205 | 39 % | 63 % | 54 % | 0.0 | 4.4 |
| pls 2 comp, original, forceZero = FALSE | n/a | 1.534 | -0.419 | 35 % | 66 % | 59 % | 0.0 | 8.8 |
| ranger, log | 1.544 (0.050) | 0.360 | 0.260 | 100 % | 0 % | 0 % | 0.0 | 15.3 |
| pls 2 comp, log | 1.131 (0.062) | 0.318 | 0.181 | 93 % | 26 % | 0 % | 0.0 | 3.4 |

## Experiment 2: three missingness patterns, log scale

| Pattern | Method (log scale) | Log-scale RMSE (SD over seeds) | of which lowest cells | of which random cells | Original-scale RMSE | Above truth | Below observed min. | Iterations | Cells not converged | Seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| left | min fill | 0.621 (0.031) | - | - | 0.065 | 100 % | 0 % | - | 0.0 | - |
| left | ranger, log | 1.544 (0.050) | - | - | 0.360 | 100 % | 0 % | 12.4 | 0.0 | 15.3 |
| left | pls 2 comp, log | 1.131 (0.062) | - | - | 0.318 | 93 % | 26 % | 6.4 | 0.0 | 3.4 |
| left | pls 3 comp, log | 1.154 (0.057) | - | - | 0.339 | 92 % | 27 % | 7.8 | 0.0 | 4.5 |
| left | pls 3 comp, log, guess median | 1.155 (0.056) | - | - | 0.339 | 92 % | 27 % | 9.6 | 0.0 | 5.8 |
| graded | min fill | 1.120 (0.027) | - | - | 0.240 | 21 % | 0 % | - | 0.0 | - |
| graded | ranger, log | 1.171 (0.038) | - | - | 0.325 | 93 % | 0 % | 12.6 | 0.0 | 16.9 |
| graded | pls 2 comp, log | 0.924 (0.045) | - | - | 0.310 | 84 % | 5 % | 5.6 | 0.0 | 3.0 |
| graded | pls 3 comp, log | 0.953 (0.043) | - | - | 0.340 | 83 % | 6 % | 7.2 | 0.0 | 4.3 |
| graded | pls 3 comp, log, guess median | 0.954 (0.043) | - | - | 0.340 | 83 % | 6 % | 7.6 | 0.0 | 4.6 |
| mixed | min fill | 1.675 (0.020) | 0.554 | 2.303 | 0.998 | 50 % | 0 % | - | 0.0 | - |
| mixed | ranger, log | 1.190 (0.042) | 1.553 | 0.648 | 0.687 | 75 % | 0 % | 12.4 | 0.0 | 17.1 |
| mixed | pls 2 comp, log | 0.918 (0.047) | 1.140 | 0.621 | 0.564 | 72 % | 11 % | 6.0 | 0.0 | 3.4 |
| mixed | pls 3 comp, log | 0.944 (0.041) | 1.158 | 0.665 | 0.591 | 72 % | 12 % | 7.0 | 0.0 | 4.2 |
| mixed | pls 3 comp, log, guess median | 0.945 (0.041) | 1.159 | 0.665 | 0.591 | 72 % | 12 % | 7.8 | 0.0 | 4.9 |

## Seeds (of 5) where the method beat the minimum fill (log-scale RMSE)

| Method (log scale) | left | graded | mixed |
|---|---|---|---|
| ranger, log | 0 of 5 | 1 of 5 | 5 of 5 |
| pls 2 comp, log | 0 of 5 | 5 of 5 | 5 of 5 |
| pls 3 comp, log | 0 of 5 | 5 of 5 | 5 of 5 |
| pls 3 comp, log, guess median | 0 of 5 | 5 of 5 | 5 of 5 |

## Caveats

- The data are simulated, low-rank (rank 3), log-normal. Their columns are strongly linearly related on the log scale, which favours PLS, a linear method. Real data with weaker or non-linear structure may favour the forest.
- `left` (exactly the lowest 5 % of every column) is the most favourable case for a minimum fill: the missing values really are the lowest ones. Real detection limits are rarely this clean.
- Only 2 and 3 PLS components were tried, and ranger is the only forest backend (no `randomForest` runs here). Other tuning values (trees, mtry, tol) were not varied.
- Five data sets per pattern, one simulated structure, one machine. The SD over seeds is shown for the headline error; a difference of a few hundredths is within it.
- Seconds are for one serial run (`workers = 1`) and include no worker pool; with several workers the ratio between methods may differ. Other jobs may have run on the machine.
- The simulation knows the truth, so these RMSEs cannot be computed on real data; there, the check is how imputed values behave (for example whether they fall below the observed range), not a distance to the truth.

