# Calibration of the default tolerances (`tol`)

## Summary

- Compared `mvImpute()` over tol1 in {0.1, 0.05, 0.025} x tol2 in {0.05, 0.025, 0.01, 0.005} with the old engine (`StatTools::mvImpWrap()`, its defaults) on three simulated scenarios, 5 replicate data sets each. The old engine ran both stages to their caps (mean total iterations (i) 75.0, (ii) 75.0, (iii) 75.0 of a possible 75), so the rule's iteration condition excluded (i) 0, (ii) 0, (iii) 0 of 12 settings and the rule reduced to the smallest mean scaled RMSE.
- Accuracy is flat over the grid: the largest mean paired difference between any setting and the scenario's pick is (i) 0.0007, (ii) 0.0016, (iii) 0.0027 scaled-RMSE units, against scaled RMSEs of about (i) 0.86, (ii) 0.49, (iii) 0.78 (old engine (i) 0.858, (ii) 0.492, (iii) 0.781). Stage 2 reached its cap in (i) 0, (ii) 0, (iii) 54 of 60 new-engine runs; in the mixed table a smaller tol does not remove the cells left.
- Proposed default: tol = c(0.05, 0.05), the one setting closest to every scenario's pick (smallest worst-case paired RMSE excess). The rule's picks per scenario: (i) tol = c(0.1, 0.005); (ii) tol = c(0.1, 0.05); (iii) tol = c(0.025, 0.025). Builder's alternative (fewest iterations within 2 standard errors of each pick): (i) tol = c(0.1, 0.05); (ii) tol = c(0.1, 0.05); (iii) tol = c(0.05, 0.01).
- Uncertainty: 5 replicates per setting, and the gaps between settings are of the size of their paired SD, so the ranking within the grid is not stable; the choice is mainly one of iteration cost. The old engine's runs are not reproducible (its parallel workers are not seeded), so its numbers are means over replicates and would move somewhat on a rerun.

## Definitions (choices for Carl to confirm)

- **Scaled RMSE** (selection metric): per replicate, the root mean squared error over all imputed cells of numeric columns, each error divided by the SD of that column's observed values (the scale the stopping rule uses), pooled over columns; then the mean over replicates. Constant columns are left out (filled exactly). **Raw RMSE** is the same without scaling. **Factor accuracy** (scenario iii) is the share of imputed factor-like cells (factor, ordered, character, logical) with the true label; it is reported, not used to select.
- **Total iterations**: stage 1 + stage 2 iterations. For the old engine, counted from the lines starting `Iteration` in its printed log.
- **Selection rule** (spec section 5.4): per scenario, the setting with the smallest mean scaled RMSE whose mean total iterations do not exceed the old engine's mean.
- **Diff. to selected**: per replicate, the setting's scaled RMSE minus the selected setting's on the same data set; mean (SD) over replicates. Positive means worse than the selected setting.
- **Builder's alternative**: among settings whose mean paired difference to the scenario's pick is at most 2 standard errors (SD / sqrt(5)), the one with the fewest mean total iterations (ties: the larger tol).
- **Stage 2 at cap**: replicates whose stage 2 ran to its iteration cap (60). **Cells left**: mean number of imputed cells still not converged at the end of stage 2 (old engine: from its log).

## Setup

- Scenarios: (i) `sim_numeric(200, 500)` (log-normal, rank-3 signal plus noise) with 5 % of cells missing completely at random (MCAR); (ii) the same data with the lowest 5 % of each column missing (left-censored, a detection-limit pattern); (iii) `sim_mixed(900, extra = 10)` (20 columns: numeric, integer, factor, ordered factor, character, logical, a constant and 10 noise columns) with 5 % MCAR over all columns.
- Replicate r (r = 1, 2, 3, 4, 5) uses data seed r, hole seed 1000 + r, and seed r for both engines (the old engine after a scoped `set.seed(r)`).
- `mvImpute()`: defaults (ranger, 100 trees, mtry 5, guess "minVar", maxIter 15 and 60) except `tol`; one 2-worker multisession pool per process.
- Old engine: StatTools 0.0.917, commit `7ee1dc09bda153de7b62a344484db53b09613516`; `mvImpWrap(method = "RF", rfMeth = "ranger", tol1 = 0.05, n1 = 15, tol2 = 0.025, n2 = 60, nCore = 2)`; first guess column minimum (its default). Its tolerance is relative to the current value, not to the column SD, so its tol values (marked *) are not on the same scale as the new ones.
- Old engine, scenario iii: character and logical columns made factors; each factor (including the ordered one) became K indicator columns (`model.matrix(~ f - 1)`, missing rows kept); imputed indicators decoded to the level with the largest value (ties to the first level).
- R version 4.6.1 (2026-06-24 ucrt); mvImputer 0.1.0; ranger 0.18.0; future 1.70.0. Sum of per-run times 638 minutes. Runs were spread over several concurrent processes, so the seconds columns include contention and are indicative only; the timing comparison belongs to the benchmark.

## Scenario i: numeric 200 x 500, 5 % MCAR

| Engine | tol1 | tol2 | Scaled RMSE (SD) | Diff. to selected (SD) | Raw RMSE | Iter. stage 1 + 2 | Total iter. (range) | Stage 2 at cap | Cells left (final stage) | Seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| old | 0.05* | 0.025* | 0.858 (0.047) |  | 2.880 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 4.6 | 695 |
| new | 0.1 | 0.05 | 0.860 (0.048) | 0.0002 (0.0008) | 2.868 | 4.8 + 5.0 | 9.8 (9-10) | 0 / 5 | 0.0 | 196 |
| new | 0.1 | 0.025 | 0.860 (0.048) | 0.0002 (0.0005) | 2.867 | 4.8 + 6.8 | 11.6 (11-12) | 0 / 5 | 0.0 | 169 |
| new | 0.1 | 0.01 | 0.860 (0.048) | 0.0003 (0.0004) | 2.869 | 4.8 + 10.8 | 15.6 (15-16) | 0 / 5 | 0.0 | 223 |
| **new (selected)** | 0.1 | 0.005 | 0.860 (0.047) | 0.0000 (0.0000) | 2.868 | 4.8 + 14.6 | 19.4 (18-21) | 0 / 5 | 0.0 | 284 |
| new | 0.05 | 0.05 | 0.860 (0.047) | 0.0003 (0.0003) | 2.867 | 6.0 + 5.0 | 11.0 (10-12) | 0 / 5 | 0.0 | 212 |
| new | 0.05 | 0.025 | 0.860 (0.047) | 0.0004 (0.0005) | 2.868 | 6.0 + 7.0 | 13.0 (13-13) | 0 / 5 | 0.0 | 188 |
| new | 0.05 | 0.01 | 0.860 (0.047) | 0.0003 (0.0004) | 2.866 | 6.0 + 11.4 | 17.4 (16-18) | 0 / 5 | 0.0 | 249 |
| new | 0.05 | 0.005 | 0.861 (0.047) | 0.0007 (0.0006) | 2.866 | 6.0 + 15.4 | 21.4 (21-23) | 0 / 5 | 0.0 | 300 |
| new | 0.025 | 0.05 | 0.860 (0.047) | 0.0001 (0.0006) | 2.869 | 8.0 + 5.0 | 13.0 (12-14) | 0 / 5 | 0.0 | 250 |
| new | 0.025 | 0.025 | 0.860 (0.047) | 0.0001 (0.0005) | 2.870 | 8.0 + 6.6 | 14.6 (13-15) | 0 / 5 | 0.0 | 180 |
| new | 0.025 | 0.01 | 0.860 (0.048) | 0.0000 (0.0005) | 2.868 | 8.0 + 11.2 | 19.2 (19-20) | 0 / 5 | 0.0 | 276 |
| new | 0.025 | 0.005 | 0.860 (0.047) | 0.0001 (0.0001) | 2.867 | 8.0 + 15.2 | 23.2 (22-25) | 0 / 5 | 0.0 | 309 |

## Scenario ii: numeric 200 x 500, 5 % left-censored

| Engine | tol1 | tol2 | Scaled RMSE (SD) | Diff. to selected (SD) | Raw RMSE | Iter. stage 1 + 2 | Total iter. (range) | Stage 2 at cap | Cells left (final stage) | Seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| old | 0.05* | 0.025* | 0.492 (0.027) |  | 0.578 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 11.4 | 746 |
| **new (selected)** | 0.1 | 0.05 | 0.493 (0.028) | 0.0000 (0.0000) | 0.590 | 4.0 + 4.4 | 8.4 (8-9) | 0 / 5 | 0.0 | 153 |
| new | 0.1 | 0.025 | 0.494 (0.028) | 0.0002 (0.0006) | 0.592 | 4.0 + 6.0 | 10.0 (10-10) | 0 / 5 | 0.0 | 147 |
| new | 0.1 | 0.01 | 0.494 (0.028) | 0.0010 (0.0010) | 0.593 | 4.0 + 9.6 | 13.6 (13-14) | 0 / 5 | 0.0 | 153 |
| new | 0.1 | 0.005 | 0.495 (0.028) | 0.0012 (0.0007) | 0.594 | 4.0 + 14.4 | 18.4 (18-19) | 0 / 5 | 0.0 | 249 |
| new | 0.05 | 0.05 | 0.494 (0.028) | 0.0003 (0.0007) | 0.592 | 5.2 + 4.4 | 9.6 (9-10) | 0 / 5 | 0.0 | 165 |
| new | 0.05 | 0.025 | 0.494 (0.028) | 0.0009 (0.0007) | 0.593 | 5.2 + 6.4 | 11.6 (11-13) | 0 / 5 | 0.0 | 149 |
| new | 0.05 | 0.01 | 0.495 (0.028) | 0.0015 (0.0007) | 0.594 | 5.2 + 9.6 | 14.8 (14-16) | 0 / 5 | 0.0 | 191 |
| new | 0.05 | 0.005 | 0.495 (0.028) | 0.0016 (0.0007) | 0.595 | 5.2 + 13.8 | 19.0 (17-22) | 0 / 5 | 0.0 | 262 |
| new | 0.025 | 0.05 | 0.494 (0.028) | 0.0007 (0.0010) | 0.593 | 7.0 + 4.2 | 11.2 (11-12) | 0 / 5 | 0.0 | 206 |
| new | 0.025 | 0.025 | 0.494 (0.028) | 0.0010 (0.0011) | 0.593 | 7.0 + 6.0 | 13.0 (13-13) | 0 / 5 | 0.0 | 166 |
| new | 0.025 | 0.01 | 0.495 (0.028) | 0.0013 (0.0010) | 0.595 | 7.0 + 10.0 | 17.0 (16-18) | 0 / 5 | 0.0 | 225 |
| new | 0.025 | 0.005 | 0.495 (0.028) | 0.0013 (0.0011) | 0.595 | 7.0 + 14.4 | 21.4 (20-23) | 0 / 5 | 0.0 | 265 |

## Scenario iii: mixed 900 x 20, 5 % MCAR

| Engine | tol1 | tol2 | Scaled RMSE (SD) | Diff. to selected (SD) | Raw RMSE | Factor accuracy | Iter. stage 1 + 2 | Total iter. (range) | Stage 2 at cap | Cells left (final stage) | Seconds |
|---|---|---|---|---|---|---|---|---|---|---|---|
| old | 0.05* | 0.025* | 0.781 (0.026) |  | 0.900 | 0.836 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 252.6 | 69 |
| new | 0.1 | 0.05 | 0.782 (0.026) | 0.0019 (0.0021) | 0.901 | 0.859 | 15.0 + 49.6 | 64.6 (23-75) | 4 / 5 | 9.2 | 69 |
| new | 0.1 | 0.025 | 0.781 (0.027) | 0.0017 (0.0016) | 0.901 | 0.859 | 15.0 + 50.4 | 65.4 (27-75) | 4 / 5 | 17.8 | 78 |
| new | 0.1 | 0.01 | 0.782 (0.027) | 0.0023 (0.0029) | 0.901 | 0.858 | 15.0 + 53.6 | 68.6 (43-75) | 4 / 5 | 25.8 | 76 |
| new | 0.1 | 0.005 | 0.782 (0.027) | 0.0027 (0.0011) | 0.902 | 0.858 | 15.0 + 57.6 | 72.6 (63-75) | 4 / 5 | 28.0 | 82 |
| new | 0.05 | 0.05 | 0.780 (0.026) | 0.0006 (0.0025) | 0.899 | 0.856 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 12.2 | 82 |
| new | 0.05 | 0.025 | 0.781 (0.026) | 0.0017 (0.0026) | 0.900 | 0.861 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 19.2 | 89 |
| new | 0.05 | 0.01 | 0.781 (0.027) | 0.0015 (0.0025) | 0.900 | 0.857 | 15.0 + 53.4 | 68.4 (42-75) | 4 / 5 | 24.2 | 78 |
| new | 0.05 | 0.005 | 0.782 (0.027) | 0.0019 (0.0030) | 0.901 | 0.862 | 15.0 + 57.2 | 72.2 (61-75) | 4 / 5 | 31.8 | 82 |
| new | 0.025 | 0.05 | 0.780 (0.026) | 0.0000 (0.0013) | 0.899 | 0.860 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 9.8 | 72 |
| **new (selected)** | 0.025 | 0.025 | 0.780 (0.026) | 0.0000 (0.0000) | 0.899 | 0.858 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 19.4 | 89 |
| new | 0.025 | 0.01 | 0.782 (0.028) | 0.0019 (0.0026) | 0.900 | 0.859 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 27.8 | 87 |
| new | 0.025 | 0.005 | 0.782 (0.026) | 0.0022 (0.0022) | 0.901 | 0.858 | 15.0 + 60.0 | 75.0 (75-75) | 5 / 5 | 37.6 | 90 |

