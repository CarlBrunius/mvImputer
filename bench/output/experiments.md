# Experiments: common random numbers and a co-missing factor pair

## Summary

- Compared the stage engine with common random numbers (crn, the current default: a column reuses one random stream in every iteration of a stage) and without (a fresh sub-stream per iteration), on numeric data (a) and on two associated factors missing in the same rows (b); 5 seeds in (a), 20 in (b), tol = (0.05, 0.025), maxIter = (15, 60), ranger, 2 workers.
- (a) Mean total iterations: 13.0 with crn, 73.8 without. Paired difference in scaled RMSE (crn minus no crn): mean 0.0031, SD over seeds 0.0014.
- (b) Final stage at its cap: 1 of 20 runs with crn, 6 of 20 without; flip-backs (2-cycle signature) summed over runs: 144 with crn, 1935 without.
- Rule verdict (spec 5.5, "kept unless it does harm"; harm = more total iterations on average, a scaled RMSE worse by more than 2 standard errors of the paired difference, or the final stage at its cap in (b) more often than without crn): harm found in rmse. Size: with crn the scaled RMSE is 0.38 % higher (worse in 5 of 5 seeds) at 5.7 times fewer total iterations. The RMSE gap is confounded with iteration count (13.0 against 73.8 iterations): the extra iterations, rather than the random streams, may explain it (not tested).
- **Carl's decision (2026-11-24): keep crn = TRUE.** Uncertainty: 5 seeds in (a), 20 in (b); seconds were measured while other jobs ran on the machine and are inflated.

## Definitions

- **Scaled RMSE**: root mean squared error over all imputed cells against the simulated truth, each error divided by the SD of its column's observed values (the stopping rule's scale), pooled over columns. **Raw RMSE**: the same without scaling.
- **Labels changed**: imputed labels of the two factors (both columns together) that differ from the previous iteration. **Flip-backs**: labels that changed and now equal their value two iterations earlier; a 2-cycle shows as a steady count of flip-backs and a final stage that runs to its cap. **Flip-backs, last iteration**: flip-backs in the last iteration of stage 2 (non-zero means a cycle still running when the stage stopped). Sequences are run-length coded: `2x13` is 2 in 13 consecutive iterations.

## (a) Common random numbers

Data: `sim_numeric(200, 300)` (log-normal, rank-3 signal plus noise), 5 % of cells missing completely at random; seed s uses data seed s, hole seed 1000 + s and engine seed s.

| seed | crn | iter. stage 1 | iter. stage 2 | final-stage cells not converged | scaled RMSE | raw RMSE | seconds |
|---|---|---|---|---|---|---|---|
| 1 | TRUE |  6 |  7 | 0 | 0.9358 | 3.0605 | 44 |
| 1 | FALSE | 15 | 54 | 0 | 0.9314 | 3.0133 | 135 |
| 2 | TRUE |  6 |  7 | 0 | 0.7175 | 1.6859 | 42 |
| 2 | FALSE | 15 | 60 | 2 | 0.7155 | 1.6763 | 145 |
| 3 | TRUE |  6 |  8 | 0 | 0.7721 | 1.8473 | 43 |
| 3 | FALSE | 15 | 60 | 3 | 0.7699 | 1.8057 | 143 |
| 4 | TRUE |  6 |  6 | 0 | 0.8487 | 2.2209 | 43 |
| 4 | FALSE | 15 | 60 | 2 | 0.8437 | 2.1762 | 139 |
| 5 | TRUE |  6 |  7 | 0 | 0.8010 | 1.9826 | 42 |
| 5 | FALSE | 15 | 60 | 1 | 0.7989 | 1.9315 | 140 |

Means over seeds:

| crn | iter. stage 1 | iter. stage 2 | final-stage cells not converged | scaled RMSE | raw RMSE | seconds |
|---|---|---|---|---|---|---|
| FALSE | 15.0 | 58.8 | 1.6 | 0.8119 | 2.1206 | 140 |
| TRUE | 6.0 | 7.0 | 0.0 | 0.8150 | 2.1595 | 43 |

Without crn the engine rebuilds each column's sub-stream from the start of the stage in every iteration (cost grows with the square of the iteration count); the seconds columns show what that costs at these sizes.

## (b) Two associated factors missing in the same rows

Data: `sim_mixed(300)` plus `fac3b`, a copy of `fac3` whose labels are permuted within a random 10 % of rows; both missing in the same 40 random rows (no other gaps). The table after iteration k is recomputed as the engine run with maxIter = k, which is exact (the first k iterations do not depend on maxIter).

Reading the crn rows: the two factors are each other's only imputed predictors and their training rows hold no imputed cells, so under crn each factor's model is the same in every iteration and an iteration is a fixed, deterministic map. A change therefore shows twice (`Nx2`: a label changes, then its partner follows) and then stops (`0x2`: two unchanged predictions freeze a cell). Seed 6's 2-cycle with crn is a genuine cycle of that map, not Monte Carlo noise.

| seed | crn | fac3 = fac3b | iter. 1 + 2 | stage 2 at cap | labels changed per iteration, stage 1 | stage 2 | flip-backs | flip-backs, last iteration | accuracy fac3 / fac3b |
|---|---|---|---|---|---|---|---|---|---|
|  1 | TRUE | 0.94 | 4 + 2 | FALSE | 11x2 0x2 | 0x2 |   0 | 0 | 0.52 / 0.55 |
|  1 | FALSE | 0.94 | 15 + 60 | TRUE | 11 15 7 5 4x11 | 4x32 3 2x27 | 240 | 2 | 0.52 / 0.52 |
|  2 | TRUE | 0.94 | 6 + 2 | FALSE | 9x2 1x2 0x2 | 0x2 |   0 | 0 | 0.55 / 0.57 |
|  2 | FALSE | 0.94 | 8 + 2 | FALSE | 9 12 6 4 2 1 0x2 | 0x2 |  12 | 0 | 0.55 / 0.57 |
|  3 | TRUE | 0.93 | 4 + 4 | FALSE | 9x2 0x2 | 2x2 0x2 |   0 | 0 | 0.68 / 0.60 |
|  3 | FALSE | 0.93 | 7 + 14 | FALSE | 9 11 2x2 1 0x2 | 4 6 5 3 2x7 1 0x2 |  32 | 0 | 0.70 / 0.62 |
|  4 | TRUE | 0.94 | 6 + 4 | FALSE | 12x2 2x2 0x2 | 2x2 0x2 |   0 | 0 | 0.62 / 0.60 |
|  4 | FALSE | 0.94 | 8 + 34 | FALSE | 12 15 7 5 2 1 0x2 | 2 4 3 2x28 1 0x2 |  73 | 0 | 0.62 / 0.60 |
|  5 | TRUE | 0.95 | 4 + 2 | FALSE | 7x2 0x2 | 0x2 |   0 | 0 | 0.45 / 0.42 |
|  5 | FALSE | 0.95 | 15 + 5 | FALSE | 7 8 2x11 1 0 | 0 1x2 0x2 |  25 | 0 | 0.45 / 0.42 |
|  6 | TRUE | 0.93 | 15 + 60 | TRUE | 10 8 2x13 | 4x2 2x58 | 144 | 2 | 0.65 / 0.70 |
|  6 | FALSE | 0.93 | 15 + 60 | TRUE | 10 9 2x13 | 5 6x59 | 375 | 6 | 0.68 / 0.72 |
|  7 | TRUE | 0.95 | 4 + 4 | FALSE | 13x2 0x2 | 2x2 0x2 |   0 | 0 | 0.62 / 0.62 |
|  7 | FALSE | 0.95 | 6 + 24 | FALSE | 13 14 0 1 0x2 | 2 3 2x19 1 0x2 |  40 | 0 | 0.62 / 0.62 |
|  8 | TRUE | 0.94 | 4 + 4 | FALSE | 12x2 0x2 | 1x2 0x2 |   0 | 0 | 0.85 / 0.82 |
|  8 | FALSE | 0.94 | 4 + 9 | FALSE | 12x2 0x2 | 1 2x5 1 0x2 |  10 | 0 | 0.85 / 0.82 |
|  9 | TRUE | 0.94 | 4 + 2 | FALSE | 12x2 0x2 | 0x2 |   0 | 0 | 0.55 / 0.62 |
|  9 | FALSE | 0.94 | 8 + 2 | FALSE | 12 13 2x3 1 0x2 | 0x2 |   8 | 0 | 0.55 / 0.62 |
| 10 | TRUE | 0.93 | 4 + 2 | FALSE | 15x2 0x2 | 0x2 |   0 | 0 | 0.47 / 0.50 |
| 10 | FALSE | 0.93 | 8 + 2 | FALSE | 15 18 5 3 2 1 0x2 | 0x2 |  11 | 0 | 0.47 / 0.50 |
| 11 | TRUE | 0.93 | 4 + 2 | FALSE | 13x2 0x2 | 0x2 |   0 | 0 | 0.65 / 0.60 |
| 11 | FALSE | 0.93 | 15 + 60 | TRUE | 13 17 8 6 4 3 2x9 | 2x60 | 159 | 2 | 0.65 / 0.57 |
| 12 | TRUE | 0.92 | 4 + 2 | FALSE | 14x2 0x2 | 0x2 |   0 | 0 | 0.52 / 0.57 |
| 12 | FALSE | 0.92 | 15 + 60 | TRUE | 14 20 11 9 8 7 6x5 5 4x3 | 4 5x2 4x57 | 320 | 4 | 0.52 / 0.57 |
| 13 | TRUE | 0.92 | 4 + 2 | FALSE | 11x2 0x2 | 0x2 |   0 | 0 | 0.75 / 0.72 |
| 13 | FALSE | 0.92 | 12 + 2 | FALSE | 11 13 4 3 2x5 1 0x2 | 0x2 |  18 | 0 | 0.75 / 0.72 |
| 14 | TRUE | 0.94 | 4 + 4 | FALSE | 11x2 0x2 | 1x2 0x2 |   0 | 0 | 0.55 / 0.55 |
| 14 | FALSE | 0.94 | 6 + 6 | FALSE | 11 12 2 1 0x2 | 1 2x2 1 0x2 |   8 | 0 | 0.55 / 0.55 |
| 15 | TRUE | 0.93 | 4 + 4 | FALSE | 10x2 0x2 | 7x2 0x2 |   0 | 0 | 0.72 / 0.60 |
| 15 | FALSE | 0.93 | 15 + 60 | TRUE | 10 17 5x2 2 4 3x4 2 3x2 2x2 | 11 8x2 6 4x2 3 2x53 | 170 | 2 | 0.75 / 0.62 |
| 16 | TRUE | 0.93 | 4 + 2 | FALSE | 12x2 0x2 | 0x2 |   0 | 0 | 0.60 / 0.55 |
| 16 | FALSE | 0.93 | 10 + 7 | FALSE | 12 13 2x5 1 0x2 | 0 1 0 1x2 0x2 |  13 | 0 | 0.60 / 0.55 |
| 17 | TRUE | 0.94 | 4 + 4 | FALSE | 8x2 0x2 | 1x2 0x2 |   0 | 0 | 0.52 / 0.45 |
| 17 | FALSE | 0.94 | 15 + 29 | FALSE | 8 11 6 5 4x3 3 2x7 | 3 4x4 3 2x20 1 0x2 |  98 | 0 | 0.55 / 0.47 |
| 18 | TRUE | 0.93 | 4 + 4 | FALSE | 10x2 0x2 | 1x2 0x2 |   0 | 0 | 0.60 / 0.55 |
| 18 | FALSE | 0.93 | 6 + 8 | FALSE | 10 11 0 1 0x2 | 1 3x2 2x2 1 0x2 |   9 | 0 | 0.60 / 0.55 |
| 19 | TRUE | 0.92 | 4 + 2 | FALSE | 10x2 0x2 | 0x2 |   0 | 0 | 0.62 / 0.60 |
| 19 | FALSE | 0.92 | 6 + 2 | FALSE | 10 11 2 1 0x2 | 0x2 |   4 | 0 | 0.62 / 0.60 |
| 20 | TRUE | 0.93 | 4 + 2 | FALSE | 14x2 0x2 | 0x2 |   0 | 0 | 0.72 / 0.68 |
| 20 | FALSE | 0.93 | 15 + 60 | TRUE | 14 19 10 9 8 7 6 5 4x7 | 4x60 | 310 | 4 | 0.72 / 0.72 |

## Setup

- R version 4.6.1 (2026-06-24 ucrt); mvImputer 0.1.0; ranger 0.18.0; future 1.70.0. One 2-worker multisession pool for all runs. Total run time 34 minutes.
- The engine is driven through its internal stage function so crn can be switched; a check in the script confirms that with crn it returns exactly what `mvImpute()` returns.

