# mvImputer (development version)

## Default settings (decided 2026-11-24, from simulation)

* `tol = c(0.05, 0.025)` stays the default. Over a grid of 12 tolerance pairs on simulated
  omics-like data (random and left-censored gaps) and a mixed covariate table, imputation
  accuracy was flat (largest mean difference 0.003 in SD-scaled RMSE), so no other setting
  was better supported. The earlier mvImp engine ran to its iteration caps in every run, so
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
