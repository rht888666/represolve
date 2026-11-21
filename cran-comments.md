## R CMD check results

Tested locally on Windows 11 x64 (R version 4.6.1, 2026-06-24 ucrt) with:

```
R CMD build .
R CMD check --as-cran --no-manual represolve_1.0.3.tar.gz
```

Result: `0 errors | 0 warnings | 1 note` (standard `New submission` note only).
Without `--as-cran`, `R CMD check --no-manual represolve_1.0.3.tar.gz` returned
`Status: OK` (`0 errors | 0 warnings | 0 notes`).

## Package scope

New CRAN submission. Pure R (`NeedsCompilation: no`, `Imports: stats, utils`),
providing fold-frozen masked-completion benchmarking against independent
replicate measurement partitions, prediction-time regime routing
(`obs_agreement`, `d_cos10`, hard Table-3 router, RBF Soft Regime Router, and
adaptive k-NN blending), paired perturbation bootstrap CIs, and fixed-size
permutation-null calibration utilities companion to the Python `represolve`
package.
