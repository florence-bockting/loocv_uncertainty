## Replicating results
- works
- using the snakemake pipeline we can replicate the elpd results from Sivula (2025)

## Extending results to other measures
- using a gaussian dgp, we can show that the three scenarios hold as well for the other measures: mlpd, mae, ic, r2, rmse, mse, rps, srps

- S1: Models make similar predictions (ie beta_t = 0): 
    + negative correlation between target and estimate (thus, underestimation of error variance)
    + does not recover for increasing N
    + see pipeline/results/full/figs/joint_gaussian_*
- S2: Model misspecification; outlier in the data
    + small diff: 
        + normal approx. is not calibrated (irrespective of N)
    + large diff:
        + normal approx. **is not** calibrated irrespective of N:
            + holds for: elpd, ic, mlpd, r2, rmse, srps
        + normal approx. **is** calibrated for large N:
            + holds for: mae, rps, mse (very large out is bimodel - complex)
    + see pipeline/results/full/figs/calibration_gaussian_*
- S3: The number of observation is small
    + normal approx. is not calibrated for small N
    + see pipeline/results/full/figs/calibration_gaussian_*

## With respect to loo diagnostics

- can we make the same claims regarding `p_worse` and `diag_diff`?
- `p_worse` can be shown for all measures (but probably not for acc, and bacc as there the results are less clear)
- diagnostic flags:
    - `N < 100`: this seems to hold for all measures
    - `|diff| < 4`: this is elpd specific. (other cutoff value are inspected from the simulation results in pipeline/results/full/coverage_cutoff_gaussian.pdf)
- blocker:
    + results differ for binomial vs. gaussian dgp; thus results are dependent on dgp in models
    + diagnostic flag `|diff| < 4` differs per measure and dgp and is dependent of the units of y for mae, rmse, rps and mse

## Current decision:

- add `p_worse` for all measures, except acc, bacc
- add diagnostic flag with `N < 100` for all measures (except acc, bacc)
- add diagnostic flag with 
    + `|diff| < 4` for elpd
    + `|diff| < 8` for ic
    + `|diff|*n_obs < 4` mlpd
    + `|diff|*n_obs < 6` for r2
    + `|diff|*n_obs < 2` for srps
    + for all others we don't have a general rule yet
