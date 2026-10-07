# LOO Measure Calibration

Reproducing Sivula et al. (2025) and extending it beyond elpd

Sep 25, 2026 · @User

## Summary

**Are the results of Sivula et al. (2025) reproducible?** Yes. All three problematic scenarios appear in the R rebuild, including the paper's sharpest claim — that the skewness of the error distribution does not fade as the data grow. For `elpd` with two equal models, the skewness of the standardised error rises from −0.32 at n = 16 to **+2.21** at n = 1024, and its standard deviation rises from 1.6 to 3.9. More data makes the normal approximation worse, not better.

**Do the other measures behave the same way?** Under Scenario 1 and Scenario 3, yes. `mse`, `rmse`, `mae`, `r2`, `rps` and `srps` track `elpd` closely, so the paper's warning is not specific to the log score. Under Scenario 2 they diverge sharply, and every one of them fails in a different way: `r2` and `rmse` produce intervals far too narrow, while `elpd` and `srps` acquire a bias of about one standard error that an inflated interval then hides. Only `mse`, `mae` and `rps` stay honest.

Two results that are not in the paper: the binomial family is uniformly worse calibrated than the gaussian one and needs a larger true difference before the approximation works; and `acc` and `bacc` appear immune to Scenario 1, which needs checking before anyone relies on it.

## What was run

The `simulated/` experiments of the paper, rebuilt in R and driven by Snakemake. The full grid ran on Aalto Triton: 3198 jobs, 2 h 40 min wall, about 269 CPU-hours.

|  |  |
| :-- | :-- |
| Cells | 378 |
| Trials per cell | 2000 |
| Families | gaussian (252 cells), binomial (126 cells) |
| `beta_t` | 0, 0.05, 0.1, 0.2, 0.5, 1 |
| `n_obs` | 16, 32, 64, 128, 256, 512, 1024 |
| `out_dev` | 0, 20, 200 (gaussian); 0, 2, 5 (binomial, logit scale) |
| Measures | `elpd`, `mlpd`, `ic`, `mse`, `rmse`, `mae`, `r2`, `rps`, `srps`, plus `brier`, `acc`, `bacc` for binomial |

The grid maps one-to-one onto the paper's three scenarios: `beta_t` is Scenario 1 (similar predictions, 0 means the two models are identical), `out_dev` is Scenario 2 (outliers), `n_obs` is Scenario 3 (small data).

**The diagnostic.** Every trial gives an estimate of the performance difference, its LOO standard error, and the true difference measured on held-out test sets. Write

> z = (estimate − target) / se

The normal approximation the `loo` package reports is the claim that z is standard normal. So `sd(z) = 1` means the reported interval is right; `sd(z) = 2` means it is half as wide as it should be; `sd(z) < 1` means it is too wide, which is conservative but not wrong. Skewness of z measures the asymmetry the paper's higher-moment analysis targets.

Two built-in consistency checks pass exactly on the binomial cells: `ic` equals 2 × `elpd` to the last bit, and `brier`, `mse` and `rps` are the same measure for a binary outcome (max difference 4 × 10⁻¹⁷).

The poisson family is built but is not in this run. See Limitations.

## Scenario 1 reproduces, including the part that does not fade

This is the paper's central and most counter-intuitive claim, and it comes out cleanly.

**sd(z) for `elpd`, gaussian, no outliers.** Rows are `beta_t`, columns are `n_obs`. The normal approximation claims 1.

| `beta_t` | 16 | 32 | 64 | 128 | 256 | 512 | 1024 |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| **0** | 1.61 | 1.79 | 1.97 | 2.33 | 2.72 | 3.31 | **3.93** |
| 0.05 | 1.63 | 1.75 | 1.95 | 2.27 | 2.53 | 2.95 | 2.68 |
| 0.1 | 1.65 | 1.71 | 1.85 | 2.03 | 1.97 | 1.64 | 1.26 |
| 0.2 | 1.67 | 1.74 | 1.76 | 1.53 | 1.22 | 1.10 | 1.05 |
| 0.5 | 1.79 | 1.58 | 1.20 | 1.11 | 1.04 | 1.03 | 1.03 |
| 1 | 1.81 | 1.32 | 1.12 | 1.09 | 1.02 | 1.02 | **1.02** |

Read the first row against the last. When the models genuinely differ, the approximation converges to correct as n grows, exactly as asymptotic theory would suggest. When the models are equal, it **diverges**: the reported standard error is 1.6 times too small at n = 16 and about 4 times too small at n = 1024. A 95% interval from the `loo` package then covers about 66% of the time.

**Skewness of z, same cells.** This is the quantity the abstract says does not fade away.

| `beta_t` | 16 | 32 | 64 | 128 | 256 | 512 | 1024 |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| **0** | −0.32 | 0.26 | 0.58 | 1.00 | 1.30 | 1.62 | **2.21** |
| 0.1 | −0.24 | 0.24 | 0.69 | 1.04 | 1.88 | 2.53 | 1.94 |
| 0.5 | 0.48 | 1.14 | 0.41 | 0.14 | 0.07 | 0.08 | −0.01 |
| 1 | 1.69 | 0.65 | 0.25 | 0.18 | 0.08 | 0.12 | **−0.01** |

Along the bottom row the skewness decays to zero, which is the ordinary central-limit behaviour. Along the top row it grows monotonically. The binomial family does the same thing, reaching 2.24 at n = 1024.

[PIT histograms for elpd, gaussian](pipeline/results/full/figs/calibration_gaussian_elpd.pdf)

PIT histograms for `elpd`, gaussian. Columns are `n_obs`, row pairs are `beta_t`. The shaded band is the uniform reference. Top row (`beta_t` = 0, no outlier): a hard spike at the left edge that does not soften as n grows across the row. Fifth row (`beta_t` = 1, no outlier): near uniform everywhere. Bottom row (`beta_t` = 1, outlier): the distribution collapses to a single narrow spike away from the centre, which is the bias described in the Scenario 2 section.

**The practical reading.** Collecting more data does not rescue a comparison between two models that predict alike. It makes the reported uncertainty more wrong and more asymmetric. This is the result a `loo` user is most likely to get wrong, because the instinct is that large n makes normal approximations safe.

## Scenarios 2 and 3 reproduce, and behave very differently from each other

**Scenario 3 (small data) is real but curable.** Hold `beta_t >= 0.5` so Scenario 1 is switched off, and the small-sample inflation stands alone: sd(z) is 1.4 to 1.8 at n = 16 for every measure, and 1.02 at n = 1024. It is a finite-sample effect that more data fixes, which is the opposite of Scenario 1. The two are separable in the grid, and the paper's framing of them as distinct problems holds up.

**Scenario 2 (outliers) reproduces, and it is where the measures stop agreeing.** With `beta_t = 1` and `n >= 256`, so Scenarios 1 and 3 are both off:

| measure | mean(z) @20 | sd(z) @20 | mean(z) @200 | sd(z) @200 | coverage @200 |
| :-- | --: | --: | --: | --: | --: |
| `r2` | −0.01 | 0.35 | **2.16** | **15.97** | **0.66** |
| `rmse` | −0.07 | 0.63 | −0.45 | 2.32 | 0.69 |
| `elpd`, `mlpd`, `ic` | 0.39 | 0.55 | **1.00** | 0.56 | 0.97 |
| `srps` | 0.31 | 0.83 | 0.74 | 0.51 | 0.97 |
| `mae` | 0.06 | 1.03 | 0.21 | 1.09 | 0.92 |
| `rps` | 0.09 | 1.02 | 0.28 | 1.01 | 0.95 |
| `mse` | 0.09 | 1.01 | 0.05 | 0.98 | 1.00 |

At `out_dev = 0` every measure sits at mean(z) ≈ 0.1 and sd(z) = 1.02, so all of this is the outlier's doing. Three distinct failures:

- **`r2` and `rmse` produce intervals far too narrow.** At the extreme level `r2` has sd(z) = 16 and skewness 11.2; its coverage of a nominal 95% interval falls to 0.66, and `rmse` to 0.69. A single extreme observation dominates the ratio `r2` is built from. These are the measures that actively mislead.
- **`elpd` and `srps` become biased, and the bias is masked.** At `out_dev = 200` the `elpd` estimate sits a **full standard error** away from the target, yet coverage stays at 0.97 because the interval is simultaneously about 1.8 times too wide. Two errors cancel. A user reading only the interval sees nothing wrong; the PIT histogram shows it immediately, as a narrow spike away from the centre rather than a hump at 0.5.
- **`mse`, `mae` and `rps` are essentially unaffected**, at both outlier levels.

This is the clearest new result in the run. The paper argues Scenario 2 should be caught by model checking before LOO-CV is used at all, so it does not dwell on it. But if a misspecified comparison does reach `loo`, **the choice of measure decides whether the user is misled, and in which direction.**

Note that the calibration figures plot `out_dev = 20`; the extreme columns above are `out_dev = 200`, which is in the grid but not in the current figure selection. Adding it to `plot.out_dev` in `config/full.yaml` would make the figures show the `r2` collapse.

[PIT histograms for r2, gaussian](pipeline/results/full/figs/calibration_gaussian_r2.pdf)

PIT histograms for `r2`, gaussian, at `out_dev` = 20. The bottom row shows the opposite failure to the numbers above: a central hump, meaning intervals far too wide, matching sd(z) = 0.35 at this outlier level. At `out_dev` = 200 it reverses into the collapse with sd(z) = 16. The figures do not yet plot that level.

## Do the other measures behave like elpd?

sd(z), gaussian family. Each column isolates one scenario with the others switched off. The S2 column is `out_dev = 200` at `beta_t = 1` and `n >= 256`; see the previous section for the bias that sits behind those numbers.

| measure | S1 `beta_t`=0 | S1 `beta_t`=0.1 | S2 outliers | S3 n=16 | S3 n=1024 |
| :-- | --: | --: | --: | --: | --: |
| `elpd`, `mlpd`, `ic` | 3.31 | 1.64 | 0.56 | 1.79 | 1.02 |
| `srps` | 3.38 | 1.66 | 0.51 | 1.65 | 1.02 |
| `rps` | 3.01 | 1.55 | 1.01 | 1.59 | 1.02 |
| `mae` | 2.89 | 1.48 | 1.09 | 1.38 | 1.02 |
| `rmse` | 2.77 | 1.50 | 2.32 | 1.63 | 1.02 |
| `mse` | 2.59 | 1.47 | 0.98 | 1.73 | 1.02 |
| `r2` | 2.58 | 1.47 | 15.97 | 1.73 | 1.02 |

**Scenario 1 and Scenario 3: yes, they behave alike.** Under Scenario 1 every measure sits between 2.58 and 3.38, and the coverage of a nominal 95% interval falls to 0.67–0.70 for all of them. Under Scenario 3 they all recover to 1.02. The paper's conclusions transfer to every measure in the `loo` package, not just the log score.

There is a consistent but small ordering. The log-score family (`elpd`, `mlpd`, `ic`, `srps`) is the worst under Scenario 1; the squared-error family (`mse`, `r2`) is the least bad. The gap is about 30% in sd(z) — real, reproducible across every cell, but not large enough to recommend one measure over another on these grounds.

**Scenario 2: no.** This is the single place where the measures separate, by more than an order of magnitude. See the previous section.

`mlpd` and `ic` are exact rescalings of `elpd` and land on identical rows throughout, which is a useful check that the pipeline is doing what it claims.

## The binomial family is worse, and needs a larger true difference

The paper works with linear, hierarchical linear, latent linear and spline models. The binomial family is new here: logistic regression, an exact LOO predictive from n penalised-IRLS refits, and a trapezoid quadrature over the linear predictor.

It is worse calibrated than the gaussian family everywhere. Median coverage of a nominal 95% interval is **0.72 against 0.85**.

The sharper difference is where the recovery boundary sits. For `elpd` with `beta_t = 0.1`, sd(z) at n = 1024:

| family | n=16 | n=64 | n=256 | n=1024 |
| :-- | --: | --: | --: | --: |
| gaussian | 1.65 | 1.85 | 1.97 | **1.26** |
| binomial | 1.72 | 1.98 | 2.59 | **2.89** |

The gaussian family has begun to recover by n = 1024 at this effect size. The binomial family has not — it is still diverging. In other words, a binary outcome needs a **larger true difference between the models** before the normal approximation starts to work at all. For practitioners comparing logistic models, the region where `loo` intervals are trustworthy is smaller than the linear case the paper studies.

The binomial outlier grid does almost nothing: sd(z) stays at 1.01 across `out_dev` = 0, 2, 5. That is expected, because `out_dev` shifts the linear predictor on the logit scale and the inverse logit saturates. Scenario 2 for a binary outcome would need a different mechanism — label noise, or a covariate outlier — and that is worth deciding before the paper draft.

[PIT histograms for elpd, binomial](pipeline/results/full/figs/calibration_binomial_elpd.pdf)

PIT histograms for `elpd`, binomial. Compare the top rows with the gaussian figure above: the left-edge spike at `beta_t` = 0 is at least as severe, and the `beta_t` = 0.2 row is still clearly non-uniform at n = 512 where the gaussian family has begun to settle.

## acc and bacc look immune to Scenario 1. Treat that with suspicion.

Accuracy and balanced accuracy are the one striking exception in the whole run. Under Scenario 1 at n >= 256, where every other measure sits at sd(z) ≈ 3.5:

| measure | S1 `beta_t`=0 | S1 `beta_t`=0.1 | S3 n=16 | S3 n=1024 |
| :-- | --: | --: | --: | --: |
| `elpd` and the rest | 3.54 | 2.89 | 1.74 | 1.03 |
| `acc` | **1.02** | **1.04** | 1.28 | 1.00 |
| `bacc` | **1.02** | **1.04** | 1.40 | 1.00 |

If this holds, it is a substantive finding: a coarse 0/1 loss would escape the failure mode that defeats every continuous measure, precisely because it discards the fine-grained information that couples the LOO estimate to its own target.

Three reasons not to claim it yet.

1. **Trials are missing.** `acc` and `bacc` produce no finite standard error when the LOO accuracy is constant across folds, and those trials are dropped before the moments are computed. At n = 16 that is 13.9% of trials; at n = 1024 it is 0.4%. The surviving sample is selected, and it is selected on exactly the quantity being measured. The Scenario 1 columns above use n >= 256, where the drop rate is 0.4–2.7%, so selection is mild there — but it is not zero.
2. **The distribution is not normal even where sd(z) is 1.** At n = 16, skewness is −0.98 for `acc`. The nominal 95% coverage of 0.94 there is an accident: an over-wide spread of 1.28 cancels a strong left skew. A calibrated standard deviation is not a calibrated interval.
3. **Discreteness has not been separated from the statistics.** With few distinct attainable values, z is lumpy, and the usual moment diagnostics are weak. A PIT histogram will look wrong for reasons that have nothing to do with the normal approximation.

**Suggested next step:** re-run the `acc` and `bacc` cells with the dropped trials retained under an explicit convention for a zero standard error, and compare. If sd(z) stays near 1, the finding is real and worth its own section in the paper.

[PIT histograms for acc, binomial](pipeline/results/full/figs/calibration_binomial_acc.pdf)

PIT histograms for `acc`, binomial. Every panel is close to uniform, including the `beta_t` = 0 rows where `elpd` spikes hard. The figure prints the dropped-trial count per panel, from 245 of 2000 at n = 32 down to 10 at n = 512, so the selection concern above is visible directly in the plot.

## Limitations

**The poisson family is missing.** It is implemented and passes its checks, but it is not in this run. The test-set score is too slow: for one cell at n = 16, fitting costs 0.0 s per trial and the LOO score 0.03 s, while the test-set score costs 14.7 s. The cause is identified — the count grid is sized from the largest predictive mean across the whole pooled test set, so one extreme row forces a 2048-wide grid onto all 800 of them. The fix is a per-block grid. A poisson run adds a third shape of predictive distribution and should follow.

**These are moments, not PIT histograms.** Everything above is sd(z), skewness and interval coverage. The paper's figures are PIT histograms, and a spike-shaped PIT can hide behind a normal-looking coverage number. The 213 figures from the run include the calibration plots and should be read against these tables before anything is written up.

**The target carries Monte Carlo error.** It is a mean over held-out test sets, so it is not the exact expected performance. With 250 test sets per binomial cell that error is about 6% of the LOO standard error, inflating sd(z) by roughly 0.2%. It does not explain the gaussian–binomial gap.

**One DGP per family.** The paper covers hierarchical, latent linear and spline models. This grid covers one linear DGP per family, with a single outlier observation as the misspecification mechanism.

**Scenario 2 is not really tested for binomial.** As noted above, `out_dev` on the logit scale saturates. The binomial Scenario 2 column is uninformative rather than reassuring.

### Reproducing this

Branch `snakemake-measures`. The full grid without poisson:

```
awk -F'\t' 'NR==1 || $2!="poisson"' config/cells_full.tsv > config/cells_nopois.tsv
sbatch snakemake_slurm.sh --config cells=config/cells_nopois.tsv
```

All numbers in this document come from `results/full/trials.rds` — 7.56 million rows, one per trial, cell and measure.
