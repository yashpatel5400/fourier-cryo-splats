# Descriptive background variation across recovered acquisition groups

1 October 2026 UTC. This is a post hoc description of the saved corner powers, using the existing source groups after the [metadata correction](BACKGROUND-SOURCE-GROUP-CORRECTION.md). It is not a new held-out study, calibrated noise model, confidence interval or independence test. Every particle and recovered group is retained.

| EMPIAR | Groups | Group mean patch energy, 5th–95th percentiles | Between-group fraction of particle-energy sum of squares | Test/pilot spectral-shape ratio range |
|---|---:|---:|---:|---:|
| 10028 | 229 | 0.35753–0.44785 | 0.159 | 0.932–1.067 |
| 10049 | 137 | 0.12022–0.13518 | 0.072 | 0.954–1.061 |
| 10076 | 351 | 0.05598–0.08123 | 0.117 | 0.829–1.081 |

Energy is the mean squared reference-adjusted non-DC DCT coefficient, pooled over four corners per particle. Group energy quantiles weight each group equally. The between-group fraction is the ordinary descriptive ANOVA decomposition of the observed particle energies: it includes finite-group sampling fluctuations and must not be called an intraclass correlation or a fraction of physical variance explained. All group sizes and means are retained.

For the last column, first divide the mean DCT power in the existing test split by the corresponding pilot-split power, then normalize the 63 ratios to mean one. The result measures spectral-shape differences after a diagonal pilot-based correction; it omits full covariance and does not validate whitening inside a molecule. The labels were fixed earlier at source-group level, but these images have been reused during development, so this is not a prospective validation experiment.

Equal-group versus particle-weighted pooled spectral shapes are also saved. These estimands differ when acquisitions contain different particle counts. No acquisition is treated as independently sampled merely because it has a distinct label. The weighted group means reproduce the pooled means to floating-point precision, and the within-plus-between sum-of-squares identity is checked independently.

This check replaces an incorrect claim of missing group metadata with an observable acquisition heterogeneity diagnostic. It does not supply a uniform noise uncertainty bound, remove signal contamination from corners, or resolve experimental calibration. The numerical gate continues to use its declared unit-noise simulator; these experimental measurements are not silently substituted into that model.

Source: `scripts/describe_acquisition_background.py`. Saved summary and all group moments: `results/uncertainty/development/acquisition-background-description-v1/`. No original spectrum arrays or splits were changed.
