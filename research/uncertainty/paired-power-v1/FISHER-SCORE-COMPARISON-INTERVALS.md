# Reporting addendum: conservative comparison intervals

Declared while the first Fisher study is in calibration, after seeing training
means/criteria but before any new held-out rejection outcome exists. This does
not change a frozen score, threshold or probability bound.

Compare Fisher minus matched rejection projections within each stack, feature
family, calibration method, density cap, sample size, deletion and amplitude.
There are 9,450 pairs. Conditional on frozen training and calibration, obtain
97.5% two-sided Clopper-Pearson intervals for the two single-image event
probabilities. Transform each by its own binomial rejection function. Subtract
opposite endpoints to bound their difference. A union bound gives pointwise
at least 95% coverage for this difference under the declared test simulator,
regardless of dependence between the two scores. This is deliberately
conservative; it is not a bootstrap, a simultaneous comparison across cells,
or an assessment of training/calibration randomness. Retain the 2x2 joint
single-image event counts so later dependence-aware diagnostics are possible.

For both variants report the regularized training Rayleigh criterion using the
same covariance and the observed separation/variance in the held-out score
arrays. Training optimality is not a claim of test-power optimality. No
multiplicity-unadjusted selection of the better variant becomes a final test.
