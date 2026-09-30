# Centered Gaussian calibration: theory and development protocol

30 September 2026. Declared after all three fixed cubic experimental intervals
failed to exclude zero. Their SD upper bounds include signal energy and
conservative exposure dependence; they are not measurements of pure noise SD.

## Elementary orthogonal-contrast extension

Let Y have n independent rows Y_i ~ N(mu_i, Sigma), with arbitrary deterministic
means and one common covariance. Fix an (n-1)-by-n matrix Q with QQ^T=I and
Q1=0. Then the rows of QY are independent Gaussian vectors with common
covariance Sigma and means (Q mu)_i. Independence follows because row j,k
cross-covariance is (QQ^T)_jk Sigma. Thus the existing noncentral trace bound
applies with n-1 observations. No equality of means is assumed. The sum of
squared contrast norms equals sum_i ||Y_i - mean(Y)||^2. Orthogonal contrasts,
Gaussian independence and the resulting degrees-of-freedom adjustment are
classical, not a new statistical principle.

For any fixed weight matrix V independent of these noises, apply this result
to QYV^T. For arbitrary jointly Gaussian dependence within each inference
exposure, use the existing sqrt(group size) row inflation of V. This changes
neither the inference estimator nor its pose/density bias. A data-fitted
projection Q other than this fixed constant-vector contrast would need a
separate independence or regression argument. No such projection is used.

The centered upper bound is not automatically smaller: its statistic decreases
but its effective sample size is n-1. If E is the uncentered projected energy
and D=n||mean(YV^T)||^2, the exact ratio is

    U_center / U_raw = (1-D/E) n t_n / ((n-1) t_(n-1)),

where m(t_m-1-log(t_m))/2=log(1/beta). Retain cases where the new bound is worse.
When E=0, report both zero bounds without defining the ratio. Do not take an
unadjusted minimum of centered and uncentered confidence bounds.

## Experiment declared before computing centered outcomes

Reuse all three noise-calibration-v1 cohorts and all twelve fixed estimators,
with every existing fixed/shift/one-/two-degree bias bound. Also retain all
three completed cubic designs on 10049 region 1. No new target, weights,
reference registration, amplitude fit, calibration samples or particle pixels
are selected. Previous data have been inspected; this is development reuse.
First replay every uncentered SD upper against the published fresh/calibration
application to 1e-9 relative/absolute tolerance. Verify all input hashes and
unchanged geometry. Use the deterministic Helmert orthogonal contrast of the
128 calibration rows, giving 127 independent rows under the stated model.

Use exactly the original per-feature allowances: Gaussian tail
(.045-1e-6)/12, calibration .005/12, spectral 1e-6/12. Each pose class has its
own twelve-feature family; the separate three-design cubic family uses the
same per-design allowance. These are alternative procedures, not a selection
rule over the historical record. Recompute only noise SD and half-width/fallback,
preserving raw observed centers, existing bias bounds and pilot-only fallbacks.
Report all 48 fixed-weight and three cubic intervals, original/registered
reference agreement, zero exclusions, raw/centered energy ratio and SD ratio.
No experimental density-coverage claim follows from these conditional results.

Tests check orthonormality, removal of a common shift, exact projected-energy
identity, and a separate known Gaussian-covariance simulation with unequal
means. The diagnostic has a 15-minute budget. Archive exact code, all input
hashes, contrast matrices and centered calibration arrays; preserve failures.
