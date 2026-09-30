# Fixed-design projected noise calibration with simultaneous selection

30 September 2026. Declared after centered-noise-calibration-v1 and before any
outcome from the projection family below. Calibration and inference pixels
have already been inspected. This is a development study, not fresh validation.

## Conditional finite-family guarantee

Condition on an independent training/design sigma-field. Calibration rows are
independent Gaussian vectors Y_i~N(mu_i,Sigma), with arbitrary fixed means
and common covariance. Let Q_j have orthonormal rows, with d_j>=1, and let V
contain the inference weights pulled back to raw Fourier coordinates and
inflated by sqrt(exposure-group size). All Q_j and V must be independent of
the calibration/inference noises. Then Q_j Y V^T has d_j independent Gaussian
rows with common covariance V Sigma V^T. The existing noncentral trace bound
U_j, calculated using d_j rows and failure probability beta/J, satisfies

    P[tr(V Sigma V^T) <= min_j U_j] >= 1-beta.

Proof: orthonormality gives independent rows for each fixed contrast, and the
union bound covers all J contrasts regardless of dependence between them.
Thus selecting the smallest U_j is permitted with this explicit multiplicity
allocation. Pair it with the unchanged density/pose bias bound, spectral
allowance and Gaussian critical value. No correct low-rank signal model is
required: projecting an arbitrary mean leaves an arbitrary mean. This is a
classical Gaussian contrast/concentration argument and union bound, not a new
general inference theorem. Common covariance and independent design remain
essential. Noise-fitted projections are not permitted by this proof.

## Seven contrasts, fixed before outcome calculation

For each of the three existing 128-exposure calibration cohorts, use:

1. Identity, rank removed zero.
2. Constant-mean complement, rank removed one.
3. Constant plus three leading centered CTF-design left singular vectors.
4. Constant plus seven leading centered CTF-design left singular vectors.
5. Constant plus three leading centered pilot-prediction left singular vectors.
6. Constant plus seven leading centered pilot-prediction left singular vectors.
7. Constant plus fifteen leading centered pilot-prediction left singular vectors.

The CTF design is the 128-by-220 known real transfer matrix on radius-12
detector coordinates, with no orientation or image values. The pilot design
is the realified 128-by-440 forward prediction of the original independently
trained, unit-L2, 24-cube pilot, using the downloaded calibration metadata's
poses, translations and CTFs. Undo the same centering phase/data sign to
return predictions to raw image coordinates. No reference map, calibration
pixel, observed interval or reconstruction outcome enters either design.
Supplied metadata were estimated upstream; their independence from noise is
an unverified condition, particularly for the pilot-informed contrasts.

Use SVD relative threshold 1e-10. If fewer than the requested singular vectors
are numerically available, clamp to that rank and retain the duplicate or
reduced-rank case. Use deterministic orthonormal complements and report all
singular values, ranks, Q Q^T errors and removed-design residuals. J remains
seven even if contrasts coincide. Do not tune rank from a noise outcome.

## Fixed estimators and reporting

Use all twelve fixed estimators with their four existing pose-class bounds,
and all three completed 10049 cubic estimators. Preserve their weights,
centers, pilot amplitudes, bias bounds, reference transforms and no-data
fallbacks. First replay each uncentered SD at its original beta to 1e-9 and
the independent pilot scores/hash. Then use beta=.005/12 shared over seven
contrasts, Gaussian alpha=(.045-1e-6)/12, and delta=1e-6/12. Each pose-class
twelve-feature family is a separate sensitivity; the cubic three-procedure
family remains separate. No joint guarantee over all historical development
is asserted. Retain all per-contrast values, selected name, observed intervals,
zero exclusions, and original/registered reference agreement, including every
disagreement. Never interpret a reference class map as population truth.

The whole diagnostic has a 15-minute budget. Archive exact code, all hashes,
predictions, contrasts, projected calibration arrays and weights. Preserve
failures and every case. Before empirical application, verify projection
identities and the finite-family procedure against independently simulated
Gaussian observations with unequal means and correlated coordinates. Any
improvement remains conditional; no pixel download or estimator refit occurs.
