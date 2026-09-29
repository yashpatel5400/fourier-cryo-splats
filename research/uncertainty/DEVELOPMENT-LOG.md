# Development log

## 2026-09-29: physical pilot and a failed coarse remainder bound

The first experiment used 96 fixed simulated poses, 40 complex Fourier samples
per particle, 40 real Gaussian coefficients, a known coefficient-error ball of
radius 0.25, and known Gaussian noise. This is an exact-dictionary feasibility
check, with a CTF-like radial transfer, not a biological reconstruction or a
practically calibrated uncertainty procedure. Its known bounds are simulation
assumptions. Exact scalar Gaussian coverage is available analytically; the rows
are not Monte Carlo coverage estimates.

The `physical-pilot.json` output is retained. For the uniform-view case, the
coarse remainder certificate at 0.5 degrees and 0.05 pixels had a half-width
0.737 times the no-data half-width. At 2 degrees it reached 0.980. Under the
preferred-view geometry it reverted essentially to the no-data interval. This
is a failure of usefulness, despite valid conditional coverage.

The first-order interaction between density and pose was subsequently retained
explicitly. The `structured-pilot.json` experiment reduced these width ratios
to 0.497 and 0.803 respectively. The missing-view feature remained largely
unidentified. These improvements justify a larger investigation, not a final
method or novelty claim. Runtime was below two seconds per small case locally.

The stress direction in this pilot was constructed to maximize the ridge
estimator's first-order bias. It is appropriate for a worst-case diagnostic but
must not be presented as representative average-case performance. A complete
comparison must also include prior-predictive calibration, many random signals,
multiple regularization choices, unregularized/sampling baselines, and each
method's own adversarial direction. A Bayesian credible interval is not designed
to promise uniform fixed-parameter coverage over a Euclidean ball.

## Structured local model

Write A_i(u)=A_i(0)+sum_a u_a D_ia+E_i(u). Keeping D_ia delta produces

    residual = A delta + J u + D(delta,u) + e + noise.

For a fixed weight block w_i, define T_i(w_i) with columns D_ia' w_i.
The density-pose interaction is bounded by

    |w_i' D(delta,u_i)| <= B ||T_i(w_i)||_op
                              <= B ||T_i(w_i)||_F,

when ||delta||<=B and ||u_i||<=1. The implemented Frobenius relaxation is a
convex group norm. It relaxes the shared rank-one delta*u_i' structure; it is
not claimed to be the exact nonlinear minimax solution. The remaining Taylor
bound is quadratic in the pose radius and includes the curvature of both the
pilot and the density perturbation. Full analytic Gaussian derivatives avoid
the old hard-cutoff discontinuity.

The multi-group solver has been compared with an independent CLARABEL conic
solution in both its small-Woodbury and direct-pixel-block regimes. Random
boundary checks verify the analytic interaction/remainder bounds. These are
numerical checks supporting the derivation, not a substitute for proof.

## Literature search status

Four broad Europe PMC queries retrieved 16,969 deduplicated records. Most are
biological applications, not relevant reconstruction methodology. API counts
and returned counts differ slightly for three queries (2, 4, and 1 records),
so the ledger explicitly marks those queries incomplete; it is not an exhaustive
systematic review. Title screening, primary full-text reading, reference
snowballing, arXiv/conference searches and version reconciliation remain ongoing.
The initial curated reading list has 60 works. Downloaded is distinct from read.

CAHRA's September 2026 preprint currently describes dataset construction; its
official website says the full competition analysis will be added in an update.
Do not describe that preprint as a completed benchmark-results paper.

## Independent-map representation audit

Downloaded EMD-2660, EMD-6487 and EMD-8434, with hashes and headers recorded.
They define known semisynthetic generators; they are not asserted to be true
densities for the corresponding real particle stacks. Reference signals use a
low-pass-resampled 64-cubed voxel map and an independently implemented NUFFT
forward projection, checked against direct summation and Cartesian FFTs.

At frequency radius six bins, 256 images per development split, and a finest
Gaussian lattice with 2,110 coefficients, the noise-free held-out image relative
RMSEs were 1.77%, 2.65%, and 7.73%. Corresponding mean FSC over shells 1–6 was
0.9989, 0.9968, and 0.9707. Coarser dictionaries had substantially larger error.
Thus representation error is observable even with zero measurement noise and
must be included in a density uncertainty analysis. These are low-frequency
development measurements, not new high-resolution reconstruction claims.

The audit also fits actual experimental pilot particles. All outputs are under
`results/uncertainty/development/representation`. Architecture diagnostics used
the split named `test` in this **development** data pool; that pool must not
subsequently be advertised as an untouched confirmatory test. Before final
evaluation, freeze the method and use fresh simulated draws and/or additional
source groups not used in this development pool. The older v0.1.0 study already
used these images as well. Keep this provenance explicit.

## Density-energy coordinates and memory reduction

An analytic whole-space Gaussian Gram matrix now permits sensitivity bounds in
integrated real-space density L2 energy, instead of a dictionary-dependent
coefficient norm. Numerically excluded null modes must be reported. A second
implementation change replaces wide density-pose derivative groups by their
pixel-space Gram factors. This preserves the group support function up to a
small conservative positive diagonal padding, reducing stored dimensions from
pixels times coefficients times pose components to pixels squared. Both changes
have independent numerical checks; larger-scale confidence experiments remain
to be run.

## Repeated fixed-pilot coverage development benchmark

`coverage-benchmark.json` contains 48 settings (12 seeds, two acquisition
geometries, two pose radii), each with 20 random bounded signals and a separate
adversarial case for each estimator. There are 2,000 independent measurement-noise
draws per specified signal, paired across methods. Exact Gaussian noise-marginal
coverage is reported alongside simulation. Bounds-violated rows are separate.

The structured-interaction certificate's minimum analytic coverage over these
in-class tests was 0.9912 at nominal 0.95. Dropping the quadratic remainder
produced a minimum of 0.9323. The coarse certificate was almost always maximally
conservative. Median widths pooled over the uniform/preferred-view designs were
0.962 and 1.000 times the no-data interval for the structured and coarse versions.
This establishes a remaining usefulness problem, especially under preferred
views; it does not establish practical superiority. Pooling geometries obscures
important differences, so final reporting must stratify them.

Narrow Gaussian priors failed some fixed-parameter boundary tests, while every
full Gaussian posterior passed its exactly matched linear prior-predictive
coverage identity at 0.95. These facts are compatible, not contradictory.
Widening the prior and unregularized inference gave much larger intervals.
Measurement-only sampling intervals missed coherent pose errors. Those controls
are necessary to interpret a favorable robust-coverage result fairly.

An initial incomplete run is preserved as
`coverage-benchmark-incomplete-covariance-factor.json`. It was interrupted to
factor the Monte Carlo correlation matrix rather than a poorly scaled covariance
matrix containing very large unregularized variances. The completed run uses
the stable factorization. Only the completed run is eligible for summaries.
