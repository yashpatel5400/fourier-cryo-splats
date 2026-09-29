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

## Real geometry and physical targets

A new development benchmark uses the three datasets' inference-half geometries
and CTFs, normalized real-particle pilot reconstructions, 128/512 images,
local Gaussian density averages, and fixed-pose/0.5-degree cases. The density
ball is expressed in whole-space L2 energy and selected as 0.25 of normalized
pilot energy before simulation truth. It remains an assumed finite-dictionary
class. Widths are generally large at averaging width 0.03 of the image field;
this scale includes appreciable target energy beyond the observed radius-five
frequency band. The implementation now exposes this physical scale explicitly.

The matrix-vector PDHG alternative passes an independent conic check across
8 orders of target scale. On the initial real-geometry case it was slower than
majorization (1.81 versus 0.34 seconds at comparable gap tolerance), so this is
not a speedup claim. Both algorithms remain available for scaling comparisons.

## Ambient-space bias audit and a method revision

The representation experiment motivated a stronger revision: compute the
confidence bias on an explicitly larger voxel space, even when a Gaussian
subspace produced the image weights. The identity separating the projected and
orthogonal adjoint residuals quantifies the uncertainty omitted by a fitting
basis. Target-directed subspace enrichment adds those omitted directions using
only design/target information, with ambient primal and restricted dual bounds.
A separate full-space matrix-free solver is implemented as a comparator.

The first EMPIAR-10028-geometry probe (fixed poses, 24^3 grid, radius-0.35 support,
128 images, prescribed SNR 0.1, unit-normalized pilot/reference densities and
radius B=2) showed severe undercoverage for dictionary-only intervals on a
known EMDB-map generator. Auditing restores the assumed-model coverage, and
enrichment reduces widths, but its gap was still large after 60 rounds in that
probe. The matrix-free comparator converged more directly. These are important
negative controls on both inference validity and optimization efficiency.
The all-three-dataset audit is being preserved separately; none of this establishes
unknown-pose or continuous-space coverage. Unit-normalized simulated densities
satisfy B=2 by the triangle inequality; this does not estimate experimental B.

The completed full-space fixed-pose batch reached a 0.5% sum-objective gap in
all 12 targets (three geometries, two locations, two smoothing scales), with
2.4–11.1 seconds of matrix-free solve time per target in this run. Median width
fractions for the broader scale were 0.050, 0.073 and 0.072 across datasets;
the finer scale gave 0.687, 0.716 and 0.719. Greedy enrichment converged in 8/12
within 100 rounds and was not consistently more efficient. Both methods and
all negative dictionary-only coverage controls are retained. This exposes a
scale-dependent information limit under the chosen class, not an atomic-resolution
claim. `paper/main.tex` has been rewritten to present the uncertainty problem,
proofs and these development results; it explicitly states unfinished work.

## Joint ambient-space and nonlinear pose experiment

The ambient voxel operator now supplies analytic phase derivatives, exact
pixel-space Gram factors for the full-density/pose interaction, and uniform
second-order remainder bounds over the declared support. The existing
majorization solver operates through forward/adjoint NUFFTs, using the exact
Fourier-column Gram diagonal as a preconditioner. This combines the two earlier
separate guarantees in one implementation. Independent finite differences,
nonlinear perturbations, direct sums and operator/dense agreement pass; the
suite now has 31 tests.

`ambient-pose.json` contains nine development cases: all three acquisition
geometries, 128 particles, the broad axial contrast, and 0.1/0.5/2-degree radii.
At 0.5 degrees the full certificate width is 0.200/0.239/0.218 times the no-data
width. At two degrees it rises to 0.866/0.996/0.910, showing substantial remaining
conservatism/information loss. Fixed-pose intervals reach analytic coverage as
low as 0.904 on tested allowed perturbations. The full certificate's minimum
is 0.999705; this conservatism is not described as optimal calibration. The
large-radius optimization gaps did not converge in 100 majorization steps;
a separately named longer refinement run is underway. Model validity does not
depend on optimization convergence, but efficiency claims do.

Stock neural cryoDRGN has also actually executed a one-epoch CPU timing run.
It is only a compute profile, not a fitted/converged neural baseline. The original
backprojection baseline remains distinct. Its timings and exact command are
in COMPUTE.md and `neural-runtime-probe/profile.json`.

## Combined audit, power and real neural baselines

The full supported-voxel and nonlinear-pose certificate is now implemented and
checked. Nine geometry/radius cases and a longer 2-degree optimization run are
preserved. One longer run (10049) misses the 0.5% gap threshold at 500 iterations;
its valid feasible interval is not called an optimum. At 0.5 degrees the broad
contrast widths are approximately 20–24% of the no-data width. This alone is not
evidence of a biologically useful inference: a subsequent sign-power experiment
finds essentially no power for the three deposited-map axial contrasts under
the prescribed B=2 class. Artificial target-directed signals are easier.
A paired central-average experiment is running to distinguish target size from
interval validity. These remain development outcomes, not a frozen final test.

Stock cryoDRGN homogeneous neural training (3 layers, width 256) has completed
20 epochs on two audited exposure-group halves for each of the three stacks.
Held-out evaluation uses identical Fourier samples for neural, Gaussian and
voxel predictions. A map-evaluation bug was caught before publication: the
image FFT helper only transforms two dimensions. An explicit volume FFT helper
and an axial-frequency regression test now prevent that mistake. Neural map
FFTs are additionally checked against direct network values; provisional FSC
files from the runner are replaced by the corrected evaluator. Neural error
still improves at 20 epochs on 10028, motivating a recorded 20-to-60 epoch
continuation. This is actual neural fixed-pose training, not ab initio inference.

A second-order pose expansion now retains pilot curvature and density/curvature
interaction in Euclidean support functions, with a uniform cubic remainder.
Mixed derivative and independent nonlinear tests pass. In the first EMPIAR-10028
contrast probe, width/no-data decreases from 0.200 to 0.143 at 0.5 degrees and
from 0.865 to 0.478 at 2 degrees. Tests remain conservative. The extension is
being checked on all three geometries and both targets; this probe alone is
not a finalized method comparison. Full-space Gaussian posterior controls also
pass their separate, exactly matched additive-linear prior-predictive identities.

An independent scalar-noise upper-bound module implements a standard noncentral
chi-square argument. Its assumptions allow fixed signal in calibration samples
but require independent Gaussian coordinates with a common variance. Synthetic
correlation violates the scale guarantee (a failure rate around 0.50 instead of
0.01 in the rho=0.5 control); it is not silently treated as valid whitening.
These are calibration experiments, not evidence that raw experimental noise
satisfies the model. The current complete unit suite has 35 passing tests.

The curvature comparison is complete for both targets, all three geometries and
three radii. At 0.5 degrees, central-reference sign power increases to about
0.999 or higher for all three geometries; small axial contrasts remain mostly
undetectable. At 2 degrees, even central averages remain difficult. The whole
regularization sweep for group bootstrap uses 999 resamples and 1,000 independent
noise replications per setting; counts are shared across noise replicates and
this conditioning is recorded. Weakly regularized full-space bootstrap often
has near-nominal reference coverage with smaller widths than the robust audit.
Strong regularization and restricted dictionaries can fail. Report both outcomes.

A physical-unit fixed-pose grid audit uses a common 64-cubed inference center and
truth generator. The 32-cubed weights change in width by less than about 0.6%
on that common fine grid for this band; 16-cubed results can differ substantially.
The new Fourier-moment adjoint computes first/second pose audits on finer grids
without full derivative Grams. The first 64-cubed probe gives the same width
ratio (0.97936768556) with NUFFT and direct sums. Profiling under concurrent load
was 289.34 versus 17.32 seconds, respectively; this is not a controlled hardware
speed benchmark. The implementations use different native threading behavior.
Larger-grid checks for all datasets/targets/radii and a higher-bandwidth finer-
target probe are in progress. Current source tests cover both adjoint backends.

## Shared-density audit, nonlinear stress search and extended comparisons

All 18 curvature fits (three geometries, two targets, three pose radii) now
meet the 0.5% optimization-gap tolerance. The subsequent Fourier-moment pose
refinement audit is complete on grids 24/32/64 for both targets and 0.5/2 degrees.
A higher-bandwidth 32-grid, radius-12, width-0.03 probe gives a no-data first-order
interval but quadratic relative width 0.44855 at 0.5 degrees for EMPIAR-10028.
Its expanded comparison is still running; the probe is not three-stack evidence.

A projected nonlinear search now eliminates the common density ball analytically
and searches 128 five-dimensional pose balls. Four starts per bias sign, 150
Adam steps each, give feasible lower bounds on worst bias. All 12 geometry/target/
radius settings are complete. The Mac MPS backend uses float32 during search;
saved candidates are reprojected and reevaluated by independent float64 Fourier
sums. Lower bounds do not imply globally worst cases. These searches expose
considerable slack in the separate density-interaction triangle bounds.

A deterministic spectral post-audit preserves the common density direction in
linear/quadratic pose terms and the nominal density residual. It uses classical
weighted block operator-norm inequalities, comparing joint, affine, separate and
triangle relaxations. Across the 12 settings, interval widths fall to roughly
0.50–0.65 of their original values without changing image weights or assumptions.
All completed feasible nonlinear adversaries remain below these tighter upper
bounds. This does not certify optimality for the new objective. A 64-grid probe
agrees closely in physical units; the full larger-grid post-audit is running.

Neural fits reached 60 epochs for 10028/10049 and 100 for 10076. Tuning-only
checkpoint selection gives epochs 60/50/100 with development test NMSE
0.829325/0.930173/0.870705. All conditional half-FSC values remain above 0.143
through the sampled limit; no threshold crossing or gold-standard resolution
is inferred. The earlier checkpoints and resume records remain available.

The complete test suite has 38 passing tests. The development manuscript now
contains the bootstrap penalty sweep, curvature/power comparison, longer neural
runs, shared-density derivation, and a more detailed target-specific literature
appendix. Selected citation metadata are recorded from Crossref/Europe PMC;
Crossref rate-limit failures are preserved rather than treated as successful
retrievals. No scientific Fable review has yet been requested.

The full shared-density 64-grid post-audit is now complete for all 12 settings.
It preserves the reduction, with finer-grid widths about 0.50–0.66 of the same
fine-grid triangle-bound interval. All 96 saved nonlinear adversaries from the
coarse class remain below their corresponding tightened coarse upper bounds.
The assumption-violation diagnostic suite is being run separately; it is not
included as completed evidence in this checkpoint.
