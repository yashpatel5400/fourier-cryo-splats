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

## Assumption failures and a continuous-density extension

The six completed assumption-stress runs (three geometries, both targets) each
retain 273 records. They include noise scale, within-particle correlation,
artificial estimator-aligned global noise, elliptical Student-t noise, coherent
defocus errors, support-only violations that retain the total L2 radius,
equal-population two-state heterogeneity, and pose-radius exceedance. The
fixed-pose estimator's own nominal bias boundary attains 0.95 analytically,
checking the control. Some support and view-coupled heterogeneity cases drive
the tighter nonlinear interval's coverage near zero. Correlated noise and scale
errors can also fail. The tested low-bandwidth defocus offsets and Student-t
cases do not produce undercoverage for the tighter nonlinear intervals. Do not
suppress those negative findings or portray constructed adversaries as typical.

An analytic continuous fixed-pose audit now treats the unknown density as any
L2 function on the unit cube. This is broader support than the prior spherical
voxel class. Exact sinc Grams and truncated-Gaussian target integrals avoid an
unjustified continuum conclusion from grid agreement. Independent quadrature
and conic optimization checks pass. Constant-cell pilots/reference generators
use exact cell Fourier transforms and target integrals, including the cell
volume, sinc and half-cell phase. They are not point-sampled volumes.

The full three-geometry/four-target dense continuous optimization is complete.
Every case meets its continuous 0.5% sum-objective gap. The first broad central
10028 probe gives width/no-data 0.08733057. Transferring the original sphere-only
weights without refitting can be much worse because the support class changed;
that unfavorable comparison is preserved. For the cell-fitted broad central
probe, using only the 64-cell projection reduces width by merely 1.7%, yet
coverage at an explicitly continuous worst-bias direction is about 0.923.
Reference-map coverage remains one. This is a worst-case audit, not a claim of
typical reference-map undercoverage.

A matrix-free Gauss–Legendre/type-3-NUFFT implementation includes an analytic
quadrature remainder in both primal and continuous dual bounds. Requested
FINUFFT tolerance and roundoff remain numerical, not interval-arithmetic,
checks. Its first 128-particle broad-target result agrees with the dense
solution to displayed precision but is slower: about 313 seconds versus 12
seconds for dense optimization/audit, under concurrent load. Stored quadrature
arrays use about 2.2 MB versus 419 MB for the dense Gram (these are not peak RSS).
A standard 256-column pivoted-Cholesky preconditioner reduces this to about
141 seconds with 23 MB stored arrays, still slower than dense. A higher-rank
preconditioner is being profiled; no speed superiority is claimed. The current
complete suite has 41 passing tests. The continuous nonlinear-pose extension in
CONTINUOUS-POSE-NEXT.md remains a hypothesis, not an implemented result.

## Completed continuous pose audit and first frozen prediction cohort

All twelve continuous-pose post-audits are complete: three real acquisition
geometries, two broad targets, and 0.1/0.5-degree joint pose budgets. Relative
half-widths range from 0.1102 to 0.4435, and all are narrower than the no-data
option. Exact reference-generator coverage is essentially one. At a coherent
nonzero pose, eliminating the full continuous density ball gives feasible bias
between 0.24 and 0.86 of the uniform upper bound. These are lower bounds, not
global pose optima; considerable slack remains. The probe encountered a JSON
serialization error for a NumPy boolean at its second setting; the full sweep
reran with explicit Python scalar conversion. Both the failed probe's first
result and exact source snapshots are retained. Each complete audit took roughly
38–40 seconds and stored 725 MB of main field/Gram arrays, not peak RSS.

The rank-1024 matrix-free fixed-pose implementation completed all twelve targets
and agrees with dense continuous optimization. It also completed broad targets
with 1024 particles on all three geometries. A higher-bandwidth continuous
radius-12 probe is still running; its first central target reaches relative
width 0.13024 but takes about 864 seconds. Scaling is not uniformly cheap.
The full unit suite now has 45 passing tests, including continuous nonlinear
fields, polynomial quadrature remainders and exact density-ball stress checks.

Commit 72a6dc0 freezes an additional-exposure experimental prediction protocol
before selected image downloads. It selects 4096 particles per stack from
115/71/185 source groups absent from the entire development pool and locks the
existing Gaussian, voxel and neural checkpoints/scales by hash. No new density
coverage or independent-pose claim follows from this conditional prediction
study. Downloading and evaluation are in progress. Uncertainty confirmation on
fresh simulated designs remains separate and is not yet frozen.

The higher-band finite-voxel comparison is now complete for both targets and all
three geometries. Every quadratic solve meets its gap tolerance. Relative widths
are 0.166–0.231 at 0.1 degrees and 0.442–0.477 at 0.5 degrees. The continuous
higher-band probe on 10028 completes both fine targets with relative widths
0.13024/0.12959 and runtimes 864/883 seconds; this is only one geometry so far.

All frozen additional-exposure prediction runs are complete. Gaussian/voxel/
neural NMSEs are 0.835319/0.841379/0.829369 (10028),
0.931579/0.932359/0.930739 (10049), and 0.874799/0.879078/0.868657 (10076).
All six paired whole-exposure bootstrap contrasts favor neural prediction even
with the prespecified six-contrast adjustment. No models were retuned. The raw
transfer totals for this new cohort are approximately 2.37/0.60/2.56 GB, with
verified range hashes. This supplies no experimental-density coverage label.

A separate continuous-v1 protocol is being committed before its simulation
outcomes. It specifies 12 new acquisition subsets, four targets and two pose
budgets (96 audits), with 2304 signal/scenario records and exact Gaussian coverage
cross-checked by 10000 noise draws each. A smoke run used only eight old
particles and unrelated seeds. The forthcoming simulation study includes all
selected procedural generators, continuous boundaries and assumption failures.


## Integrated moments, wider poses and exact nonlinear stress checks

The complete 1/2/5-degree broad-target development sweep has 18 original
no-data fallbacks. Integrating Fourier-column derivative envelopes over the
cube and retaining actual estimator weights reduces cubic remainder bounds by
6.79--14.74 times across 30 small/large-angle settings. The refined audit has
10 no-data fallbacks; all six one-degree relative widths lie in 0.300--0.532.
The proof uses classical moment comparison, Holder and Taylor bounds, not a new
general concentration theorem. Three independent numerical tests bring the full
suite to 48 passing tests. No experimental pose calibration follows.

The second frozen study was committed as 8d66785 before its outcomes, while v1
continued unchanged. It uses disjoint particle subsets and fresh simulation
seeds, broad targets and one/two-degree budgets, with 48 audit settings and six
matched procedures. Its protocol explicitly discloses that the refinement was
motivated by development outcomes and that partial v1 outcomes were already
known. Both studies remain running at this entry.

Continuous nonlinear adversaries now cover all 18 target/geometry/angle settings
at 0.5/1/2 degrees, three starts and both bias signs (108 candidates). The search
uses an order-16 quadrature surrogate and MPS float32, with final feasible poses
projected and evaluated in float64 using independent sinc-integral formulas.
All candidates stay below the integrated-moment upper bound. Best feasible
biases improve coherent-pose controls by 1.06--2.54 times but attain only
0.16--0.57 of the upper bounds. Direct-field discrepancies are at most 5.35e-6.
This does not establish global sharpness or validated arithmetic. Total recorded
setting time was approximately 606 seconds under concurrent load.

The survey now has 92 curated candidates, with new targeted readings on pose
quality/transfer, image-level calibration, empirical-Bayes resampling and
heterogeneity benchmarks. Access failures and version distinctions remain in
the ledger. A stock cryoDRGN ab initio 32-particle timing probe completes locally
in 49.4 seconds; it is explicitly not a converged reconstruction baseline.
No cloud compute was rented. The provisional first GPU batch cap remains $100,
subject to profiling and an actual offer including storage/transfer charges.


## Higher-band completion, background diagnostics and fixed-length efficiency

All six radius-12 continuous fixed-pose targets on the three stacks now meet
0.5% sum-objective gaps. Width/no-data ranges from 0.120 to 0.238, but reference
sign power is below 1.1e-9 throughout. Solve/audit times range from 483 to 1223
seconds under concurrent load. Cell64 projections of the same residuals give
continuous-boundary coverage only 0.791--0.894. This extends the bandwidth check
without claiming high-resolution biological utility or isolated speed rankings.

Post-hoc annular-background diagnostics use all 12288 already-evaluated fresh
images. Mean nearest-neighbor pair products are 0.197--0.292, compared with near
zero in identically centered/scaled white-Gaussian controls. Molecular signal,
ice and preprocessing can contaminate the background: this diagnoses an
experimental assumption gap, not a calibrated noise covariance. All axis/lag
summaries and whole-exposure bootstrap aggregates are retained.

A classical two-point Gaussian-testing bound now evaluates the necessary
half-length of every uniformly valid deterministic-length interval, even with
nonlinear centers. Across all 18 completed fixed-pose continuous settings, its
numerical lower/upper ratios are 0.840--0.882. Thus the reported widths lie within
1.20 times the optimal fixed length in the stated density/noise class, subject
to numerical integration accuracy. This is standard lower-bound theory, not a
new general theorem, and does not establish sharpness with uncertain poses or
for variable-length intervals. The proof and post-hoc design-only computations
are independent of frozen-study coverage outcomes.

The main paper now emphasizes the continuous method; finite-representation
methods and early experiments remain in the appendix. The literature ledger
adds Gold-standard local validation (93 candidates total) and targeted EMMIVox
reading. No scientific reviewer has yet been invoked.


## Completed frozen studies and first review candidate

Both locked continuous studies are complete. V1 contains 96 settings, 2304
signal/scenario records and 48 converged fits; minimum coverage across its 1536
in-class records is 0.993772 for the pose audit, 0.939215 for fixed poses and
zero for the 24-cell audit. The severe projection-audit failures occur at
constructed continuous boundaries. V1 has no no-data fallbacks, although fine
feature widths nearly reach the no-data bound. Deliberate violations did not
expose failure of the pose audit at their selected severities.

V2 contains 48 settings, 1152 records and 24 converged fits. The integrated
pose audit covers all 768 in-class records with analytic values rounded to one
in float64. Fixed-pose and cell-projected minima are 0.671143 and 0.125960.
Its 16 no-data fallbacks improve on the spatial-maximum correction's 48.
One-degree relative widths are 0.303676--0.500927; at two degrees only 10076
remains narrower than no data (0.865642--0.986669). In-class cell-generator
sign power is at most 2.2843e-5: useful structural detection is not established.

All raw JSON records, fitting gaps, widths, sign-power values, scalar Monte
Carlo checks and stratum summaries are retained. Reporting scripts and plots
were developed after outcomes; experimental algorithms and locked protocols
were not altered. The manuscript now integrates all frozen outcomes, including
unfavorable utility results, and is a candidate for independent review.
The classical fixed-length efficiency proof also makes explicit its generic
1.19756 factor from Gaussian testing and a 0.5% primal-dual gap; this is not
a novel minimax theorem. No scientific Fable review has occurred at this entry.


## Authentic Fable 5.1 review, round 1

The complete candidate at a4117ab was sent through the authenticated Claude
CLI to the exact claude-fable-5-1 model. The provider response confirms that
canonical model. The review took approximately 770 seconds and is preserved
verbatim with the full input packet and source/PDF hashes. It says reject,
confidence 4/5, not a strong ICML contender. Its main objections are useful
precision, experimental calibration, pose-aware weight optimization, pose
scaling and incremental novelty. A finite substantive response plan is saved;
no claim of reviewer approval is made. The goal remains active.

While the review ran, a development probe was launched to select a frequency
band using the pose-aware width, holding physical noise and geometry fixed.
It uses the old 10028 development geometry and does not alter either frozen
study. This is an exploratory design check, not resolution of the review.

## Post-review numerical and experimental development

The bandwidth-only probe completes all twelve settings. At 0.5/1/2 degrees its
best relative widths are 0.209/0.402/0.945; reference sign power is negligible.
Changing bandwidth alone does not solve the usefulness problem.

New matrix-free polynomial pose fields use ten NUFFT moments without a global
dense Gram. Dense/adjoint/weight-gradient checks pass. A Gaussian-start power
upper bound includes an explicit numerical failure probability; a Ritz value
is only a lower diagnostic. Audits at 128 and 1,024 particles complete in
50.3/353.3 seconds with peak resident memory 390/434 MB. Their widths are
0.438/0.740 of no data. The radius-12, 10,000-distinct-particle matched-filter
audit remains in progress, with its exact source snapshot saved at launch.

A convex fixed-scale pose-aware objective now optimizes estimator weights.
An independently assembled tiny CLARABEL problem brackets its optimum with
relative gap 0.000260. The first attempt lost JSON serialization because two
comparison flags were NumPy booleans; the failure log is retained and the rerun
casts these to Python booleans. The 20-iteration larger probe has relative width
0.255 but a 97.8% optimality gap. Its six independent feasible nonlinear searches
reach 0.432 of the bias upper bound and minimum reference sign power is
2.95e-103. It is not a converged or useful-feature result. Longer optimization
and all three geometries/two targets/three angles at 0.5 Angstrom shifts are
running, with best-weight checkpoints and source snapshots.

The new critical-value helper retains both Gaussian tails. Its audit of 288
old widths changes them by at most 1.42e-13 and changes no fallback decisions.
Frozen dependencies remain byte-identical: 23 original and 26 follow-up source
hashes rechecked. The old shortcut remains only in the preserved frozen code.

An experimental conditional-noise attempt now runs on all three stacks. It
uses a common-covariance Gaussian trace envelope that permits nonzero signal,
with one particle per fresh calibration exposure. The grouped inference runs
use 58/30/93 particles and 57/35/92 calibration groups. All eighteen intervals
contain the approximate deposited-map value and none excludes zero. All 10
Angstrom targets switch to no data at nonzero pose budgets. Unknown pose-bound
calibration, physical density radius, covariance commonality and full-data
consensus-pose dependence are still unresolved. The first 128-particle attempt
is retained with its additional within-exposure independence limitation.

Uniform nonlinear CTF/gain/envelope sensitivities are implemented without
dividing by CTF zeros. On all six radius-12 fine-target estimators, a 100 Angstrom
defocus envelope expands relative widths to 0.647--1 and causes two no-data
fallbacks. All six fall back at 500 Angstrom and under combined sensitivities.
The triangle bound loses cancellations; these outcomes do not establish
impossibility for tighter methods.

The expanded suite passes 58 tests in 4.34 seconds. The manuscript now has
30 pages including references and appendices, with the main text ending on
page 7. The revised first page, final main-text page and new methods/results
pages were visually inspected; no overflow or undefined-reference warnings
remain. The survey adds targeted primary reading of Cai--Low and the
Kuczynski--Wozniakowski technical-report introduction (95 ledger candidates).
Low 1997 full text remains inaccessible from the attempted primary URL.
No second Fable review or positive verdict is claimed.

### Larger audit, conditional Fourier baseline and visual review input

The 10,000-distinct-particle radius-12 audit completed with 220 independent
Fourier pairs per particle, measured runtime 4,965.23 seconds and peak RSS
1,948,549,120 bytes. Relative interval half-width is 0.6551 at one degree and
0.5 Angstrom translation. Its randomized spectral upper is 1.1359 times the
Rayleigh lower diagnostic. Weights are deterministic matched-filter weights
for a broad target; this is a numerical scale result, not an optimized
reconstruction or useful fine-feature result. No GPU rental was needed.

The first completed pose-aware 100-iteration grid case (10049, center, 0.5
degrees, 0.5 Angstrom) has width 0.22086 of no data but a 25.27% optimization
gap. Six independent nonlinear searches reach 0.60809 of the bias upper;
minimum reference power is 8.44e-5. Remaining grid cases continue to run.

An independent implementation of the Ullrich fixed-pose Gaussian variational
objective completed on all three geometries, with 1,024 particles, radius 12,
10 Angstrom targets and three declared prior scales. Both likelihood and prior
terms are included through the exact Gaussian optimum. CTF and pilot-mean
extensions are stated explicitly. All 18 diagonal feature intervals are wider
than the corresponding full posterior intervals, by factors 1.073--1.415.
Full posterior prior-predictive coverage is exactly 0.95, whereas fixed-map
coverage often fails. Interpolation/representation forward errors of 16--25%
and shrinkage are material; these failures must not be attributed exclusively
to a diagonal covariance approximation. The new plots retain the full sweep.

An exact continuous support isometry and outside-density sensitivity bound are
implemented and independently tested. The initial 10049 half-side support case
has hard-support relative width 0.417 but cropped-reference power below 1e-150;
an outside radius of 0.1 makes it vacuous. The original reference's outside
norm is 0.320, so it is not an in-class example. The 0.75-side case is also
vacuous at all tested outside radii; the full-cube control is still running.

The actual requested Fable model correctly described the test figure. The
future review runner includes all rendered manuscript pages and explicitly
compacted iteration histories. A local-only 31-page packet preflight passed;
no round-2 scientific review has been invoked. The paper now has 31 pages,
with the main text still ending on page 7. Revised mathematical/results pages,
the new table and the six-panel baseline figure were visually inspected. The
expanded suite passes 62 tests in 5.06 seconds. The scientific readiness
objections remain open, and the authentic first review still recommends reject.

The same-estimator continuous audit of all 18 Fourier baselines also completed.
It charges the known pilot/interpolation center offset explicitly. Widths inflate
56--352 times, remain 0.257--0.815 of no data, and have numerically zero reference
sign power even after recentering. These are still fixed-pose, known-noise
conditional checks. The survey now has 96 entries after targeted reading of
sections 1--5.4 of the September version of an adjacent radiative-CT uncertainty
preprint; the earlier search title was stale. Its full-text download, version
and read scope are recorded without claiming a complete critical review.

CPU profiling preserved both outcomes. Reusing FINUFFT plans gives only 1.032x
median speedup. An isolated CPU-only OpenMP wheel gives 1.85x/2.21x/2.63x at
2/4/6 threads for the tested pose Gram matvec, with maximum relative differences
below 4.5e-15. The profiling process never imports Torch; the existing Mac
runtime remains unchanged. No running experiment has been silently switched
to the alternate numerical binary. Threaded execution is a possible next
optimization, not a completed accelerated end-to-end study.

### Checkpointed CPU continuation and a negative dual diagnostic

The original 10028 center/0.5-degree optimization failed with ARPACK reporting
seven of eight requested eigenvectors converged after 101 iterations. Its raw
failure and best-weight checkpoint are preserved. The other original grid
processes were then deliberately interrupted after validating and backing up
all saved checkpoints. They continue in isolated two-thread CPU processes,
with the original fixed-pose weights still defining block scales and entirely
fresh final-audit seeds. Additional iteration budgets are explicit. The exact
source/runtime hashes and process migration are recorded in CPU-CONTINUATION.json.
Completed 10049 center cases at 0.5 and one degree are reused, not rerun silently.

The optimizer now accepts available unit modes from a partial eigensolver
convergence only as approximate optimization directions. Its final upper bound
still uses separate randomized probes; Ritz values are never upper bounds.
A failure with no available mode remains a failure. A regression test forces
partial convergence and verifies the independent certificate path.

The newly completed original 10049 one-degree center case has relative width
0.33558, an 85.33% objective gap, and minimum reference sign power 7.07e-25.
Six nonlinear stress searches reach 0.46458 of its bias upper. It is not a
converged optimization or a useful detection. A separate feasible dual-mixture
postprocess reduces the gap only to 85.12%, without changing any estimator,
interval or upper bound. Thus an inadequate softmax mixture alone does not
explain the large gap. Twelve unit spectral directions, their feasible mixture
and the unchanged source fit are retained in the postprocess record.

### Directional noise calibration and the first completed continuation

All eighteen directional-noise re-audits are complete across the three stacks.
Weights selected using the first fresh-exposure half are held fixed; only the
second half calibrates their scalar noise variance. Raw-frame weights correctly
undo each supplied centering phase before projection. The direct trace lemma
allows correlated projected coordinates and arbitrary signal means. This does
not establish the remaining common-covariance, pose-independence or density
assumptions, and this exploratory diagnostic reuse is not a new frozen study.

Noise standard-deviation upper bounds are 0.13250--0.28481 times the original
ones. Four of eighteen cases still use the no-data fallback. Exactly one broad
10076 fixed-pose interval excludes zero; none does at one or two degrees. All
approximate deposited-map values remain inside the intervals, which is not a
coverage label. Full rows, unchanged-weight hashes, proofs, table and plots are
retained under directional-noise-audit/summary and DIRECTIONAL-NOISE-CALIBRATION.md.

The completed 10076 center/half-degree CPU continuation has relative width
0.16056 and a certified sum-objective gap 0.05023. Its optimization stopping
message alone is not the convergence certificate. Six independently rechecked
nonlinear candidates reach 0.64983 of the bias upper. Minimum analytic sign
power on the three reference-generator scenarios is 0.67032. This conditional
broad-feature result is more useful than previous cases, but does not establish
performance at one/two degrees or at a finer structural scale.

The complete expanded test suite passes 68 tests in 7.75 seconds. Tests include
exact correlated/noncentral calibration controls and the raw Fourier-coordinate
pullback identity; passing them cannot validate archive noise assumptions.

The revised PDF has 33 pages including references and appendices, with the
main text still ending on page 7. The directional calibration derivation, all
six plot panels and both experimental tables were rendered and visually
inspected. The build has no overflow or undefined-reference warnings. No new
Fable scientific review has yet occurred; the first rejection remains active.

The clustered-mode profile on a preserved 10028 checkpoint recovered all eight
modes with subspaces of size 33/49 (246/237 Gram calls); size 17 returned only
three converged modes after its bounded 15-restart run (121 calls). The partial
mode count changes the smoothed surrogate and can disrupt line search. The
10028 continuation was therefore checkpointed and restarted with size 49,
unchanged block scales and new final-audit seeds. Other running cases retain
their original source/runtime snapshots. RITZ-CONTINUATION.json records the
migration; no hidden replacement of prior output occurs.

The legacy smoothed 10028 one-degree probe (translation about 0.201 Angstrom,
not the current 0.5 Angstrom grid) failed near its 100-iteration budget when
ARPACK again returned only seven of eight requested modes. This older loaded
code had neither partial-mode recovery nor intermediate weight checkpoints.
The initial record, raw failure-log hash and last evaluation are preserved in
its failure companion JSON. No certified final result or unsaved best weights
are claimed. The current realistic-shift grid continues separately.

### Covariance-guided fits and restricted spectral exchange

All eighteen new experimental covariance-guided sensitivity rows are complete.
The first exposure half supplies a 20%-isotropic-shrunk second moment as a
fitting guide; the independent second half supplies the final variance bound.
All six fixed-pose widths improve over the preceding same-weights directional
audit. Five one-degree widths worsen and the sixth stays at no data. All
six two-degree cases now use no data, for nine fallbacks overall. Only the
same broad 10076 fixed-pose interval excludes zero. All six proxy fits have
objective gaps below 0.00321; those are NOT calibrated-noise objective gaps.

The first 10049 two-degree spectral-cut probe reduces relative width from
0.76117 to 0.69246 after 30 rounds. Its full-space dual gap remains 0.97307;
feasible stress reaches 0.45878 of the bias upper and reference power is below
1e-240. A second 30-round probe exporting the cubic cone's dual vectors reaches
width 0.68883, with its own gap 0.96929. Joint optimization over spectral support
mixtures and cubic dual balls reduces the original fit's gap from 0.99269 to
0.63110, and the follow-up's own gap to 0.93909. These postprocesses do not alter
weights or intervals. Different feasible lower bounds must not be confused
with changes in the primal estimator, or with a global optimum.

The full suite passed 73 tests after the conic-dual update; a subsequent
independent conic test also verifies joint pose-hull/cubic-ball projection.
The curated survey now has 97 entries after targeted reading of the accepted
cryoSENSE paper and supplement, with exact primary PDF hashes. Experimental
scope now explicitly names the shared-density assumption: allowing arbitrary
means in calibration does not extend inference to arbitrary structural mixtures.
The known heterogeneous nature of 10076 is documented from its primary archive.

### Full-weight inner projections and joint density/pose auditing

The rank-256 six-round full-weight conic trial finished in 201.50 seconds.
Its projected candidates failed the ambient comparison, so input weights
were retained. A fresh final spectral audit gives width 0.69024 and gap
0.76228. This is not an estimator improvement; adaptive continuous residual
enrichment is being tested separately. The initial dual-sign test failure
was corrected before this experiment and is preserved.

Six same-weight cross-term post-audits are complete. They account for the
relative directions of the nominal continuous residual and polynomial pose
field, while retaining the original spectral event and cubic remainder.
All six widths decrease; the 10076 half-degree broad-target conditional
reference power rises from 0.6703 to 0.9680. The two-degree restricted-cut
fit's width falls from 0.68883 to 0.59821, with saved feasible stress reaching
0.53939 of the new upper and reference power still below 1e-160. Sharper exact
cube fourth/sixth directional moments reduce its width further to 0.57162.
These are supplied-bound simulated-reference checks, not measured coverage
on experimental particles. The common-covariance, density and pose assumptions
remain unresolved. JOINT-DENSITY-POSE-BIAS.md gives the complete argument.

The expanded full suite passes 82 tests in 7.42 seconds. The Gaussian boundary
moment, Hilbert cross-term, adaptive PSD projection and nonlinear remainder
tests use independent calculations. A scalar-shape error in the remainder
test fixture was corrected; its failure log remains separate from the passing
run. No second full scientific Fable review has yet occurred.

The 30-round adaptive full-weight trial then completed in 720.99 seconds.
It retained the same input weights but raised the full-space lower bound to
9.200106, reducing the OLD triangle sum-objective gap to 0.181913. This lower
bound must not be used as a lower bound for the tighter joint-bias objective.
Current new inner candidates still have larger ambient objectives than the
input, so a follow-up checks recent-candidate convex averages with fresh
final spectral probes. Averages are only candidate-generation heuristics.

The full-cube 10049 support control also completed (6,324.04 seconds).
Its nominal sum-objective gap is 0.001867, but after pose auditing every
tested tail allowance uses the no-data interval; the cubic bias is 115.36.
This completes the three support-size sweep. Smaller-support improvements
must still be interpreted against the declared outside-density allowance,
and the original uncropped reference violates the smaller-support classes.

The 34-page revision was built without final-pass overflow/undefined-reference
warnings, and new cross-term formulas and the six-row result table were
rendered and visually inspected. Source/results checkpoint 1b91132 was pushed
publicly; a separate v0.4 archive packages completed numerical array owners
and excludes active/failed post-review checkpoints.

### Known-pilot refinement, larger-band fit, and ambiguity diagnostics

The 60-round adaptive full-weight average changes weights, reaches an old
triangle-objective upper/lower pair 10.5226/9.84465 (gap 0.06443), and has
relative width 0.64566. Joint and sharp-moment post-audits reduce this to
0.58103 and 0.55686. Using the actual known pilot pairing instead of only its
norm reduces it further to 0.51016, on the same density class, weights and
spectral event. Saved feasible nonlinear bias reaches 0.66665 of the last
upper bound, but reference sign power remains below 5e-122. The elementary
pilot-pairing derivation and analytic constant-cell moments are in
PILOT-POSE-PAIRING.md. All prior values remain available.

Focused Fable audits 03 and 04 prompted a real legacy-fallback regression,
explicit missing-center handling, expanded attempt accounting, selected-column
direct exponential-sum checks, and a source-bias guard. The archive-scale pilot
check has maximum relative discrepancy 3.4818e-13 across sixty columns per
case in the initial ten cases. This is not an all-column numerical error pad.
The original full-paper rejection remains unchanged.

The physical 10 Å, 1,024-particle fixed-pose fit on 10049 completes at a
sum-objective gap 0.000992 and relative half-width 0.013777. The reference-map
central value is 1.92452 and its conditional sign power is one. This is known-
noise/fixed-pose development; its one-degree, 0.5 Å pose post-audit is still
running and may remove the apparent precision. The 128- and 1,024-particle
subsets are not nested, so their difference is not an isolated count effect.

The wide-Ritz 10028 half-degree optimization was interrupted after its best
checkpoint at evaluation 131 plateaued near objective 2.37284. Its original
record, checkpoint and traceback are preserved. A new bounded run resumes
with maxfun 100, ftol 1e-5 and a fresh final certificate seed; its first case
has width 0.14610 but objective gap 0.66745. A valid final upper bound is not
convergence. Other source-grid and conic-continuation cases remain active.

A new prospective two-pose ambiguity protocol supplies constructive lower
bounds for arbitrary deterministic-length confidence intervals by classical
Gaussian testing. Independent small-Hilbert-space conic and direct-Fourier
checks pass; the full suite now has 98 passing tests. The first development
attempt hit a NumPy-boolean JSON serialization bug after its first fit. That
failed attempt and traceback are retained; the scalar conversion is fixed
and the identical prespecified grid reruns in `two-pose-modulus-v2`.
The grid and its focused mathematical review are in progress. Neither its
selected pose pairs nor any eventual numerical convergence should be called
a global optimization over the unknown-pose class.

The larger-band pose audit completed in 1,491.42 seconds using 1.08 GB peak
memory. At one degree and 0.5 Å it returns no data: the polynomial pose bias
is 6.3760, cubic bias 56.4927, final unselected half-width 61.2227 and no-data
half-width 34.3900. The target's Gaussian standard deviation is 10 Å; the
sampled Fourier band ends at 19.68 Å, so this is not a 10 Å reconstruction
resolution claim. The saved fixed-pose reference power of one does not survive
this pose bound.

The 60-round sharp-cubic conic continuation also completed, in 4,666.18 seconds,
with sum-objective upper/lower 9.89606/9.59349 (3.057% gap) and width 0.60696.
The joint and known-pilot post-audits give 0.54750 and 0.50264. The strongest
saved feasible stress reaches 0.63008 of the last upper; reference power is
9.32e-119. This is progress on a conditional bound, not practical detection.

The two-pose v3 grid is complete: thirty cases, all fixed-pair brackets below
0.4831%. The focused review confirms its real-arithmetic argument; its tests
now exercise inconsistent norms, nonzero quadrature pads, nonlinear pilot
projections, and distinct best-upper/lower iterates. All prior attempts remain.
The full suite passes 102 tests in 7.78 seconds. The two main pose-optimization
grids are still active; no full second-round acceptance assessment has occurred.

The enclosing-domain remainder refinement has now been implemented and checked
against independent spherical/cube integration of derivative tensors and direct
nonlinear pose fields, including distorted embeddings. Its ball probe reduces
the 1,024-particle cubic bias from 56.4927 to 32.6233 but still returns no data.
The subsequent expanded-cube probe reduces it to 21.1103, with half-width
25.8403 (0.75139 of no data), still zero reference sign power. Each probe takes
about five seconds and preserves the existing polynomial certificate. They do
not narrow the unknown density class. Source snapshots and both outcomes remain;
this is exploratory refinement, not a frozen confirmation. All 109 tests pass
in 6.48 seconds. The exact requested Fable model is conducting a focused proof
and code audit; no second full-paper review has occurred.

The enclosing-domain focused Fable audit is complete. Its conditional pose-set
objection does not invalidate the actual joint-ball source case, but revealed a
missing guard. New guards require the source pose/density class and recompute
both original bias branches before substitution. Guarded v2 replays reproduce
both bounds. New tests cover a near-tight closed-form translation, antipodal
frequencies, scaled detector embeddings, product directions, large angles and
finite differences. Six independent 60-decimal derivative-norm checks have
maximum relative discrepancy 4.501e-16. The suite passes 116 tests in 6.57 s.
These fixes leave the first full-paper rejection and usefulness limitations intact.

A bounded pose-optimized ambiguity pilot follows a written prospective protocol.
The 24/48/72-cell projections are feasible subclasses of the original continuous
ball, with the original pilot prolonged exactly. Local pose fitting improves
each projected pair but yields relative lower widths 0.09624/0.12169/0.12438,
all below the source continuous witness (0.12634). No larger fixed-density grid
was launched. The higher-resolution projected pair already has testing-distance
slack, so a subsequent continuous-density refit at those poses is a distinct,
prospectively specified next step. That first refit reaches 0.13884 (fixed-pair
bracket gap 0.2694%). Two further bounded alternations are running. None is a
global pose maximum. The witnesses contain substantial negative-density energy;
they demonstrate the declared signed class, not biological admissibility.
An initial unit-test mismatch arose because the old forward NUFFT uses a looser
requested tolerance than the analytic-moment path. The final distance is now
checked with direct physical-cell sums for every particle, with recorded NUFFT
discrepancies and magnitude pads. The current suite passes 118 tests in 6.45 s.

Two additional targeted primary readings, CryoPROS and CoCoFold, bring the curated
candidate list to 100. The exact XML sources and selected sections are recorded;
this is not a full-reading count. Their validation controls motivate acquisition-
intervention and structural-prior tests without supplying calibrated pose sets.

The bounded alternating ambiguity pilot has completed its three planned density
refits. Relative lower widths are 0.138839, 0.145941 and 0.149617, with final
fixed-pair bracket gap 0.1422%. No further alternation was launched. The manuscript
retains all projection losses and the signed-density limitation, and depicts the
first continuous witness. The latest 10076 contrast/two-degree optimization also
completed; its same-weight pilot refinement is 0.51237 of no data with zero
reference sign power. The completed 10028 center/two-degree refinement is
0.43581, also with zero power. These improvements do not establish practical
feature detection or global optimizer convergence.

A new pilot-only target-selection lock fixes three positive smoothed regions
per stack and a center control, all at 20 Å Gaussian standard deviation. It uses
only independent pilot coefficients, with at least 40 Å between the selected
regions. All 12 targets are retained. A subsequent written experiment allocates
a family error budget before any new reference outcomes; the first fixed fit
is running. This is a development utility check, not a new confirmation or a
10 Å resolution claim. Its optional experimental application is also specified
before feature outcomes, using independent exposure representatives and a
group-size variance inflation to allow arbitrary within-exposure Gaussian
dependence. Other experimental assumptions remain unverified. Target-lock,
score-replay and group-noise tests pass; no new target outcomes are claimed yet.
