# Revision evidence following full review round 1

**Development checkpoint, not a completed resubmission or an acceptance claim.**
The unmodified first Fable 5.1 review remains the latest full acceptance
assessment: reject, confidence 4/5, not a strong contender. Focused mathematical
audits are not full reviews. All twelve fixed-pose feature fits, nonlinear pose
audits and the frozen fresh noise-calibration application are complete. Their
results remain conditional and do not resolve experimental calibration. The eighteen-case
pose-optimization grid and its post-audits are complete.
This file also supports local packet-size/layout checks without invoking a model.

## Registered-reference follow-up, 30 September

A pilot-only frame registration was frozen before the eleven completed reconstruction FSC
comparisons and before the subsequent target sensitivity. That sensitivity
retains all original results: the unchanged cubic estimator's simulated power
on the previously selected 10049 region rises from .00654 to .95725. This
supersedes the interpretation of the old generator as correctly located
relative to the consensus, not its numerical output. The subspace fit remains
the same estimator and the large optimization gap remains unresolved.
All twelve fixed-pose target powers are at least .9563 after registration;
the original one-/two-degree fixed-weight audits remain uninformative.

A separate predeclared application reuses the existing experimental calibration
pool for all three cubic estimators. None excludes zero. The original cubic
interval is centered at 2.46178 with half-width 4.30107, and its noise SD bound
is 4.55146 times the simulation SD. The registered approximate reference value
2.48846 lies inside, which is not a density coverage label. Thus the frame audit
corrects a material interpretation error but does not resolve R1/R2/R7 or
constitute an acceptance assessment. Details and all alternatives are in
REGISTERED-TARGET-SENSITIVITY-RESULTS.md. Earlier paragraphs below retain the
historical reference-frame outcomes with this explicit superseding context.

A subsequent deterministic Helmert-contrast calibration uses 127 rather than
128 independent Gaussian rows and permits arbitrary means. It reduces the
original cubic half-width to 3.07893, which still contains zero; fixed/shift-only
exclusions rise to ten/six, rotational exclusions remain zero. This is a
classical Gaussian projection argument and reused-data sensitivity, not a
new coverage principle or evidence that experimental assumptions hold. Five
centered fixed-weight intervals exclude the registered reference, all on
heterogeneous 10076 (regions 1 and 3 at fixed/shift-only poses, center at fixed
poses). These discrepancies are retained; a single Class A reference is not
a ground-truth coverage label for that population.

The seven-contrast projected-calibration follow-up allocates beta/7 to each
noise-independent design and retains every result. Its SD bounds are .981–1.030
times the separate mean-only bounds; counts and five reference discrepancies
are unchanged, and all three cubic intervals become slightly wider. This is
negative calibration evidence, not a resolution of R1.

Two unknown-pose RELION continuations now finish with verified exposure halves.
The 10028 registered-reference mean FSC is .702; 10049 is much weaker at .202,
despite reporting convergence. Its native half FSC crosses .143 at 18.95 A.
Every curve and the earlier timed-out initializers remain. The 10076 continuation
is running. These are conventional reconstruction comparisons, not density
coverage or calibrated pose-radius evidence.

The adaptive triangle-objective enrichment completes in 3,476.62 seconds.
Its width improves by only .79 percent to .178365 of no data; its full-space
gap remains .927173. Registered-frame simulated power is .964129, with no
experimental application yet. A separately declared joint density/pose design
starts from the original weights and thirteen-column basis and is running. Its
residual-controlled trust-region proposition and independent small robust-norm
SDP agree, but no empirical improvement is claimed. A focused Fable audit of
this mathematics finds no validity error and seven numerical/design notes.
Its corrections precede the first empirical joint fit; it is not full review
round 2. The latest full regression run gives 246 passed, one skipped and one
intentional warning;
all 116 locked fresh-calibration files and 47 earlier source files verify.

## Major concerns

**R1 — Experimental inputs: still open.** We implemented the noncentral Gaussian
trace bound, directional calibration for a fixed estimator, and an extension
allowing arbitrary jointly Gaussian dependence within an inference exposure.
The latter inflates each weight row by the square root of its group size.
These are elementary conservative results, not data-driven proofs of a common
covariance or supplied-pose independence. Completed exploratory applications
are in experimental-noise-grouped, directional-noise-audit and noise-metric-design.
In those early applications only one broad fixed-pose feature excludes zero;
none at one/two degrees.
The already inspected second calibration pool is explicitly exploratory reuse.
The new cohort reserves 128 unused exposures per stack. All twelve estimators
and 116 files were locked and published in commit 60efd9b before downloading
any reserved pixels; all three downloads and the frozen recalibration are complete.
The full twelve-feature family now has six fixed-pose and four shift-only zero
exclusions, but none at one/two degrees. Fresh noise SD bounds are 0.729–0.955
of the original bounds; the exclusion counts are unchanged. All 48 intervals,
including sixteen no-data fallbacks, are retained. Old inference images
remain development data. The all-row metadata inventory finds zero-valued
pose/shift ESS fields throughout; these are not calibrated zero-error bounds.
The paper explicitly identifies 10076's heterogeneous population as incompatible
with treating homogeneity as experimentally established. No end-to-end coverage
claim follows from the new covariance algebra or deposited-map inclusion.

A separate known-map information calculation retains all 768 particle/map
cases at the simulation's limited radius-twelve band. Its local variance
scales do not calibrate pose balls or diagnose deposited full-data alignment
errors. Targeted moment and orbit-likelihood readings also make clear that
our negative local-pose audit is not a general reconstruction impossibility.

**R2 — Useful fine-scale inference: still open.** The 1,024-particle/radius-12
10 A sigma fit is precise at fixed poses but fails sign detection under a
one-degree/0.5 A audit. An enclosing-cube remainder improves its width to 0.751
of no data, without restoring power. The continuous support ablation retains
outside-class references; the smallest masks exclude them. A further exact
support-function/Jensen sign-class ablation also retains all excluded-reference
cases and all grid sizes. It does not rescue the inherited pose remainder.
The new locked sigma-20 A target study selects three pilot regions and one
matched center per stack, before evaluating their feature outcomes. All twelve fixed-pose fits
are complete, with maximum relative sum-objective gap 0.00437. The third 10049
region has fixed-pose reference sign power 0.0158; all other features exceed
0.99999. Nonlinear audits and the full experimental application are complete.
Every one-/two-degree audit lacks useful reference sign power; 10028 uses no
data at one degree and every stack does so at two degrees. A declared two-case higher-order
remainder diagnostic reduces one component, but is not a new interval. The
separate full cubic audit has completed: it reduces one selected width from
0.473 to 0.252, yet minimum reference sign power remains 1.48e-10. Its actual
polynomial/density contribution remains too large for useful inference. A new
single-case cubic weight-optimization protocol and code were committed before
fitting; prerequisite gradient/support checks and an independent conic solve
passed. That empirical fit completes in 7,191 seconds, reducing the relative
width to 0.180. Minimum reference sign power remains 0.00654; the iteration
limit and 0.940 relative surrogate gap remain explicit failures of usefulness
and tight optimization. The final-weight two-tolerance checks are stable but
do not validate all operator error. Its gap triggers the previously declared
coordinate-metric follow-up from the same original weights; that fit completed without improvement: relative width 0.223, minimum reference
sign power 7.73e-7, relative surrogate gap 0.977 and no convergence at the
thirty-iteration budget. A separate declared thirteen-dimensional convex
design is complete and selected the original weights to relative difference
4.61e-16; it supplies no estimator improvement. The fresh full-space gap remains
0.940164.
Focused Fable mathematics checks and their responses do not constitute a new
acceptance review. Its sigma
is not its sampled frequency-band resolution. These tests do not satisfy the
review's requested fine-scale usefulness criterion yet.

**R3 — Pose-aware weights and bound looseness: partially addressed, not resolved.**
There is now a continuous matrix-free pose objective, numerical spectral upper
certificate, joint residual/pose cross term, sharper cube remainder, and exact
pilot-cell moment pairing. All comparisons retain failed/nonconverged fits and
their gaps. The eighteen-case grid is complete, including nonconverged fits. Its six
two-degree pilot-refined widths are 0.427--0.633 and feasible/upper bias ratios
0.274--0.360; the largest minimum-reference power is 3.27e-72. This does not
meet the reviewer's all-stack usefulness or bound-tightness criteria. The summary
also retains three separately declared 10049 follow-up optimizations. A full-weight two-degree probe reaches
a relative width around 0.503 but has negligible sign power. Feasible adversary
ratios and optimization gaps are distinct and are both reported. The coordinate
frame is declared externally: coherent rotations belong to the specified
sensitivity class. We do not delete that mode while keeping the same target.
All joint versus product rotation/translation balls are labeled explicitly.

The additional two-pose construction supplies a deterministic-length lower bound
on the same continuous class. The 30-case revised grid has fixed-pair brackets
below 0.5 percent. Three prespecified alternating pose/density refits raise one
selected feasible lower width from 0.1263 to 0.1496 of no data. These signed
witnesses are not asserted to be molecular structures, a global pose search,
or sharp limits for variable-length procedures. Classical testing/duality is
credited; this is an empirical diagnostic of the remaining gap.

**R4 — Scaling: a completed computational demonstration.** The matrix-free
10,000-particle/radius-12 audit completes on the Mac at 1.95 GB peak memory.
The randomized upper event, quadrature errors, exact source snapshots and
small dense-operator checks are preserved. Ordinary floating-point/NUFFT guards
are distinguished from validated interval arithmetic. This resolves the earlier
dense-storage obstacle for this run; its broad feature still has a wide interval.
Scaling alone does not address R1, R2 or R7.

**R5 — Verification versus calibration: rewritten.** The abstract states the
conditional simulation assumptions and does not present record counts as
experimental calibration evidence. Analytic coverage and original Monte Carlo
cross-checks remain in the immutable records as implementation verification;
widths and sign power carry the usefulness conclusions. No extra simulation
count is advanced as evidence of scientific generalization.

**R6 — Delivered method: reframed.** The title and core framing now identify
continuous Fourier-slice auditing. Gaussian dictionaries are numerical tools
and motivating representation controls, not a claimed new reconstruction family.
The fresh-exposure prediction comparison remains supplementary context, explicitly
distinguished from density uncertainty. It remains favorable to the stock neural
baseline in all six prespecified prediction contrasts.

**R7 — Novelty: remains contingent on substantive utility.** No new general
statistical coverage principle is claimed. Optimal recovery, folded-normal
critical values, two-point testing, support functions and elementary concentration
arguments are attributed or identified as standard. The operator-specific
continuous integration, scalable nuisance bounds and diagnostic implementation
must be judged on their actual value; extra theorems, tables or passing tests
do not by themselves answer this concern.

**R8 — Nuisance scope: expanded sensitivities, not calibrated inputs.** New pose
studies use 0.5 A translations with the stated joint-ball convention. The code
supports particle-specific radii, but no archived per-particle confidence radii
have been established. Uniform defocus/astigmatism/phase/gain/B-factor envelopes
are evaluated at the higher band. Their negative results and conservatism are
retained. Consensus-pose dependence, heterogeneity and gain/density calibration
remain scientific limitations.

## Separate likelihood-validation development

A full-orientation structural-compatibility candidate applies established
universal inference and Lindsay mixture duality. Its first discrete-view,
oracle-numerator screen detects full local removal on two of three geometries;
all half-removal failures and a third-geometry failure remain. Shared noise
profiling preserves similar counts. This is neither a practical learned test
nor a new general statistical principle. Three continuous-SO(3) computations
and two diagnostics are now complete, but their best recorded global gaps
remain 17,077--20,253 log units. They do not resolve the practical objection.
The declared envelope-guided refinement completes but leaves continuous gaps
of 12,898--15,892 log units. The subsequent reviewer-motivated disk diagnostic
does not improve any of the 285 sampled coarse cells.
The new paper appendix includes these outcomes and explicit assumptions; R1,
R2 and R7 remain open. Source and failed-attempt provenance is retained.

## Baseline fairness and smaller items

The expanded Fourier Gaussian baseline independently implements Ullrich et al.'s
equation-14 objective with explicit project extensions; it does not pretend to
run the original full pipeline. Exact full and diagonal posteriors are compared
under the same finite Hermitian model. Diagonal intervals are wider on the
reported targets. Matched prior coverage, fixed-reference coverage, interpolation
error and continuous same-estimator auditing remain different results.

For the locked targets, the original prior sweep had poor fixed-reference
coverage at several features. A separately declared broader-prior sensitivity
improves nominal full-posterior coverage to 0.9872–0.9998 with relatively narrow
intervals. Every broader-prior result is retained; no favorable prior replaces
the original record. This is meaningful favorable baseline evidence. The new
continuous feature family is now complete and is reported alongside it, without
claiming superiority from its different uniform-class guarantee.
The favorable wider-prior behavior also persists in the declared one-/two-degree
coherent and random-boundary checks: full-posterior fixed-generator coverage at
coordinate SD one ranges from 0.9861 to 0.9999. All five priors and every pose
pattern remain reported; these selected-pattern probabilities are not uniform
pose coverage, and some fall below the 0.995833 nominal marginal level.

A separately declared local Gaussian pose-marginalization comparison is now
complete on the same twelve features and both broader priors. All 48 solves
complete and fixed-pose controls exactly reproduce their original widths.
Pose marginalization widens these local-model intervals by 0.0066–0.2875 percent;
favorable and unfavorable reference probabilities remain reported. Its pilot
linearization error at one degree is 0.71–2.41 percent of signal norm, and it
omits density/pose products. This is a classical local Gaussian control, not
the full Rangan Hessian method or an independently calibrated pose posterior.
It does not resolve R1, R2 or R7.

- M1: the paper is a public development preprint, not an anonymized submission.
- M2: the current curated count is 119 candidates; historical counts and targeted
  reading versus retrieval are distinguished.
- M3: the manuscript uses the source-summary value 0.339 for the strong-ridge
  bootstrap control, alongside all favorable weak-ridge comparisons.
- M4: the revised critical-value implementation removes the asymptotic shortcut;
  frozen old code and an audit of its negligible numerical effect are retained.
- M5/M6/M7/M11: endpoint neural checkpoints, physical shift units, the meaning
  of signed witnesses and supplied simulation assumptions are explicit.
- M8: Cai–Low is cited for the relevant modulus/expected-length discussion;
  citations must not be added as substitutes for checking their precise claims.
- M9: the user's requested broad survey is separately maintained. The development
  appendix remains long and will need an editorial pass for a final submission.
- M10: new runners record outcome/failure status rather than aborting before
  saving scientific failures. Frozen historical evaluators remain reproducible.
- M12: the round-1 rendered main text ended on page 7, despite the source-only
  review's concern. Future full packets include rendered pages and figures.

No proposed or running experiment in this response is represented as completed.
Before a full new review, update this response, finish the declared studies,
rebuild/inspect the paper, and attach exact completed-artifact provenance.


The original-code CryoLike comparison completes all 24 cases on the three
experimental stacks, retaining both viewing grids, both metrics and all paired
contrasts. Its rankings differ by metric. This supplies an external scoring
comparison but does not resolve R1/R2/R7 or provide matched uncertainty power.
The focused mixture audit and response remain separate from the full review.


A new reference-frame audit, motivated by the first converged RELION run,
finds a large registration effect on 10028/10049. Three transforms were selected
only against the old pilot and published before comparing inference maps.
All original curves remain. On 10076, the Class A reference is not truth for the
heterogeneous consensus and registered agreement remains weak. Historical
reference-generator power/coverage numbers are conditional studies in their
original frame, not biological localization evidence. The subsequently completed registered feature sensitivity is summarized at the
top of this response; these changes do not close the experimental-input objection.
