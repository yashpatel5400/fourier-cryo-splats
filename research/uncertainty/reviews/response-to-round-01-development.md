# Revision evidence following full review round 1

**Development checkpoint, not a completed resubmission or an acceptance claim.**
The unmodified first Fable 5.1 review remains the latest full acceptance
assessment: reject, confidence 4/5, not a strong contender. Focused mathematical
audits are not full reviews. All twelve fixed-pose feature fits are complete; nonlinear pose audits and the
fresh noise-calibration application are running;
they must be completed and incorporated before the next full review. The eighteen-case
pose-optimization grid and its post-audits are complete.
This file also supports local packet-size/layout checks without invoking a model.

## Major concerns

**R1 — Experimental inputs: still open.** We implemented the noncentral Gaussian
trace bound, directional calibration for a fixed estimator, and an extension
allowing arbitrary jointly Gaussian dependence within an inference exposure.
The latter inflates each weight row by the square root of its group size.
These are elementary conservative results, not data-driven proofs of a common
covariance or supplied-pose independence. Completed exploratory applications
are in experimental-noise-grouped, directional-noise-audit and noise-metric-design.
Only one broad fixed-pose feature excludes zero; none at one/two degrees.
The already inspected second calibration pool is explicitly exploratory reuse.
The new cohort reserves 128 unused exposures per stack. All twelve estimators
and 116 files were locked and published in commit 60efd9b before downloading
any reserved pixels; all three downloads are complete and verified. The frozen
calibration application awaits the complete original twelve-feature results. Old inference images
remain development data. The all-row metadata inventory finds zero-valued
pose/shift ESS fields throughout; these are not calibrated zero-error bounds.
The paper explicitly identifies 10076's heterogeneous population as incompatible
with treating homogeneity as experimentally established. No end-to-end coverage
claim follows from the new covariance algebra or deposited-map inclusion.

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
0.99999. Nonlinear audits and the full experimental application remain pending. The first
region's completed one-degree audit has no useful sign power on any of the three
stacks, including one no-data fallback. A declared two-case higher-order
remainder diagnostic reduces one component, but is not a new interval. The
separate full cubic audit has completed: it reduces one selected width from
0.473 to 0.252, yet minimum reference sign power remains 1.48e-10. Its actual
polynomial/density contribution remains too large for useful inference. A new
single-case cubic weight-optimization protocol and code were committed before
fitting; prerequisite gradient/support checks and an independent conic solve
passed. That empirical fit is running, without an asserted successful outcome.
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
continuous feature family must be compared with it after completion.
The favorable wider-prior behavior also persists in the declared one-/two-degree
coherent and random-boundary checks: full-posterior fixed-generator coverage at
coordinate SD one ranges from 0.9861 to 0.9999. All five priors and every pose
pattern remain reported; these selected-pattern probabilities are not uniform
pose coverage, and some fall below the 0.995833 nominal marginal level.

- M1: the paper is a public development preprint, not an anonymized submission.
- M2: the current curated count is 105 candidates; historical counts and targeted
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
