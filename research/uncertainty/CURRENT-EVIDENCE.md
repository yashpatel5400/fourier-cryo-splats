# Current evidence and unresolved research decisions

1 October 2026 UTC, post-round-2 development checkpoint. The latest full independent review
is **reject**, confidence 4/5; the work is not yet an acceptance-level result.
Historical protocols, failed attempts and earlier source snapshots remain
unchanged. This page identifies current outcomes rather than replacing them.

| Question | Completed evidence | Practical limit |
|---|---|---|
| Does a finite reconstruction dictionary hide uncertainty? | Continuous adjoint-residual checks and explicit excluded directions; matched finite-model Bayesian and bootstrap controls | The strongest excluded directions need not resemble molecular density. |
| Can nonlinear pose terms be audited at useful scale? | A 10,000-particle quadratic audit completes at 1.95 GB; eighteen pose-aware designs and higher-order follow-ups are retained | Scaling is demonstrated for that audit, not every cubic optimizer or high-resolution reconstruction. |
| Does weight optimization improve the selected cubic case? | Original width .179785 of no data; coordinate .222993; fixed-subspace .179695; enrichment .178365; joint .175192 | Enrichment/joint gains are only .79%/2.55%. Enrichment triangle full-space gap is .927173; the joint restricted guide gap .000885 does not certify full-space convergence. |
| Did the map frame affect simulated usefulness? | Pilot-only registration changes original cubic minimum sign power .00654 to .95725 without changing weights | This repairs the interpretation of the generator location; it is not a new estimator or experimental coverage finding. |
| Do the intervals transfer to experimental pixels? | All five cubic designs were applied using raw and centered calibration; all ten intervals contain zero | Original centered interval 2.46178 ± 3.07893; joint 2.46896 ± 3.01088. SD bounds include signal energy and conservatism; they are not pure-noise measurements. |
| What did fresh calibration establish? | Twelve estimators and 116 files were locked before 128 new calibration exposures per stack; all three downloads/applications finish | The raw procedure excludes zero for 6 fixed-pose and 4 shift-only features; none with rotations. It does not establish common covariance, independent metadata or valid density/pose radii. |
| What changes with centering or more contrasts? | Centering gives 10/6/0 exclusions and five Class A reference disagreements on 10076. The seven-contrast family gives the same counts and slightly worse cubic widths | Reused-data development; alternatives are not combined by an unadjusted minimum. Class A agreement is not a truth label for a heterogeneous consensus. |
| Were reconstruction baselines actually run? | Supplied-pose Gaussian/voxel/stock-neural fits on all three stacks; unknown-pose RELION finishes on all three; 24 original-code CryoLike scores finish | RELION registered-reference mean FSC is .702/.202/.122, with reported convergence on all three. The 10076 reference is one Class A assembly state. FSC and scoring are distinct from density coverage. |
| Does global pose marginalization currently solve validation? | Discrete-view oracle screens and continuous-mixture bounds with retained failed/nonconverged outcomes | Continuous global gaps remain 12,898–15,892; a practical learned independent predictor is absent. |
| Is the survey comprehensive enough to delimit novelty? | A 119-candidate ledger, targeted primary readings, source/version/access records, direct Bayesian pseudo-atom and selection-bias prior art | Candidate retrieval is not full reading. No novelty claim follows from not finding a matched paper. |
| What has independent review established? | Authentic full Fable 5.1 reviews 1 and 2 both reject at confidence 4/5; 246 numerical tests pass, one skips | Focused audits and passing tests do not establish scientific utility or acceptance. Round 2 finds no end-to-end coverage experiment and an inadequately matched continuous Gaussian-prior baseline. |

All currently declared fitting and application runs are complete. The
55-page manuscript and 128-member artifact increment are published as
v0.6.2-dev. Round 2 is complete and preserved unchanged. Its priority is a frozen study that re-estimates poses from every noisy replicate, with matched continuous Gaussian-prior and mixed pose-error comparisons. A mixed model must explicitly state the independence assumptions; substituting random pose errors does not automatically justify inference with estimated poses. No new coverage result exists yet. Further calibration-contrast or single-feature optimizer variants are not planned.

Detailed current reports: [round-2 response plan](reviews/response-to-round-02-development.md),
[enrichment](CUBIC-ENRICHMENT-RESULTS.md),
[joint design](JOINT-CUBIC-DESIGN-RESULTS.md),
[enriched/joint experimental application](JOINT-ENRICHED-EXPERIMENTAL-RESULTS.md),
[experimental target/frame checks](REGISTERED-TARGET-SENSITIVITY-RESULTS.md),
[centered calibration](CENTERED-NOISE-CALIBRATION-RESULTS.md),
[projected calibration](PROJECTED-NOISE-CALIBRATION-RESULTS.md), and
[RELION 10049](RELION-10049-CONTINUATION-RESULTS.md) and
[RELION 10076](RELION-10076-CONTINUATION-RESULTS.md).
