# Current evidence and unresolved research decisions

1 October 2026 UTC, post-round-3 development checkpoint. The latest full independent review
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
| Is the survey comprehensive enough to delimit novelty? | A 131-candidate ledger, targeted primary readings, source/version/access records, direct Bayesian pseudo-atom and selection-bias prior art | Candidate retrieval is not full reading. No novelty claim follows from not finding a matched paper. |
| What has independent review established? | Authentic full Fable 5.1 reviews 1, 2 and 3 all reject at confidence 4/5; 279 numerical tests pass, one skips, with one expected conic warning | Focused audits and passing tests do not establish scientific utility or acceptance. The requested simulation and matched baseline are now complete; the proposed pose intervals fail the practical informativeness criterion. The third review also rejects and prioritizes alignment-bias reanalysis. |

The pre-round-2 fitting and application runs are complete. Their historical
manuscript and 128-member artifact increment are published as v0.6.2-dev.
The revised study is complete: all 600 newly simulated datasets finish, with 12,000 converged weight solves and all 216 raw/selected procedure cells retained. Every cell has observed coverage 200/200, including the simpler fixed-pose and Gaussian baselines; individual exact 95% Monte Carlo intervals are [.9817, 1]. All three pose-audit procedures fall back to no data on every trial. Raw median relative widths span 642–1,242, and nonlinear remainders contribute 99.31–99.53% of every raw mixed-zero width. The coarse local optimizer worsens near-truth poses, with rotation RMS medians 12.0–18.5 degrees. This does not establish a coverage advantage, useful pose-aware inference, or performance of production refinement. The main paper reports these outcomes directly.

Round 2 is preserved unchanged. The mixed-pose theorem still requires conditional centering that local alignment does not establish. Further calibration-contrast or single-feature optimizer variants are not planned.

Detailed current reports: [round-2 response plan](reviews/response-to-round-02-development.md),
[enrichment](CUBIC-ENRICHMENT-RESULTS.md),
[joint design](JOINT-CUBIC-DESIGN-RESULTS.md),
[enriched/joint experimental application](JOINT-ENRICHED-EXPERIMENTAL-RESULTS.md),
[experimental target/frame checks](REGISTERED-TARGET-SENSITIVITY-RESULTS.md),
[centered calibration](CENTERED-NOISE-CALIBRATION-RESULTS.md),
[projected calibration](PROJECTED-NOISE-CALIBRATION-RESULTS.md), and
[RELION 10049](RELION-10049-CONTINUATION-RESULTS.md) and
[RELION 10076](RELION-10076-CONTINUATION-RESULTS.md).

Post-review-2 additions: [matched continuous Gaussian derivation](CONTINUOUS-GAUSSIAN-EQUIVALENCE.md), [mixed-pose assumptions and proof](MIXED-POSE-DERIVATION.md), and [frozen end-to-end protocol](END-TO-END-LOCAL-POSE-PROTOCOL.md). The corrected manuscript states the prior equivalence and reporting limitations. These additions received the third full review, which rejects.

The 384-dataset local pose calibration and all three 200-replicate test batches are complete; see the [full refitting report](END-TO-END-LOCAL-POSE-RESULTS.md). Reporting additions are complete: [registered dictionary replay](REGISTERED-DICTIONARY-RESULTS.md), [density/noise breakdowns](BREAKDOWN-RADIUS-RESULTS.md), and a [phase-only independence control](PHASE-SPLIT-CONTROL-RESULTS.md). The latter confirms the familiar distinction between independent measurement noise and alignment attenuation; it is not a full cryo-EM method. The Gaussian radius-12 rerun is complete: all 48 fits converge, with all 672 conditional scenarios retained. Matched-scale widths are nearly identical on two stacks and 34–35% larger on 10028; the minimum tested coverage across both prior scales is .995970. The original four nonconverged fits and incomplete flag are preserved. See [matched comparison](CONTINUOUS-GAUSSIAN-V2-RESULTS.md). The focused manuscript now includes all 600 refitting outcomes. A [raw-movie acquisition pilot](RAW-MOVIE-PILOT-RESULTS.md) finds spatially correlated frame differences but does not establish pure-noise independence. All five [v0.7.0-dev numerical bundles](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.7.0-dev) and the rewritten 26-page PDF are published, with independently verified archive/member hashes and matching remote asset digests. Page QA and the third full review are complete. Its unchanged rejection and the [response plan](reviews/response-to-round-03-development.md) are public. The new [post hoc reanalysis](REFITTING-REANALYSIS-PROTOCOL.md) uses saved arrays to diagnose alignment bias and exact realized class envelopes.

The [post hoc replay](REFITTING-BIAS-REANALYSIS-RESULTS.md) is complete: 31,200 rows, 18 converged true-pose fits, exact reproduction of saved centres and all 2,400 realized class envelopes. R13 is now reported in the main paper. R14 remains a construction problem: actual class sensitivity and bound slack both matter. The smaller-radius diagnostics are all outside the full density class, which is explicitly flagged. No fourth full review has occurred.

The [direct folded-width search](FOLDED-RIDGE-REVIEW3-RESULTS.md) completes all
six known-pose targets and 119 converged linear solves, with only .062--.098%
width reduction. All finite-range gaps are below .005; several global gaps
remain loose. The updated 25-page manuscript reports this classical ridge
comparison without changing the frozen noise trials.

The [paired-exposure candidate gates](paired-power-v1/METHOD-GATES-RESULTS.md)
are also complete: 36 larger-view diagonal witnesses, 18 compressed matrix
cases, three full-frequency removal fits plus controls, twelve common-shift
fits, and 80 adversarial pose searches. The latter find 42 positive verified
violations on the strongest stack; four worst witnesses independently replay
within 8e-14. The new authentic focused Fable audit prompted these checks,
not a favorable acceptance verdict. The saved candidate directions do not
satisfy the required continuous-pose condition. Experimental covariance,
paired-frame mean agreement and practical density calibration remain open.
