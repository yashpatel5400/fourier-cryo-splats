# Independent review record

Round 1 is complete through the authenticated Claude CLI using the exact
`claude-fable-5-1` model. The unmodified review says **reject**, confidence 4/5,
and **not a strong contender for ICML acceptance**. It checked the supplied
mathematics favorably but identified major gaps in usefulness, experimental
calibration, pose-aware optimization and scaling. See [the full review](round-01/review.md).
The raw provider metadata confirms the requested model; no proxy reviewer was used.

After the complete frozen studies and current manuscript have been integrated,
`scripts/review_uq_candidate.py --round 1 --invoke` creates an immutable review
packet and invokes that exact model. It saves the original prompt, every input
file's hash, Git status, PDF hash, command, raw response, stderr and model-usage
metadata. Tools/MCP/local customizations are disabled, and the reviewer receives
the self-contained manuscript source, proofs, all project Python scripts/modules,
tests, protocols, critical
literature notes and result summaries. It is explicitly told that figures are
not supplied as rendered images and that it must state verification limitations.

The review requests a critical verdict with no desired outcome, stable concern
IDs and a finite prioritized revision plan. It may reject the work. Neither a
favorable verdict nor repeated reviews establish actual conference acceptance.
No alternative model may be represented as Fable 5.1.

Subsequent rounds require a response document and include every preceding
unaltered review. Responses should address each concern with evidence, distinguish
implemented fixes from argued disagreements, and identify unresolved issues.
New experiments must preserve previous outcomes and state whether their designs
were motivated by feedback. Never omit an unfavorable review, rewrite its verdict,
or ask for approval regardless of scientific merit.

The research goal is not complete merely because a review was requested or a
fixed batch of computations finished.

Four subsequent focused mathematical audits are separate from full review
rounds. `bound-audit-01` contains an authentic but incomplete terminal response;
its missing beginning is disclosed in `COMPLETENESS.md`. The complete
`bound-audit-02` checks the cross-term, moment and conic refinements and flags
the essential joint-pose-set definition plus implementation/reporting issues.
Its response plan links the fixes and limitations. Neither audit is an ICML
acceptance assessment. Future invocations preserve every assistant text message
in the event stream, including any continuation, rather than only the terminal
result string.

`bound-audit-03` found a legacy fallback broadcasting bug, which is fixed and
checked by replaying an actual archived case. `bound-audit-04` confirms the
known-pilot moment refinement and requests an updated inventory, an archive-scale
direct-sum numerical check, and a fixed-pose source guard. Its response records
the fixes and explicitly limits the numerical check to selected columns.
None of these reports changes the full round-1 verdict.

The separate `modulus-audit-01` checks the later two-pose ambiguity lower bound.
It confirms the real-arithmetic construction and testing argument, then requests
failure guards, exact-continuous integration-pad tests, distinct upper/lower
witness records and clearer numerical disclosures. Its evidence-linked response
preserves the original report and distinguishes that check from a full review.

`enclosure-audit-01` independently checks the later Fourier-cancellation Taylor
remainder. Its report confirms the identities but flags missing source-class
and bias-decomposition guards plus test gaps. The response documents the new
fail-closed checks, unchanged guarded replays, deterministic and near-tight
controls, and six independent 60-decimal spot checks. It is also separate from
full review rounds and does not alter the rejected scientific assessment.
