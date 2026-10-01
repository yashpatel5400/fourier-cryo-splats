Round 3 is complete: **reject, confidence 4/5, not a strong ICML contender**. The [unchanged report](round-03/review.md), raw response, exact prompt and evidence hashes are preserved. The exact model is `claude-fable-5-1`; all 26 pages were visible and all 2,437 evidence copies were unchanged. The reviewer executed no tests and read no third-party papers or numerical arrays. The [response plan](response-to-round-03-development.md) prioritizes saved-array alignment-bias analysis and diagnosis of the structurally loose pose bounds before new experiments.

The later [paired-statistics focused audit](paired-statistics-audit-01/review.md)
checks Gaussian/cone algebra but identifies serious missing information and
continuous-pose gates. Its [response](paired-statistics-audit-01/response.md)
records full-frequency bounds, common-shift checks and independently verified
continuous-pose counterexamples. This is not a fourth full ICML review or a
favorable acceptance assessment.

Round 2 is now complete: **reject, confidence 4/5, not a strong ICML contender**. The [unaltered report](round-02/review.md), raw response/event stream, prompt, evidence index and manifest are preserved. The provider reports exact `claude-fable-5-1`, and the runner found no changed evidence. All 55 pages were sent; the reviewer states that pages 1–19 were not visible. It read code but executed no tests and read no third-party papers. These are limits of this review, not evidence that the omitted material passed. The [response plan](response-to-round-02-development.md) records the substantive next work.

# Independent review record

Round 1 is complete through the authenticated Claude CLI using the exact
`claude-fable-5-1` model. The unmodified review says **reject**, confidence 4/5,
and **not a strong contender for ICML acceptance**. It checked the supplied
mathematics favorably but identified major gaps in usefulness, experimental
calibration, pose-aware optimization and scaling. See [the full review](round-01/review.md).
The raw provider metadata confirms the requested model; no proxy reviewer was used.

The original round-1 packet contained manuscript source, proofs, project code,
protocols, literature notes and result summaries. Tools, MCP and local
customizations were disabled. Rendered figures were not supplied, and that
limitation remains in the unmodified review. Later focused visual preflights
verified the requested model's image input; subsequent full packets include
every manuscript page, uncropped, at the explicitly recorded render resolution.
The runner supports 150 dpi (default) and 125 dpi for the larger revision packet.

The revision runner saves the original prompt, input hashes, Git status, PDF
and rendered-page hashes, exact command, raw response, stderr and provider
model-usage metadata. `--preview-dir` constructs a local packet without invoking
a model or consuming a full-review round. An actual revision uses a new round
number, a response document, and `--invoke`; the completed round 1 must never
be overwritten or rerun in place.

The optional `--read-only-evidence` mode keeps manuscript, notes and
summary-level outcomes in the initial text. All implementation modules, tests,
detailed case records and runner scripts remain available in an exact, checksum-indexed `evidence/`
copy, together with referenced archived source snapshots. Selection uses file
roles, not scientific outcomes. Only Read, Glob and Grep are enabled, confined
to the review directory; shell execution, writes, web access, MCP and local
customizations are disabled. The runner verifies the evidence hashes again
after the response. Merely making a file available is not evidence that the
reviewer inspected it. The prompt requires the reviewer to state that limit.
This mode was used for full round 2.

Current manuscript sources are discovered from `main.tex` and its literal
inputs/bibliographies. Historical reconstruction drafts remain indexed and
available in the evidence copy, explicitly distinguished from current claims.
All non-summary result records, including top-level case records, are deferred
by document role rather than outcome. The page PNGs may be recompressed
losslessly for transport; every decoded pixel, color mode and dimension is
checked unchanged, and both original/rendered byte hashes and pixel hashes are
recorded. No page is cropped or omitted; the chosen DPI is recorded explicitly.

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
None of these focused reports changes a full-review verdict.

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

`cubic-audit-01` and `cubic-design-audit-01` check the higher-order field and
its weight-design objective. They support the real-arithmetic constructions
while identifying numerical and test gaps. Their unchanged reports and
evidence-linked responses remain separate from acceptance reviews. The final-
weight two-tolerance numerical check is complete and stable;
an operator tolerance check is not a validated global arithmetic bound.

`mixture-audit-01` is a complete focused audit of the separate pose-marginal
likelihood development. The exact requested model supports the real-arithmetic
identities, identifies statistical caller assumptions and coarse-cell slack,
and suggests two computational diagnostics. Its response distinguishes confirmed
issues, a clarified result-label misunderstanding, implemented tests and work
still pending. It is not a full-paper acceptance assessment. The later full round 2 remains a rejection.

The focused `joint-trust-audit-01` checks the new residual-controlled
trust-region upper bound and reduced design. It found no validity error in
that proposition, while identifying seven numerical/design concerns. The
response documents fallback, selection bookkeeping, root evaluation, test and
guard changes made before the empirical joint fit began. This focused audit
does not rederive every inherited Fourier bound, has not reviewed its own
corrections, and supplies no acceptance verdict.
