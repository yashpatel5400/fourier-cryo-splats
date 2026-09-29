# Independent review record

No scientific review has yet been performed. The exact requested model,
`claude-fable-5-1`, succeeded in an availability-only invocation through the
authenticated Claude CLI. An availability response is not a research review.

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
