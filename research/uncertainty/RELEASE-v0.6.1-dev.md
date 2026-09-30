# Development release v0.6.1-dev

Active research, not an acceptance-ready submission. The first authentic full
Claude Fable 5.1 review remains reject; no later full acceptance assessment is
claimed. The cubic enrichment fit and two RELION continuations remain ongoing
and contribute no completed result to this snapshot.

This incremental release preserves the registered-target sensitivity, all three
experimental cubic applications, and the separately declared centered-noise
alternative. Correct registration raises one unchanged cubic estimator's
simulated minimum sign power from .00654 to .95725, but both its uncentered
and centered experimental intervals contain zero. Centering changes the
fixed/shift-only zero exclusions from six/four to ten/six; rotational exclusions
remain zero. Five narrower intervals disagree with the Class A reference on
heterogeneous 10076. These comparisons are not ground-truth coverage labels.

The corrected 52-page ICML-format paper includes these outcomes and a
117-candidate critical survey. It corrects the stale v0.6 main-text statement
that no converged ab initio comparison had run: the completed 10028 RELION
result was already reported in that release's appendix. Earlier release files
remain immutable.

The seven-array incremental archive is 46,220,151 bytes, SHA256
`de423a17ec36be407d05348bcaabce539c85e3611b1c2d7ed2cadf06067f7054`. Every member has been independently streamed and
verified. Extract at repository root after the v0.6.0-dev array archive.
Unlike a weights-only archive, this increment includes saved transformed
Fourier calibration and inference observations derived from public EMPIAR
particles. Original full stacks and third-party paper PDFs are not included.
Source records identify acquisition, array and exact implementation hashes.

The complete suite before the later enrichment module reports 230 passed,
1 skipped and 1 intentional warning; enrichment has two separate passing
prerequisite tests. Neither count establishes experimental calibration.
