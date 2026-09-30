# Development release v0.6.0-dev

This is an active research preprint, not an acceptance-ready submission. The
first full Claude Fable 5.1 review remains reject (confidence 4/5). The later
focused mathematical audit is preserved and does not replace a full review.

This snapshot adds all 24 original-code CryoLike score comparisons, the completed
continuous-orientation and cubic-design follow-ups (including unsuccessful
ones), and the first converged unknown-pose RELION continuation on EMPIAR-10028.
It also retains the separately declared pilot-only reference registration on
all three stacks and both registered/unregistered reference comparisons.
Other RELION continuations remain in progress; they have no published result
in this snapshot. Continuous likelihood bounds remain too loose for a useful
test, and experimental uncertainty assumptions remain unresolved.

The 50-page ICML-format development paper includes the new outcomes, definitions,
proofs, all reporting limitations, and the 116-candidate critical survey.
The archive contains 905 arrays, including 686 explicitly indexed post-review
owners. Every archive member and its byte hash were independently streamed and
verified. Its compressed size is 1,219,119,797 bytes; SHA256:
`bd1a2f02c59f6d9c1c6e87a9a202ae1d9cfa9c4a3f69ebae4758fd7f041bd0d1`.

Extract the array archive at repository root. The manifest identifies exact
paths and owner records; the source repository retains provenance, per-case
statuses, protocols, FSC curves and review reports. The original trained
checkpoints still require release v0.2.0-dev and original-data downloads use
the recorded acquisition protocols. No raw experimental particle pixels or
third-party paper PDFs are bundled. Earlier releases are unchanged.
