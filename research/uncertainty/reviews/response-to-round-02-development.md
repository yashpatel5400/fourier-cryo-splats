# Response to full independent review 2 — active development

1 October 2026 UTC. The unaltered review recommends rejection, confidence 4/5. We accept that the current work is not an ICML contender. No new acceptance assessment has occurred. This is a work plan, not a claim that the concerns have been resolved.

## Priorities and acceptance checks

1. **R10/R7: operator-matched baseline.** Implement the continuous isonormal Gaussian-prior calculation with prior directional SD B and B/2, retaining fixed-pose and pilot first-order pose-marginal variants. Run all twelve locked targets. Distinguish a generalized white prior from an L2-valued density prior; verify and explain the frequentist bias/variance relation rather than assuming all credible intervals lack a uniform guarantee.
2. **R9: mixed pose model.** Derive a conditional bound for independent centered local pose errors and a declared common-mode component. Treat this as an explicit alternative assumption, not a consequence of a bounded pose radius. First-order bounded variables require a Hoeffding/sub-Gaussian critical value, not an unjustified Gaussian quantile. Estimated designs/weights may violate the assumed independence.
3. **R11/A: end-to-end validation.** Build a simulation with a fixed local pose refinement and recomputed weights on every replicate; use a separate tuning batch, then freeze sources and protocol before 200 replicate datasets per each of three geometries. Report failures, empirical coverage with Monte Carlo intervals, widths and pose-error distributions. Include same-image and independent alignment/inference noise controls where computationally feasible. This is a local-refinement simulation, not ab initio or experimental coverage.
4. **R1/D/R2/R8:** one homogeneous experimental stack needs evidence-based nuisance inputs and useful inference. Raw-frame independence and covariance transfer are still unresolved. No relaxation of density/noise/pose assumptions will be presented as empirically verified. Existing 10076 work remains heterogeneous stress evidence.
5. **R12/M:** rewrite the paper around these substantive outcomes. Report prior observed centres, registered-reference disagreements, relative-to-estimate widths and breakdown radii. Recompute the excluded-dictionary example in the registered frame, correct protocol and survey counts, incorporate resolving-kernel/GP prior art with primary-source checks. Move historical development out of the main paper.

## Preserved facts

Both full reviews reject. The positive mathematical checks do not resolve usefulness or calibration. All five cubic designs' experimental raw/centered intervals contain zero. Existing analytic coverage tables remain conditional implementation checks. No new experiment or new theorem is claimed in this response yet. Frozen calibration files and historical source snapshots will not be modified; revised evaluators will save failed outcomes before signalling an error.

## Review access limits

Provider identity is exact Claude Fable 5.1; evidence immutability verified. The transport records all 55 uncropped pages, while the reviewer reports seeing only pages 20–55. It did not run tests, inspect all individual frozen records, or read third-party sources. This does not change or invalidate its substantive critique.
