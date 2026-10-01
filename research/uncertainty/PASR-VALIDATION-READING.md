# Acquisition processing and perturbation validation

1 October 2026 UTC. One additional candidate; targeted reading, not baseline execution.

[Burton Smith and Murata (2026)](https://journals.iucr.org/m/issues/2026/05/00/rq5016/) study pixel doubling before movie alignment for datasets already limited by sampling. Their validation combines reconstruction comparisons, sub-pixel flipping, map-to-model FSC and Q-scores. The discussion distinguishes processing movies from already aligned micrographs and notes unresolved noise-model questions. These controls address acquisition/reconstruction fidelity, not confidence coverage for a fixed density functional.

Our inference from this reading is that a validation study must retain its processing order and specify what a perturbation is intended to detect. Successful processing or plausible FSC does not establish the conditional noise distribution used by our theorem. We do not apply this preprocessing to the frozen movie pilot.

Read the primary abstract, method overview and Figure 2, perturbation-control discussion, model-validation and data-availability sections. Supplementary figures and author code were not inspected. Source/access hashes are in `pasr-reading-source.json`. CryoDiff full text remained blocked on this access attempt; no new calibration claim is attributed to it.
