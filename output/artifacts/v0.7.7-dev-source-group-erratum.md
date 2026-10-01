# Correction: source acquisition groups are available

1 October 2026 UTC. The first background-spectrum protocol and reports wrongly
inferred incomplete micrograph identities from the numbers of `.cs` blob paths
(229 / 2 / 1). Those are stack-storage paths. The project had already recovered
source groups from deposited STAR/Frealign metadata in `prepare_uq_splits.py`.
This was an interpretation error in the new narrative, not a new data discovery.

| EMPIAR | Selected particles | Source groups | Source field | Median particles/group |
|---|---:|---:|---|---:|
| 10028 | 8,192 | 229 | STAR rlnMicrographName | 36 |
| 10049 | 8,192 | 137 | STAR rlnMicrographName | 64 |
| 10076 | 8,192 | 351 | Frealign FILM, column 8 | 21 |

For 10028/10049, the join uses original stack path and image index. For 10076,
the Frealign particle IDs and all defocus U/V values match the deposited STAR;
the STAR's sequential micrograph numbering is not used as the acquisition ID.
The complete source tables have 1,081 / 349 / 2,071 groups respectively.
`scripts/verify_background_source_groups.py` replays every selected join,
checks the original split hashes, and confirms that no recovered group crosses
the saved pilot/tune/inference/test splits. The audit summary is retained in
`results/uncertainty/development/background-source-group-correction-v1/summary.json`.

No spectra, optimizer outcomes, simulation observations, splits or numerical
gates are changed. Acquisition grouping permits cluster-level descriptive
checks and resampling under additional sampling assumptions. It does not prove
independence between groups, pure noise in corners, stationarity inside a
particle, independence of published consensus poses, or sub-percent noise
calibration. The original background analysis did not use these groups and
cannot retrospectively be described as acquisition-level inference.

The four v0.7.7-dev numerical archives were frozen before this correction. Their
bytes and manifests are preserved. Their embedded background protocol/results
and current-status narrative must be read together with this erratum; the live
repository documents link it. The focused post-round-4 Fable consultation was
also sent the uncorrected narrative. Its immutable prompt is retained, and a
factual correction must accompany interpretation of its response. This is not
a new full manuscript review or a change to any existing rejection verdict.
