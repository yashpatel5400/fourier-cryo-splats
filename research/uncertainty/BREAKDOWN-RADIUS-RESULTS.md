# Density and noise breakdowns of the existing experimental intervals

1 October 2026 UTC. All 96 fixed-weight combinations replay the archived B=2 raw half-widths within the declared tolerance. The calculation varies the assumed density radius while holding the observed center, pilot, poses, image weights and noise calibration fixed. It does not choose a new radius or validate the physical class.

The table gives the centered-noise procedure at fixed poses. B* is the boundary where the raw data-based interval first includes zero as B increases. A zero B* means it already includes zero at B=0. The original declared radius remains B=2. Width/center uses the raw interval, avoiding the different no-data center.

| Stack | Feature | Half-width / absolute center | B* | Critical SD / reported SD |
| --- | --- | ---: | ---: | ---: |
| 10028 | pilot_region_1 | 0.655 | 4.832 | 1.836 |
| 10028 | pilot_region_2 | 0.649 | 4.868 | 1.867 |
| 10028 | pilot_region_3 | 0.861 | 2.877 | 1.257 |
| 10028 | matched_center | 0.832 | 3.055 | 1.326 |
| 10049 | pilot_region_1 | 0.674 | 153.172 | 1.483 |
| 10049 | pilot_region_2 | 1.501 | 0.000 | 0.666 |
| 10049 | pilot_region_3 | 2.502 | 0.000 | 0.399 |
| 10049 | matched_center | 0.711 | 133.654 | 1.406 |
| 10076 | pilot_region_1 | 0.256 | 186.381 | 3.921 |
| 10076 | pilot_region_2 | 0.347 | 126.048 | 2.889 |
| 10076 | pilot_region_3 | 0.267 | 178.140 | 3.763 |
| 10076 | matched_center | 0.281 | 161.521 | 3.569 |

At the one-degree pose class, every feature has a zero noise-SD breakdown: its recorded bias bound already exceeds the observed center, even with measurement SD zero. At B=0, neither homogeneous stack has a one-degree exclusion under either calibration; centered 10076 regions 1 and 3 have only tiny density-radius thresholds (.054 and .016). This diagnoses the existing fixed-weight audit, not intrinsic impossibility or a physical error radius. The separately optimized cubic estimators are not part of these twelve frozen estimators.

The complete CSV retains raw and centered calibration, fixed/shift/one-/two-degree classes, selected and raw widths, fallbacks, replay errors and thresholds. Large B* values can reflect a small fixed-pose adjoint residual; they do not validate unknown-pose inference. The five centered 10076 registered-reference disagreements remain unchanged.

Declared protocol: [BREAKDOWN-RADIUS-PROTOCOL.md](BREAKDOWN-RADIUS-PROTOCOL.md).
