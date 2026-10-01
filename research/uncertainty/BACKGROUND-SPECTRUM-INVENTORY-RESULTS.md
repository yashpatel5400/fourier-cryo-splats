# Descriptive background spectra of the three particle stacks

1 October 2026 UTC. All 8,192 downloaded particles per stack are included, with four 8×8 corner patches per particle. These are the already processed extracted particles, not new raw movie acquisitions. No particle or corner is removed. The [protocol](BACKGROUND-SPECTRUM-INVENTORY-PROTOCOL.md) was committed before outcomes.

The known Fourier-crop operation, including removal of the Nyquist row and column, is incorporated in a white-noise reference covariance. Patch means are removed and DCT coefficients transformed using that reference covariance. The real data are **not** empirically whitened. Consequently a white-noise reference would have equal expected normalized powers; the observed ratios below describe departures from that reference.

| EMPIAR | Patches | Scale | Normalized coordinate-power range | Second-moment eigenvalue range | Off-diagonal norm fraction | Median marginal kurtosis |
|---|---:|---:|---:|---:|---:|---:|
| 10028 | 32768 | 0.402278 | 0.366–2.101 | 0.358–2.360 | 0.224 | 3.133 |
| 10049 | 32768 | 0.127311 | 0.783–3.000 | 0.750–3.227 | 0.170 | 3.015 |
| 10076 | 32768 | 0.0688689 | 0.551–2.162 | 0.542–2.540 | 0.231 | 3.152 |

The nonflat patterns persist in both recorded halves. The largest coordinate-power ratios span approximately 3.8–5.7-fold within a stack after correcting the known preprocessing. These are descriptive departures, not p-values. The column labelled covariance eigenvalues uses the second-moment matrix; the centered covariance is also retained. Mean-vector energy is below 0.01% in each pooled half, so that distinction is small here. Scales are in the downloaded images’ intensity units and are not directly comparable to the unit-noise simulator without its whitening convention.

## Recorded-half stability

| EMPIAR | Cohort | Scale | Power range | Off-diagonal fraction |
|---|---|---:|---:|---:|
| 10028 | all | 0.402278 | 0.366–2.101 | 0.224 |
| 10028 | source_half0 | 0.401991 | 0.362–2.107 | 0.226 |
| 10028 | source_half1 | 0.402564 | 0.371–2.121 | 0.227 |
| 10049 | all | 0.127311 | 0.783–3.000 | 0.170 |
| 10049 | source_half0 | 0.127434 | 0.778–2.985 | 0.177 |
| 10049 | source_half1 | 0.127188 | 0.788–3.015 | 0.171 |
| 10076 | all | 0.0688689 | 0.551–2.162 | 0.231 |
| 10076 | source_half0 | 0.0690024 | 0.549–2.126 | 0.235 |
| 10076 | source_half1 | 0.0687364 | 0.550–2.197 | 0.233 |

## Separate corners

Corner order is top-left, top-right, bottom-left, bottom-right in array coordinates.

| EMPIAR | Corner | Scale | Power range |
|---|---|---:|---:|
| 10028 | 0 | 0.390239 | 0.381–2.064 |
| 10028 | 1 | 0.402875 | 0.361–2.127 |
| 10028 | 2 | 0.40153 | 0.367–2.085 |
| 10028 | 3 | 0.414468 | 0.357–2.132 |
| 10049 | 0 | 0.126444 | 0.782–2.928 |
| 10049 | 1 | 0.127215 | 0.780–3.043 |
| 10049 | 2 | 0.126558 | 0.775–3.040 |
| 10049 | 3 | 0.129026 | 0.767–3.095 |
| 10076 | 0 | 0.0668573 | 0.557–2.150 |
| 10076 | 1 | 0.068317 | 0.567–2.153 |
| 10076 | 2 | 0.0695377 | 0.539–2.191 |
| 10076 | 3 | 0.0707634 | 0.540–2.204 |

## Verification and interpretation

The explicit 63-frequency Fourier-basis covariance agrees with the analytic reference within 3.78e-15. Thirty-six selected patch transforms reproduce using an independent cosine matrix and linear solve within 1e-12. The first attempt stopped before reading particle outcomes because a 2e-15 reference-check tolerance was smaller than ordinary accumulation error. Its log and original source commit are retained; the second attempt uses 1e-13 and changes no data or procedure.

Background can contain particle tails, other particles, ice, and preprocessing artifacts. These patches are correlated within each particle. Existing metadata joins recover 229 / 137 / 351 source acquisition groups; the earlier inference that identities were unavailable was incorrect (see the [dated correction](BACKGROUND-SOURCE-GROUP-CORRECTION.md)). This analysis pooled particles and did not use those groups for inference. The recorded halves are not independent noise exposures. The analysis therefore cannot identify pure-noise covariance inside a molecule, establish stationarity, assign valid acquisition-level error bars, or certify the sub-percent noise precision demanded by the earlier moment scores. It does show that treating these archived image coordinates as exactly known white noise needs a measured and validated whitening model. A future likelihood study must retain this limitation instead of transferring unit-noise simulation guarantees directly to experiment.

Full second moments, centered covariances, marginal moments and source hashes are in `results/uncertainty/development/background-spectrum-inventory-v2/summary.json`; per-particle patch energies and powers are saved in the three NPZ files.
