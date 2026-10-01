# Final catalogue repair of orientation integration

1 October 2026 UTC. **All three stacks pass the unchanged numerical gate.** This completes the single permitted substantive repair. All 384 image cases, both banks, original optimizer failures and thirteen remaining out-of-tolerance images are retained. The cohort, maps, 8,192 draws per bank, tolerance and median-ESS thresholds were fixed before outcomes; the catalogue proposal uses continuous rotations and a full normalized mixture density.

| EMPIAR | Images within .01 | Four median ESS range | Four lower-decile ESS range | Stack gate | Minutes |
|---|---:|---:|---:|---|---:|
| 10028 | 128/128 (1.0000) | 2245.0–2254.7 | 2130.5–2138.4 | pass | 22.65 |
| 10049 | 124/128 (0.9688) | 2233.8–2262.4 | 1323.2–1426.5 | pass | 22.71 |
| 10076 | 119/128 (0.9297) | 1264.7–1346.8 | 512.0–549.8 | pass | 22.67 |

Each stack must have at least 90% of its images within .01 log-ratio units between banks and each of its four median ESS values at least 256. Lower-decile ESS is secondary: the later consultation suggested 900/450 on 10049/10076, but these were not substituted into the frozen gate. Timings are concurrent Mac wall times, not a GPU or matched-baseline speedup. No external compute was rented.

| EMPIAR | RMS ratio difference | Maximum ratio difference | RMS absolute log-integral difference | RMS prefix ratio change |
|---|---:|---:|---:|---:|
| 10028 | 0.000886931 | 0.00441846 | 0.0287539 | 0.000714897 |
| 10049 | 0.00451915 | 0.0175564 | 0.0295617 | 0.00455028 |
| 10076 | 0.00507446 | 0.0176302 | 0.0424311 | 0.00390426 |

Two proposals can miss the same mass; ratio error also benefits from common-mode cancellation between the two nearby maps. Passing this per-image gate does not certify all parameter values, a product likelihood over a large dataset, or real-image marginal likelihoods. It does not validate Haar experimental views, white noise, known amplitudes, correct maps or uncertainty coverage. The [population materiality screen](POPULATION-MATERIALITY-RESULTS.md) fails independently, so this numerical pass does not authorize the proposed scientific method.

## Independent checks

- 10028: all 128 integral/ESS/gate summaries and observations replay. Independent full-4D covariance densities check 8,192 points (32 fixed points in each bank per image), maximum log-density error 1.35e-13. Six direct 64³-cell residual calculations cover both models at all 220 frequencies, maximum discrepancy 6.93e-11.
- 10049: all 128 integral/ESS/gate summaries and observations replay. Independent full-4D covariance densities check 8,192 points (32 fixed points in each bank per image), maximum log-density error 1.32e-13. Six direct 64³-cell residual calculations cover both models at all 220 frequencies, maximum discrepancy 4.64e-11.
- 10076: all 128 integral/ESS/gate summaries and observations replay. Independent full-4D covariance densities check 8,192 points (32 fixed points in each bank per image), maximum log-density error 1.56e-13. Six direct 64³-cell residual calculations cover both models at all 220 frequencies, maximum discrepancy 2.4e-11.

The replay recomputes every stored importance integral and verifies identical observations to attempt one. The density check uses explicitly transformed inverse 4D covariance matrices and the physical check explicitly sums cells; neither calls the runner’s corresponding shortcut. Only the stated density/physical subset is checked independently, not every integrand point.

## Post hoc saved-draw autopsy

| EMPIAR | Proposal | Original .01 criterion | Images | Mean matching-state Haar weight mass | Largest planted log-integral jump | Median bank discrepancy / delta SE |
|---|---|---|---:|---:|---:|---:|
| 10028 | adaptive | within | 128 | 3.83422e-06 | 0.000475271 | 0.6073 |
| 10028 | catalog | within | 128 | 6.79125e-06 | 0.000951043 | 0.6714 |
| 10049 | adaptive | within | 103 | 0.00639017 | 0.0602882 | 0.7452 |
| 10049 | adaptive | outside | 25 | 0.0983724 | 0.529422 | 1.274 |
| 10049 | catalog | within | 124 | 0.0012939 | 0.00547181 | 0.5873 |
| 10049 | catalog | outside | 4 | 0.0030117 | 0.00106534 | 1.6 |
| 10076 | adaptive | within | 79 | 0.0803211 | 1.49593 | 0.4681 |
| 10076 | adaptive | outside | 49 | 0.149842 | 1.79773 | 1.283 |
| 10076 | catalog | within | 119 | 0.00691675 | 0.062589 | 0.5899 |
| 10076 | catalog | outside | 9 | 0.0113052 | 0.002712 | 1.415 |

The planted diagnostic replaces draw zero with the saved generating pose under the matching state only. Its [reciprocal-unbiased identity and assumptions](PLANTED-TRUTH-INTEGRATION-AUDIT.md) are proved by exchangeability and checked by exact finite-space enumeration. Large jumps expose sensitivity to posterior mass absent from ordinary saved draws; they do not reveal the exact integral. In particular, the old 10076 within-tolerance group contains a jump of 1.50 log units, illustrating why bank agreement is insufficient. After repair its largest within-tolerance jump is .0626. No original mode has a clipped Hessian eigenvalue, so clipping does not explain these failures. Delta SEs remain tail-sensitive plug-in diagnostics.

| EMPIAR | Proposal | Four 64-image stochastic bracket widths |
|---|---|---:|
| 10028 | adaptive | 7.38248–7.38566 |
| 10028 | catalog | 7.39543–7.3981 |
| 10049 | adaptive | 7.56229–8.06717 |
| 10049 | catalog | 7.40416–7.41509 |
| 10076 | adaptive | 9.84761–10.2598 |
| 10076 | catalog | 7.50787–7.51896 |

These are pointwise 95% stochastic brackets from one forward and one planted reciprocal estimate per independent image, with .025 allocated to each tail. Each bracket concerns one generating state and one bank. The paired states are never multiplied together, and no simultaneous coverage across the displayed brackets is asserted. The additive floor 2 log(40) is about 7.378 log units. These intervals are simulation-only and too coarse to act as a general numerical certificate.

## Decision and chronology

The numerical branch ends here as planned: no further proposal, images, draws or favorable subset are selected. The method consultation was complete at the provider at 11:38:56 UTC but had not been retrieved or read when the 11:40 scheduling amendment and run were made. Its later recommendation of a local-majority proposal differs from the frozen .45 local/.50 catalogue/.05 Haar choice. The [response correction](reviews/post-round04-method-consultation/response.md) preserves that distinction. Original protocols and provider text remain unchanged.

All numerical points, labels, catalogue weights, physical kernels, source hashes and unfavorable cases are preserved. This is a classical integration repair and a development audit, not a new uncertainty theorem, an experimental population estimate or a fifth full-paper review.
