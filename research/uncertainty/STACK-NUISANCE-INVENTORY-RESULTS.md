# Recorded imaging nuisances on the three stacks

1 October 2026 UTC. All three declared metadata inventories finish. This responds to review 4 R18/R19 at the descriptive level; it does not establish the true viewing law, physical amplitude bounds or noise spectra. Source poses are estimates and can reflect symmetry conventions, shared processing and particle selection. The [frozen protocol](STACK-NUISANCE-INVENTORY-PROTOCOL.md) states all grids and subsets.

## Recorded viewing directions

For 128 equal-area bins of unoriented plane normals, the peak ratios and cross-half chi-square statistics are well above those compatible with a recorded distribution near uniform. Both are resolution-dependent summaries. A binned maximum is not a continuous-density upper bound. Half-set agreement does not certify pose accuracy or independent processing.

| Stack/cohort | Particles | Peak ratio | Cross-half chi-square | Half-set TV |
|---|---:|---:|---:|---:|
| 10028 / all_metadata | 105,247 | 4.5327 | 1.0738 | 0.0252 |
| 10028 / published_filter | 93,852 | 4.5839 | 1.0910 | 0.0246 |
| 10049 / all_metadata | 108,544 | 6.2300 | 2.0067 | 0.0207 |
| 10076 / all_metadata | 131,899 | 2.7396 | 0.2478 | 0.0254 |
| 10076 / published_filter | 87,327 | 3.3859 | 0.4429 | 0.0299 |

All 32/128/288 normal-bin and 64/512/1,728 full-rotation-bin counts, including both recorded halves, are retained in the JSON. On all-metadata cohorts, normal-bin peak ratios range from 2.62 to 6.07 (10028), 4.87 to 7.33 (10049), and 2.06 to 3.13 (10076) as bin count increases. The uniform finite-count expectation for plug-in chi-square is recorded separately; the observed deviations are much larger. The published filter changes 10076 appreciably, illustrating that source selection matters.

## CTFs and recorded scale fields

| Stack, all metadata | Defocus mean, 5th–95th percentile (micrometres) | Alpha, 5th/50th/95th percentile | Alpha variation associated with view bins |
|---|---|---|---:|
| 10028 | 1.349–2.783 | 0.899 / 1.204 / 1.394 | 2.62% |
| 10049 | 1.415–2.192 | 1.000 / 1.000 / 1.000 | constant recorded field |
| 10076 | 1.053–2.438 | 0.590 / 1.154 / 1.559 | 19.44% |

All CTF scale entries equal one and all alignment weight entries equal zero. Those constants are metadata values, not physical confidence statements. Alpha varies on 10028 and 10076, but its exact legacy processing interpretation has not been independently calibrated. The 10076 view association falls from 19.44% to 0.74% after the published filter. No causal claim or valid amplitude interval follows. The constant alpha on 10049 cannot establish absence of physical amplitude variation.

The nontrivial per-particle defocus spread is incompatible with interpreting the earlier one-CTF simulator as an experimental forward model for the full cohort. These files contain scalar residual/image powers, not frequency-resolved pure-noise covariance estimates. They therefore cannot justify the exactly known white-noise assumption. Existing image-based calibration still has the previously reported signal-contamination and dependence limits.

## Verification and consequence

An independent direct-atan2 / `histogramdd` implementation reproduces all 90 full/half histogram count vectors exactly, and verifies every input hash and all scalar finite counts/extrema. The main runner checks every stored rotation against the transposed cryoSPARC rotation-vector convention and checks orthogonality. Neither verifier supplies outward-rounded numerical certification, all-quantile replay, or a statistical guarantee for latent poses.

The earlier kappa=1.1 Haar comparison is a controlled simulation assumption, not a nuisance level supported by these metadata. This diagnosis supports freezing that branch and comparing latent-pose likelihood information before selecting a replacement. It does not prove that cryo-EM uncertainty inference is impossible or that a more flexible viewing model will work.
