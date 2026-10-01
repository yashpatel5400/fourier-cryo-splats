# Population materiality screen: result

1 October 2026 UTC. **The primary gate fails on all three stacks.** All thirty prelisted pair/assignment/grid comparisons complete. The frozen criterion requires at least two stacks to have modeled absolute population bias at least .02 and harmonic amplitude-tangent surrogate at most 8.8125, for the same pair on the primary 8-bin-per-axis grid. No sensitivity grid substitutes for a primary outcome.

| EMPIAR | Largest absolute primary bias | Primary harmonic surrogate range | Stack gate |
|---|---:|---:|---|
| 10028 | 0.00021411 | 1.290924–1.291933 | fail |
| 10049 | 0.00012183 | 1.038862–1.039476 | fail |
| 10076 | 0.00349347 | 3.871904–3.916364 | fail |

The failure is materiality, not the harmonic threshold. The largest primary bias, .003493 on 10076, is about 0.35 percentage points versus the required 2 percentage points. No sensitivity-grid bias reaches .02 either. Under the stopping rule, we do not build the proposed population-feature method on these compact regions. This does not show that pooled likelihood is generally adequate, that latent poses are harmless, or that other conformations have negligible population bias.

## What was calculated

A is each previously fitted normalized map; B deletes its fixed 20-Angstrom peak region. The population is .75. The calculation knows the pose, one CTF and white noise, and uses 220 Fourier frequencies. We compare the two published source halves on all stacks and published filtering versus its complement on 10028/10076, in both assignments to A/B. These metadata groups are not biological states. Their estimated rotations define a piecewise-Haar proxy through fixed 4/8/12-bin ZYZ grids. No actual biological class-specific law or experimental population truth is available in these files.

The exact-in-model population-score root uses one-dimensional Gaussian integration; its odds relation is implicit. The weak-signal approximation and every signed outcome are retained in the CSV. The training draws are reused development data, not an independent confirmation set. Within-cell Monte Carlo errors in each JSON condition on the fixed histograms and omit pose/model/metadata uncertainty.

## Amplitude separation diagnostic

| EMPIAR | Haar mean region energy | Full-map tangent residual | Arithmetic reduction | Haar inverse-energy mean | Haar projected harmonic surrogate |
|---|---:|---:|---:|---:|---:|
| 10028 | 0.819637 | 0.779043 | 4.95% | 1.225737 | 1.291212 |
| 10049 | 1.494523 | 0.939452 | 37.14% | 0.714794 | 1.120616 |
| 10076 | 0.280687 | 0.257457 | 8.28% | 3.603656 | 3.918368 |

Real-amplitude projection removes about 5%, 37% and 8% of arithmetic mean separation, respectively. This is a local Gaussian-mean diagnostic. Its harmonic substitution is not the variance theorem for a population estimator with unknown amplitudes; that theorem in the algebra note requires amplitude one and known poses. A small projected harmonic value neither certifies latent-pose identifiability nor provides an uncertainty interval.

## Verification

- 10028: all 65,536 information/projection values and every metadata histogram replay; 12 expected population scores checked by direct density-ratio integration with a separate 160-node quadrature. Maximum score residual 1.84e-11; maximum amplitude-residual discrepancy 1.33e-15; twelve direct cell sums across the first, boundary and last training batches differ by at most 1.08e-11.
- 10049: all 65,536 information/projection values and every metadata histogram replay; 6 expected population scores checked by direct density-ratio integration with a separate 160-node quadrature. Maximum score residual 3.07e-11; maximum amplitude-residual discrepancy 2.66e-15; twelve direct cell sums across the first, boundary and last training batches differ by at most 5.59e-12.
- 10076: all 65,536 information/projection values and every metadata histogram replay; 12 expected population scores checked by direct density-ratio integration with a separate 160-node quadrature. Maximum score residual 3.36e-12; maximum amplitude-residual discrepancy 5e-16; twelve direct cell sums across the first, boundary and last training batches differ by at most 2.85e-12.

Four targeted tests compare the formulas with direct density integration, separately optimized expected likelihood and explicit real-coordinate least squares. The later-batch projection checks also address the external reviewer’s concern about orientation regeneration beyond its first sample. No new images, regions, seeds, extra rotations or relaxed thresholds were introduced.

See the [frozen protocol](POPULATION-MATERIALITY-PROTOCOL.md), [algebra audit](POPULATION-SCREEN-ALGEBRA.md) and [unchanged focused consultation](reviews/post-round04-method-consultation/critique.md). This negative development result supplies neither a new uncertainty method nor an ICML acceptance claim.
