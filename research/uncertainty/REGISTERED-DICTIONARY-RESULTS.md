# Registered replay of all original dictionary examples

1 October 2026 UTC. All twelve estimators exactly reproduce the archived original widths, biases and analytic coverages (maximum stored replay discrepancy zero). Original-versus-two-stage reference resampling differs by at most 2.15e-7 relatively. Registration is the previously frozen pilot-only transform; no estimator is refitted to improve coverage. These are conditional finite-voxel simulations, not experimental truth labels.

| Stack | Target | Width / field | Original coverage | Registered coverage | Registered audited coverage |
| --- | --- | ---: | ---: | ---: | ---: |
| 10028 | center | 0.03 | 0 | 0 | 1 |
| 10028 | z_contrast | 0.03 | 0.916352 | 0.054327 | 1 |
| 10028 | center | 0.07 | 4.18037e-05 | 3.73035e-14 | 1 |
| 10028 | z_contrast | 0.07 | 0.79958 | 0.367381 | 1 |
| 10049 | center | 0.03 | 0.814392 | 0.75855 | 1 |
| 10049 | z_contrast | 0.03 | 0.00176813 | 0.506725 | 1 |
| 10049 | center | 0.07 | 0.87122 | 0.927568 | 1 |
| 10049 | z_contrast | 0.07 | 4.89913e-16 | 0.785984 | 1 |
| 10076 | center | 0.03 | 0.00129502 | 0.249517 | 1 |
| 10076 | z_contrast | 0.03 | 0.605057 | 0.379971 | 1 |
| 10076 | center | 0.07 | 0.917992 | 0.942364 | 1 |
| 10076 | z_contrast | 0.07 | 0.669595 | 0.898948 | 1 |

The highlighted 10028 central .07 example changes from 4.18037e-5 to 3.73035e-14 coverage. Other targets improve or worsen; none are omitted. The ambient-audited intervals cover these two specific generators with probability numerically one, which is conservative and does not establish that an experimental density belongs to the class. The 10076 generator remains one Class A map, not heterogeneous consensus truth.

Declared protocol: [REGISTERED-DICTIONARY-PROTOCOL.md](REGISTERED-DICTIONARY-PROTOCOL.md). Saved weights and both reference generators accompany the three geometry records.
