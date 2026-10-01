# Completed phase-alignment independence control

1 October 2026 UTC. Twenty declared cases, 2,000 independently simulated datasets per case, five paired procedures, and all 40,000 sets of interval endpoints are retained. This is a one-frequency phase model with known Gaussian coordinate noise, not a cryo-EM reconstruction benchmark. Numerical tests check the product-normal MGF/moments and confidence-set inversion.

| Amplitude / noise | Particles | Oracle coverage | Same-image Student | Independent nominal | Independent Student | Invariant bound |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 16 | 0.9545 | 0.0005 | 0.7175 | 0.7750 | 0.9990 |
| 0.5 | 64 | 0.9585 | 0.0000 | 0.2230 | 0.2640 | 0.9975 |
| 0.5 | 256 | 0.9505 | 0.0000 | 0.0000 | 0.0000 | 0.9990 |
| 0.5 | 1024 | 0.9510 | 0.0000 | 0.0000 | 0.0000 | 0.9950 |
| 0.5 | 4096 | 0.9510 | 0.0000 | 0.0000 | 0.0000 | 0.9940 |
| 1 | 16 | 0.9510 | 0.2320 | 0.5655 | 0.6935 | 0.9980 |
| 1 | 64 | 0.9475 | 0.0000 | 0.0780 | 0.1285 | 0.9955 |
| 1 | 256 | 0.9475 | 0.0000 | 0.0000 | 0.0000 | 0.9975 |
| 1 | 1024 | 0.9485 | 0.0000 | 0.0000 | 0.0000 | 0.9960 |
| 1 | 4096 | 0.9475 | 0.0000 | 0.0000 | 0.0000 | 0.9970 |
| 2 | 16 | 0.9445 | 0.8095 | 0.7505 | 0.8405 | 0.9970 |
| 2 | 64 | 0.9450 | 0.3525 | 0.3185 | 0.4210 | 0.9970 |
| 2 | 256 | 0.9510 | 0.0025 | 0.0030 | 0.0045 | 0.9955 |
| 2 | 1024 | 0.9495 | 0.0000 | 0.0000 | 0.0000 | 0.9940 |
| 2 | 4096 | 0.9430 | 0.0000 | 0.0000 | 0.0000 | 0.9940 |
| 4 | 16 | 0.9530 | 0.9225 | 0.9175 | 0.9325 | 0.9955 |
| 4 | 64 | 0.9560 | 0.8290 | 0.8230 | 0.8335 | 0.9965 |
| 4 | 256 | 0.9605 | 0.4650 | 0.4325 | 0.4505 | 0.9935 |
| 4 | 1024 | 0.9410 | 0.0175 | 0.0160 | 0.0175 | 0.9935 |
| 4 | 4096 | 0.9510 | 0.0000 | 0.0000 | 0.0000 | 0.9945 |

At 4,096 particles, every naive aligned interval has zero observed coverage in all four SNR cases, including the independent-image intervals with sample variance. The oracle remains near .95 and the invariant finite-sample bound is conservative. These outcomes demonstrate the stated attenuation/noise-selection counterexample; they do not establish a new inverse method or a universal failure of independent data splitting.

The invariant bound uses a different, phase-invariant estimand before taking a square root, and requires the known equal-variance Gaussian model. It does not calibrate experimental noise, CTFs or physical heterogeneity. Empty sets and their counts are retained in the JSON; binomial intervals and widths are available for every case.

Protocol: [PHASE-SPLIT-CONTROL-PROTOCOL.md](PHASE-SPLIT-CONTROL-PROTOCOL.md). Proof: [PHASE-SPLIT-DERIVATION.md](PHASE-SPLIT-DERIVATION.md).
