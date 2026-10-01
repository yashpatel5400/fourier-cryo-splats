# Fixed-bank statistic separation and sampling error

1 October 2026 UTC. Post hoc reporting addition to the existing matched ledger. It uses every saved bank/level/deletion (96 rows), with no new simulated images, orientation draws or score tuning. A finite template-bank log ratio is a statistic even when it is an inaccurate Haar marginal likelihood. Its empirical squared separation D=(E[T1-T0])²/Var(T0) can be described without calling it a KL divergence or information ceiling.

| EMPIAR | Bank | Paired mean gap ± MC SE | Squared separation D ± delta MC SE |
|---|---|---:|---:|
| 10028 | 0 | 0.049656 ± 4.54e-05 | 0.0475194 ± 0.000764 |
| 10028 | 1 | 0.0496285 ± 4.67e-05 | 0.0471416 ± 0.000757 |
| 10049 | 0 | 0.0714897 ± 0.000318 | 0.0685766 ± 0.0012 |
| 10049 | 1 | 0.0716322 ± 0.000317 | 0.0681974 ± 0.00121 |
| 10076 | 0 | 0.00705424 ± 5.82e-05 | 0.00673921 ± 0.000146 |
| 10076 | 1 | 0.00704881 ± 5.83e-05 | 0.00674133 ± 0.000145 |

The table shows the pre-existing 25%-deletion, 32,768-template entries; all 96 entries remain in the CSV. SEs concern 8,192 independent simulated source images conditional on the fixed bank, maps, CTF and unit-white-noise/Haar model. The null and alternative of a source image share noise and pose. Bank errors are correlated and are not independent replications of a clinical or biological finding. No experimental model-error uncertainty is included.

For h_i=T1_i-T0_i, mean gap mu and sample null variance v, the first-order influence used for D is 2 mu/v (h_i-mu) - mu²/v² [n/(n-1)(T0_i-mean(T0))²-v]. Its sample standard deviation divided by sqrt(n) is the reported asymptotic delta SE. This accounts for estimated variance and covariance with the gap. It is not a finite-sample confidence bound or a normal-theory power guarantee. Three independent finite-difference contamination derivatives per row verify the implemented influence. The gap SE alone is not relabelled as an SE of D.

This reporting correction accepts D1 of the [focused consultation](reviews/post-round04-method-consultation/critique.md) only as a descriptive separation comparison. The existing template statistic shows substantially larger separation than the moments under the matched simulator. Its instability still prevents treating it as a converged marginal likelihood, a statistical impossibility bound or a likelihood certified uniformly over parameters. The improved continuous integration gate uses different image-specific proposals and does not retroactively change these frozen template-bank scores.
