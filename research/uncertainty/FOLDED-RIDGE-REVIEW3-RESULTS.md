# Direct folded-width comparison after review 3

All six declared targets completed in 45.57 seconds. All 119 CG solves converged. The selected widths are .062–.098% smaller than evaluating the folded critical value on the previous sum-objective weights. This small improvement does not resolve pose estimation, density-class calibration, or novelty concerns.

| Stack | Target | Old width | Selected width | Width ratio | Restricted gap | Global gap | Solves |
|---|---|---:|---:|---:|---:|---:|---:|
| 10028 | center | 1.408308 | 1.407218 | 0.999226 | 0.004751 | 0.238005 | 20 |
| 10028 | contrast | 2.068084 | 2.066314 | 0.999144 | 0.004850 | 0.289202 | 20 |
| 10049 | center | 2.006489 | 2.004840 | 0.999178 | 0.004920 | 0.004920 | 20 |
| 10049 | contrast | 2.788483 | 2.785761 | 0.999024 | 0.004687 | 0.033191 | 24 |
| 10076 | center | 1.723184 | 1.722114 | 0.999379 | 0.004582 | 0.004582 | 17 |
| 10076 | contrast | 2.466517 | 2.464952 | 0.999365 | 0.004491 | 0.057139 | 18 |

The search covers [initial ridge / 64, initial ridge * 64]. Its monotonicity-based relative gap is .00449–.00492 on that finite range. The global outer-ray diagnostic remains .238 and .289 on 10028, .00492 and .0332 on 10049, and .00458 and .0571 on 10076. Restricted convergence is not global optimality. Root evaluation, quadrature and CG errors are accounted for as described in the protocol; floating-point/NUFFT pads are diagnostic rather than validated enclosures.

The known-pose observation design, targets, B=2 and alpha=.05 are unchanged. No inference images choose the regularization. The source record retains every ridge evaluation and both matched Gaussian prior-scale widths. Frozen refitting results are unchanged.
