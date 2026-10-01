# Classical risk comparators and Monte Carlo slack: complete results

Post-outcome analysis of the same two calibration archives, motivated by focused audit V1--V5. All eight score/candidate cases are retained. The grouped methods require amplitude conditionally independent of noise given view. The individual envelope allows a larger null. None of the four new procedures is selected by taking an unadjusted minimum with another.

## Removed-candidate probability bounds at kappa=1.1

| Stack | Score | Paired variance | CVaR DKW | CVaR split EB | Mean only | Unpaired variance |
|---|---|---:|---:|---:|---:|---:|
| 10028 | power_range09_11 | 0.486165 | 0.495290 | 0.491135 | 0.526781 | 0.499463 |
| 10028 | power_bispectrum_range09_11 | 0.478401 | 0.487299 | 0.483099 | 0.517880 | 0.491474 |
| 10049 | power_range09_11 | 0.430953 | 0.441040 | 0.436137 | 0.466863 | 0.444721 |
| 10049 | power_bispectrum_range09_11 | 0.447611 | 0.457221 | 0.452201 | 0.484699 | 0.461019 |

## Variance-bound decomposition (all candidates)

| Stack | Score | Candidate | Raw product mean | Square-root term | Range term | Distance subtraction | Final upper V |
|---|---|---|---:|---:|---:|---:|---:|
| 10028 | power_range09_11 | true | 3.19114e-05 | 0.000175438 | 0.000299959 | 3.33044e-08 | 0.000507275 |
| 10028 | power_range09_11 | removed | 3.16225e-05 | 0.000174808 | 0.000307394 | 2.90753e-07 | 0.000513534 |
| 10028 | power_bispectrum_range09_11 | true | 2.86195e-05 | 0.000177149 | 0.000300446 | 2.46495e-06 | 0.00050375 |
| 10028 | power_bispectrum_range09_11 | removed | 7.73654e-05 | 0.000176711 | 0.000311333 | 3.85246e-06 | 0.000561557 |
| 10049 | power_range09_11 | true | -8.60264e-05 | 0.000174302 | 0.000323941 | 0 | 0.000412216 |
| 10049 | power_range09_11 | removed | -9.84214e-05 | 0.000171225 | 0.000340086 | 0 | 0.000412889 |
| 10049 | power_bispectrum_range09_11 | true | 1.73772e-05 | 0.000173294 | 0.000312459 | 0 | 0.00050313 |
| 10049 | power_bispectrum_range09_11 | removed | -3.22728e-05 | 0.000172574 | 0.00033156 | 0 | 0.000471861 |

All 480 new projections and paired critical-count differences are in CSV. For removed 10049 candidates, raw centered-product estimates are negative; adding the finite Monte Carlo concentration terms makes the upper bounds positive. These results confirm slack-dominated variance bounds and do not measure negligible population view variance.

The paired method gives smaller bounds than these implemented CVaR and unpaired comparators at kappa=1.1, for this particular budget. This does not refute exact population CVaR optimality: CVaR is applied to noisy view-group means, whereas the paired product targets their conditional mean variance. Nor does it prove superiority under optimized allocation, other replication budgets, or independent calibration repetitions. Those comparisons remain unperformed. At fixed saved frequencies, method ordering follows the critical counts directly; overlapping single-image probability intervals do not test a paired method difference.

Six targeted tests check the CVaR density linear program (including fractional boundary atoms), explicit DKW breakpoint evaluation, conditional Jensen ordering, a noise-adaptive-amplitude counterexample, near-degenerate cubics and the variance decomposition. They do not prove scientific novelty. The primary reading scope and all formulas are in the frozen protocol.
