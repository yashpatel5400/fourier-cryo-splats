# Matched information and nuisance diagnosis

1 October 2026 UTC. All three frozen stage-A/stage-B comparisons complete. These are post hoc diagnostics on the same candidate, 20 Å region, 220 complex frequencies, transfer profile and 8,192 paired images per stack. The [protocol](MATCHED-INFORMATION-LEDGER-PROTOCOL.md) preceded their outcomes. The original simulator still assumes Haar orientations, one CTF, amplitude one and known independent unit real/imaginary noise. No experimental calibration or new test is claimed.

## Matched separation at 25% deletion

Define D(T) = (E_alt T − E_null T)² / Var_null T. Moment means/conditional variances use exact Gaussian formulas, with Monte Carlo integration over views. The known-pose column is twice Gaussian KL averaged over these views; it also equals D for a conditionally centered known-pose linear score. The numerical likelihood columns are D of approximate Haar mixture log ratios, not certified KL values or information ceilings. Sampling standard errors for the paired likelihood mean gap and known-pose energy are retained in the full JSON; they do not include quadrature bias.

| EMPIAR | Known-pose 2 KL | Power matched / Fisher D | Combined matched / Fisher D | Numerical likelihood D, banks 0 / 1 |
|---|---:|---:|---:|---:|
| 10028 | 0.0511916 | 0.00133703 / 0.00137886 | 0.00102825 / 0.00161733 | 0.0475194 / 0.0471416 |
| 10049 | 0.0934722 | 0.0015915 / 0.00160227 | 0.00173507 / 0.0017262 | 0.0685766 / 0.0681974 |
| 10076 | 0.0175716 | 0.000118677 / 0.000118666 | 0.000122378 / 0.000120901 | 0.00673921 / 0.00674133 |

## Integration is not converged

| EMPIAR | Null ESS median, bank 0 / 1 | Bank RMS / score SD, null | Last doubling RMS, bank 0 / 1 |
|---|---:|---:|---:|
| 10028 | 1.23 / 1.226 | 0.07166 / 0.22779 | 0.061573 / 0.061523 |
| 10049 | 7.752 / 7.674 | 0.065906 / 0.273 | 0.051568 / 0.05232 |
| 10076 | 28.72 / 28.76 | 0.023613 / 0.08593 | 0.017105 / 0.017212 |

The mean separation agrees reasonably between banks, while individual log ratios and the last quadrature doubling remain materially unstable. In particular, typical 10028 images have nearly single-orientation support. A shared ratio can partly cancel density errors, but that does not certify accurate mixture densities or calibrated inference. More uniform random orientations alone is not the next statistical-method variant; the unresolved numerical integration must be addressed directly. These diagnostics justify investigating a likelihood approach, not reporting its asymptotic efficiency or treating failure as impossibility.

## Amplitude and noise sensitivity

| EMPIAR | Frozen score | Sum of power weights | White-noise variance error matching deletion gap | View share of null variance | Alternative outside amplitude envelope: global / view-dependent |
|---|---|---:|---:|---:|---|
| 10028 | power_matched | 10.9504 | 0.5136% | 15.1% | yes / no |
| 10028 | power_fisher | 10.8992 | 0.5071% | 13.9% | yes / no |
| 10028 | power_bispectrum_matched | 7.73215 | 1.03% | 13.5% | yes / no |
| 10028 | power_bispectrum_fisher | 10.3575 | 0.6388% | 14.2% | yes / no |
| 10049 | power_matched | 0.671769 | 6.864% | 2.88% | yes / yes |
| 10049 | power_fisher | 0.591083 | 7.756% | 2.71% | yes / yes |
| 10049 | power_bispectrum_matched | 0.635435 | 7.678% | 2.87% | yes / yes |
| 10049 | power_bispectrum_fisher | 0.570831 | 8.426% | 2.7% | yes / yes |
| 10076 | power_matched | 4.46231 | 0.2634% | 0.456% | yes / no |
| 10076 | power_fisher | 4.31758 | 0.2718% | 0.464% | yes / no |
| 10076 | power_bispectrum_matched | 4.37733 | 0.2739% | 0.467% | yes / no |
| 10076 | power_bispectrum_fisher | 4.21943 | 0.2815% | 0.468% | yes / no |

The amplitude interval is [.9,1.1]. These are **mean-score** envelopes, not the original event-probability null bounds. For 10028 and 10076, view-dependent amplitude can cover the 25%-deletion mean gap for every original direction. Even without that ambiguity, roughly 0.3% white-noise variance error matches the 10076 deletion gap. The bispectrum mean has zero noise shift only under the specified independent nonredundant Gaussian coordinates. The calculation does not quantify arbitrary colored or correlated noise. The experimental [metadata inventory](STACK-NUISANCE-INVENTORY-RESULTS.md) does not establish this precision.

## Existing combined score components

These are the unchanged coefficients of each combined score, with one block zeroed for diagnosis. Neither block is refit or promoted to a new test. Their separations do not add because their noise and view covariances are nonzero.

| EMPIAR | Component | D at 25% deletion |
|---|---|---:|
| 10028 | power_bispectrum_matched_power_component | 0.00133703 |
| 10028 | power_bispectrum_matched_bispectrum_component | 0.000380396 |
| 10028 | power_bispectrum_fisher_power_component | 0.00137027 |
| 10028 | power_bispectrum_fisher_bispectrum_component | 0.000387632 |
| 10049 | power_bispectrum_matched_power_component | 0.0015915 |
| 10049 | power_bispectrum_matched_bispectrum_component | 0.000172497 |
| 10049 | power_bispectrum_fisher_power_component | 0.00160074 |
| 10049 | power_bispectrum_fisher_bispectrum_component | 0.000148738 |
| 10076 | power_bispectrum_matched_power_component | 0.000118677 |
| 10076 | power_bispectrum_matched_bispectrum_component | 4.53958e-06 |
| 10076 | power_bispectrum_fisher_power_component | 0.000118651 |
| 10076 | power_bispectrum_fisher_bispectrum_component | 3.17119e-06 |

## Verification and remaining work

All 60 original amplitude-one score prefixes replay within 1.01e-11. Independent complex Gaussian monomial enumeration checks 216 conditional variances across three noise levels and three fixed image indices; the largest difference is 1.78e-15. Direct Euclidean Gaussian residuals plus SciPy log-sum-exp reproduce 576 saved likelihood ratios within 8.53e-14, including quadrature ESS and maximum weights. All saved summary separations are replayed. These are floating-point checks, not rigorous integration error bounds.

The full ledger still lacks the original **event** envelope decomposition, finite-calibration slack at matched draws, measured experimental noise spectra and realistic-nuisance likelihood validation. The surviving scientific question is whether a tractable nuisance-aware likelihood can support a clearly specified uncertainty target. No new score fitting, calibration repetition, acceptance-level claim or full-review request follows from these numerical comparisons.

Sources: `scripts/analyze_matched_information.py`, `scripts/analyze_matched_haar_likelihood.py`, `scripts/verify_matched_information.py`; complete summaries and raw arrays are under `results/uncertainty/development/matched-*-v1`. The reporting CSV includes all four deletion fractions, both banks and every declared prefix. Earlier arrays and failed method branches remain unchanged.
