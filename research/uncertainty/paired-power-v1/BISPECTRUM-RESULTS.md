# Translation-invariant moment diagnostics: complete outcomes

This is post-review method development, not experimental calibration. All 24 declared cells finish in 38.124 seconds: twelve exact true-map controls and twelve removal fits. Every removal fit reaches the 2,500-iteration limit; none meets the 1e-5 distance-gap criterion. Feasible mixtures and separating directions bound the finite-hull distance independently of convergence.
The moments use the first saved CTF/noise profile and an oracle alternative averaged over 4,160 orientations. Power plus 512 complex bispectra has 1,244 real coordinates. Amplitude assumptions differ from the earlier unrestricted-amplitude bilinear e-value model.

| Stack | Features | Amplitude | Distance lower | Distance upper | Separation after other 10,000 views | After continuous searches | Finite-catalog sufficient n |
|---|---|---|---:|---:|---:|---:|---:|
| 10028 | power | 1 | 0.172361 | 0.181022 | 0.152398 | 0.14652 | 8813.52 |
| 10028 | power | [.9,1.1] | 0.0230997 | 0.0401629 | 9.81945e-05 | -0.00369244 | 2.26859e+10 |
| 10028 | power_bispectrum | 1 | 0.214318 | 0.222109 | 0.127162 | 0.102879 | 16653.3 |
| 10028 | power_bispectrum | [.9,1.1] | 0.0473788 | 0.0675974 | -0.0376435 | -0.0410645 | none |
| 10049 | power | 1 | 0.0541664 | 0.0554906 | 0.0456263 | 0.0436558 | 23163.8 |
| 10049 | power | [.9,1.1] | 0.024014 | 0.0272112 | 0.018051 | 0.017434 | 139784 |
| 10049 | power_bispectrum | 1 | 0.0693988 | 0.0697392 | 0.0525673 | 0.0507924 | 18221.1 |
| 10049 | power_bispectrum | [.9,1.1] | 0.0418072 | 0.0429289 | 0.0305114 | 0.0296193 | 56011.2 |
| 10076 | power | 1 | 0 | 3.74227e-05 | 0 | no direction | none |
| 10076 | power | [.9,1.1] | 0 | 3.10761e-05 | 0 | no direction | none |
| 10076 | power_bispectrum | 1 | 0 | 0.00270974 | 0 | no direction | none |
| 10076 | power_bispectrum | [.9,1.1] | 0 | 0.00115923 | 0 | no direction | none |

**The last column is not a continuous-pose sample-size guarantee.** It is the declared two-Cantelli sufficient bound for the enlarged finite rotation catalog, continuously profiled over the amplitude interval, before the adversarial searches. It is not a necessary sample size, lower complexity bound or empirical power measurement. The negative 10028 margins supersede any tempting favorable interpretation of finite-catalog diagnostics.

## Continuous orientations and independent calculations

All 160 declared searches complete (159 optimizers report success) in 25.654 seconds. Four all-zero 10076 directions are explicitly skipped. Each selected pose is recomputed by the physical cell Fourier operator. The eight worst selected views also agree with independent direct cell sums within 3.79e-12. The interpolation guide never serves as a continuous certificate.

A separate polynomial-moment enumerator, independent of the real-derivative implementation, checks all 48 declared conditional variance comparisons across eight nonzero contrasts, three means and two amplitudes. The maximum absolute discrepancy is 3.55e-15. All twelve saved removal mixtures and distance brackets also replay. These checks establish numerical agreement under the assumed noise law, not its physical calibration.

Both ranged-amplitude contrasts on 10028 have null orientations whose expected score exceeds the alternative expectation. The fixed-amplitude 10028 contrasts and all four 10049 contrasts retain positive margins at the searched poses, which supplies no upper bound over all orientations. On 10049, the ranged bispectrum margin is about .0296 after search, compared with .0174 for power only. This is an oracle directional comparison with different noise variance and nonzero fitting gaps, not demonstrated testing power or a new reconstruction result. No positive separator was found for 10076; its nonzero upper residuals and fitting limits preclude a claim that the full third-moment family cannot distinguish the structures.

## Decision and open requirements

Do not promote any saved direction to a calibrated continuous-view test. 10049 warrants a genuinely global continuous-orientation bound or an explicitly narrower justified model before noise/power trials. The first profile, assumed independent circular Gaussian noise, map normalization and amplitude interval still require physical validation. The present Cantelli calculation needs no exponential moment; a cubic Gaussian statistic cannot simply replace a bilinear term in the earlier e-value. The exact Hermite variance, convex projection and moment invariants are classical mathematics, not established novelty. All three full Fable reviews remain rejections; no fourth full acceptance review occurred.
