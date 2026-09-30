# Known-map pose information: completed diagnostic

30 September 2026. Protocol, source and derivative checks were published in
commit `12a55e8` before execution. All six cases complete: two specified known
signals and 128 particles on each of the three stacks. All 768 five-parameter
information matrices meet the declared numerical full-rank threshold. Every
matrix, inverse and coordinate scale is retained, not only these medians.

| EMPIAR | Reference rotation RMS (degrees) | Reference shift RMS (Å) | Pilot rotation RMS (degrees) | Pilot shift RMS (Å) |
|---|---:|---:|---:|---:|
| 10028 | 2.443 | 2.090 | 3.299 | 3.303 |
| 10049 | 9.121 | 3.216 | 4.737 | 1.856 |
| 10076 | 6.214 | 5.089 | 7.576 | 5.105 |

These are square roots of traces of covariance blocks of `(JᵀJ)⁻¹`, where J
uses three rotation coordinates in degrees and two translations in Å. Each
particle's other coordinates are unknown in this joint local calculation.
The reference and Gaussian-pilot signals are separately normalized to unit
continuous L2 norm, with the previously supplied Gaussian simulation noise.
Neither is the unknown experimental density. The raw records also retain
maximum directional standard deviations and the known-translation calculation.

**Scope of the conclusion.** The radius-twelve frequency band ends at 40.2,
19.68 and 34.93 Å, respectively. It contains limited information about a single
particle's pose in this simulation, even when its map is known. The original
consensus alignments used other information and higher frequencies. These
numbers cannot estimate, upper-bound, or invalidate their actual errors.
The classical Cramér–Rao variance comparison requires regularity and local
unbiasedness; the inverse information also describes the linear Gaussian
tangent experiment. Neither interpretation supplies a nonlinear confidence
ball, guarantees an efficient estimator, or resolves global orientation modes.
This diagnostic therefore does not calibrate our one-/two-degree sensitivity
budgets. It is not an uncertainty-method improvement or evidence of coverage.

Input, code, protocol and per-case output identities are in
`results/uncertainty/development/oracle-pose-information-v1/`. Analytic cell
derivatives include cell extent and agree with independent nonlinear finite
differences in the two dedicated tests. The earlier test-collection failure
and repaired test run remain recorded in the verification manifest. The new
tests were run separately after the previous full 204-test suite; there is no
claim that a 206-test full suite has been run.
