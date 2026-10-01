# Finite-view oracle gate: complete with numerical failures

All 288 planned cases were attempted. Two CLARABEL cases fail and three report `optimal_inaccurate`; all remain in the records. Completed cases have full-domain dual gaps no larger than 5.4e-7 per profile. The 64-view true-map control has no positive selected expected log growth. Earlier formatting failures and the interrupted 139-case numerical attempt are retained. No continuous-orientation validity or measured rejection power is claimed.

| Stack | Candidate | Variance upper | Expected log, 128 particles: feasible lower / dual upper | Failed profiles |
|---|---|---|---|---|
| 10028 | true_map | 1 | 0 / 2.87709e-08 | 0 |
| 10028 | true_map | 2 | 0 / 3.34511e-13 | 0 |
| 10028 | true_map | 4 | 0 / 2.85284e-13 | 0 |
| 10028 | region_half_removed | 1 | 0.252475 / 0.252475 | 0 |
| 10028 | region_half_removed | 2 | 0.0850745 / 0.0850746 | 0 |
| 10028 | region_half_removed | 4 | 0.0264558 / 0.0264558 | 0 |
| 10028 | region_removed | 1 | 1.06779 / 1.06779 | 0 |
| 10028 | region_removed | 2 | 0.362146 / 0.362146 | 0 |
| 10028 | region_removed | 4 | 0.113318 / 0.113318 | 0 |
| 10028 | zero_signal | 1 | 4454.02 / 4454.02 | 0 |
| 10028 | zero_signal | 2 | 1431.35 / 1431.35 | 0 |
| 10028 | zero_signal | 4 | 419.033 / 419.033 | 0 |
| 10049 | true_map | 1 | 0 / 8.08204e-17 | 0 |
| 10049 | true_map | 2 | 0 / 2.71683e-13 | 0 |
| 10049 | true_map | 4 | 0 / 8.02817e-14 | 0 |
| 10049 | region_half_removed | 1 | 0.0434357 / 0.0434358 | 0 |
| 10049 | region_half_removed | 2 | 0.0114463 / 0.0114464 | 0 |
| 10049 | region_half_removed | 4 | 0.00294242 / 0.00294247 | 0 |
| 10049 | region_removed | 1 | 0.22992 / 0.229921 | 0 |
| 10049 | region_removed | 2 | 0.0520214 / unavailable | 1 |
| 10049 | region_removed | 4 | 0.013447 / unavailable | 1 |
| 10049 | zero_signal | 1 | 63.7135 / 63.7135 | 0 |
| 10049 | zero_signal | 2 | 15.9799 / 15.9799 | 0 |
| 10049 | zero_signal | 4 | 3.99824 / 3.99824 | 0 |
| 10076 | true_map | 1 | 0 / 7.6949e-14 | 0 |
| 10076 | true_map | 2 | 0 / 2.33683e-14 | 0 |
| 10076 | true_map | 4 | 0 / 6.30038e-14 | 0 |
| 10076 | region_half_removed | 1 | 0.0167517 / 0.0167604 | 0 |
| 10076 | region_half_removed | 2 | 0.00434154 / 0.00434161 | 0 |
| 10076 | region_half_removed | 4 | 0.00110543 / 0.00110551 | 0 |
| 10076 | region_removed | 1 | 0.0771069 / 0.0771071 | 0 |
| 10076 | region_removed | 2 | 0.0200324 / 0.0200325 | 0 |
| 10076 | region_removed | 4 | 0.00510886 / 0.00510893 | 0 |
| 10076 | zero_signal | 1 | 137.126 / 137.126 | 0 |
| 10076 | zero_signal | 2 | 36.7059 / 36.7059 | 0 |
| 10076 | zero_signal | 4 | 9.38213 / 9.38213 | 0 |

Local-removal growth is small even with a 64-view null and an oracle alternative. At unit covariance bound the full-removal upper values are 1.068, .230 and .0771, all below log(20) at 128 particles. Increasing the particle count scales these *expectations*, not a demonstrated rejection probability. Enlarging the null view catalog may remove the separation entirely. The next bounded diagnostic should check whether the alternative mean power belongs to a larger candidate cone before attempting a continuous angular certificate.
