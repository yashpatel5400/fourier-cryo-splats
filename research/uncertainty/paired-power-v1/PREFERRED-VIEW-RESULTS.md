# Fresh preferred-view controls: complete results

All 18 law/stack runs finish (nine laws on each of two stacks). Each law has 65,536 new independent views, with density kappa on an axis-specific Haar cap union. There are 1,179,648 independently drawn rotations in total; maps, amplitudes and scores share draws. The physical amplitude is independent of the measurement noise in this generator. All cases in the frozen protocol are retained, including both correct-null maps and all .9/1/1.1 amplitudes.

## Removed-candidate rejection on true-map data at n=10,000

The table reports the range of point projections across all three spatial axes, using the paired view-variance procedure. These ranges are not confidence intervals or power guarantees over all viewing laws. Pointwise 95% intervals and every method are in the complete CSV. Actual rejected groups of size 1,024 are reported separately.

| Stack | Score | kappa | a=.9 | a=1 | a=1.1 |
|---|---|---:|---:|---:|---:|
| 10028 | power_range09_11 | 1.1 | 0.0001069--0.0003568 | 0.001116--0.002363 | 0.006741--0.01542 |
| 10028 | power_range09_11 | 2 | 2.938e-13--3.075e-12 | 4.266e-11--4.559e-10 | 2.382e-09--3.112e-08 |
| 10028 | power_range09_11 | 5 | 6.339e-36--1.484e-29 | 2.807e-32--1.783e-26 | 1.69e-28--6.973e-24 |
| 10028 | power_bispectrum_range09_11 | 1.1 | 8.551e-06--2.678e-05 | 0.003536--0.008267 | 0.04966--0.09577 |
| 10028 | power_bispectrum_range09_11 | 2 | 4.828e-15--3.977e-14 | 5.955e-10--7.595e-09 | 1.077e-07--3.481e-06 |
| 10028 | power_bispectrum_range09_11 | 5 | 8.063e-39--1.94e-34 | 2.718e-29--2.615e-26 | 2.039e-23--1.479e-20 |
| 10049 | power_range09_11 | 1.1 | 0.0005069--0.007669 | 0.1507--0.3806 | 0.9016--0.9758 |
| 10049 | power_range09_11 | 2 | 1.096e-09--7.498e-08 | 2.514e-05--0.001088 | 0.04046--0.2551 |
| 10049 | power_range09_11 | 5 | 6.656e-28--1.93e-17 | 5.027e-19--2.197e-10 | 9.098e-12--0.0001022 |
| 10049 | power_bispectrum_range09_11 | 1.1 | 0.003565--0.00579 | 0.3123--0.3859 | 0.9836--0.9881 |
| 10049 | power_bispectrum_range09_11 | 2 | 7.729e-10--4.34e-07 | 6.047e-05--0.007735 | 0.07391--0.5218 |
| 10049 | power_bispectrum_range09_11 | 5 | 2.772e-30--1.686e-19 | 1.078e-20--2.008e-11 | 1.602e-12--5.06e-05 |

## Correct-null and numerical controls

Across all 648 correct-null method/law/amplitude/score cells, the largest observed rejection count is 2/64 groups of 1,024 particles. The associated exact 95% interval is wide; absence of an observed excess is not a proof of unconditional 5% control. Calibration was held fixed, and cells share draws, so pooling these outcomes as independent trials would be wrong. At n=10,000, the largest correct-null point projection across all procedures is 0.000258882.

The full CSV contains 1944 binomial projections and a separate CSV contains 648 pairs of actual group-rejection counts. Every 64-vector of observed counts remains in JSON. The independent verifier replays 216 event counts, 2592 critical values and 1296 group-count vectors; the first 216 direct scores agree within 7.96e-12. Thirty-six physical-cell first-view checks agree within 1.02e-11 per Fourier coordinate. Frequencies are unique, nonzero and have no antipodal pairs, and the first transfer profile matches its saved source exactly. Only these declared score/operator checks were independently replayed; this is not interval arithmetic.

The two overlapping Mac runs take 417.51 and 413.02 seconds. Their sum is not elapsed wall time.

## Interpretation

The prescribed preferential viewing does not trigger an observed type-I failure in these controls. Power changes materially with amplitude, so amplitude-one Haar results do not summarize this experiment. The tests remain weak under strong view concentration. Axis caps are not the extremal event-probability tilt; an adversarial-law study remains outstanding. The known proper noise, exact transfer profile, supplied viewing-density cap and oracle full-removal score are unchanged limitations. No experimental viewing law, uncertainty coverage or three-stack success is established.
