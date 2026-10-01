# Complete local-pose refitting study

Exactly 200 independent noisy datasets were attempted for each of three fixed truth/acquisition geometries. All methods, templates, targets and image controls within a dataset are paired. Every failed trial remains in the denominator; interval width summaries describe available results. Exact binomial Monte Carlo intervals and all paired discordance counts are retained in the JSON/CSV summaries.

## EMPIAR-10028

Completed 200/200; failures 0. Runtime 6718.2 seconds; peak resident 0.681 GB.

| Template | Target | Method | Same raw / selected coverage | Independent raw / selected coverage | Median raw / selected relative width | Fallback rate | Independent sign exclusion |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| O | C | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.0867/0.0867 | 0.000 | 1.000 |
| O | C | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.0887/0.0887 | 0.000 | 1.000 |
| O | C | Bounded pose | 1.000/1.000 | 1.000/1.000 | 775/1 | 1.000 | 0.000 |
| O | C | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 761/1 | 1.000 | 0.000 |
| O | C | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 922/1 | 1.000 | 0.000 |
| O | C | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.136/0.136 | 0.000 | 0.890 |
| O | C | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.136/0.136 | 0.000 | 0.890 |
| O | C | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0729/0.0729 | 0.000 | 1.000 |
| O | C | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0729/0.0729 | 0.000 | 1.000 |
| O | Z | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.105/0.105 | 0.000 | 0.000 |
| O | Z | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.108/0.108 | 0.000 | 0.000 |
| O | Z | Bounded pose | 1.000/1.000 | 1.000/1.000 | 827/1 | 1.000 | 0.000 |
| O | Z | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 811/1 | 1.000 | 0.000 |
| O | Z | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 982/1 | 1.000 | 0.000 |
| O | Z | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.168/0.168 | 0.000 | 0.000 |
| O | Z | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.168/0.168 | 0.000 | 0.000 |
| O | Z | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0892/0.0892 | 0.000 | 0.000 |
| O | Z | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0893/0.0893 | 0.000 | 0.000 |
| P | C | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.0872/0.0872 | 0.000 | 1.000 |
| P | C | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.0893/0.0893 | 0.000 | 1.000 |
| P | C | Bounded pose | 1.000/1.000 | 1.000/1.000 | 848/1 | 1.000 | 0.000 |
| P | C | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 834/1 | 1.000 | 0.000 |
| P | C | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 1.01e+03/1 | 1.000 | 0.000 |
| P | C | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.137/0.137 | 0.000 | 0.820 |
| P | C | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.137/0.137 | 0.000 | 0.830 |
| P | C | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0733/0.0733 | 0.000 | 1.000 |
| P | C | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0733/0.0733 | 0.000 | 1.000 |
| P | Z | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.107/0.107 | 0.000 | 0.000 |
| P | Z | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.109/0.109 | 0.000 | 0.000 |
| P | Z | Bounded pose | 1.000/1.000 | 1.000/1.000 | 907/1 | 1.000 | 0.000 |
| P | Z | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 891/1 | 1.000 | 0.000 |
| P | Z | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 1.08e+03/1 | 1.000 | 0.000 |
| P | Z | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.169/0.169 | 0.000 | 0.000 |
| P | Z | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.169/0.169 | 0.000 | 0.000 |
| P | Z | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.09/0.09 | 0.000 | 0.000 |
| P | Z | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0901/0.0901 | 0.000 | 0.000 |
## EMPIAR-10049

Completed 200/200; failures 0. Runtime 6648.4 seconds; peak resident 0.660 GB.

| Template | Target | Method | Same raw / selected coverage | Independent raw / selected coverage | Median raw / selected relative width | Fallback rate | Independent sign exclusion |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| O | C | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.113/0.113 | 0.000 | 1.000 |
| O | C | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.118/0.118 | 0.000 | 1.000 |
| O | C | Bounded pose | 1.000/1.000 | 1.000/1.000 | 1e+03/1 | 1.000 | 0.000 |
| O | C | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 988/1 | 1.000 | 0.000 |
| O | C | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 1.2e+03/1 | 1.000 | 0.000 |
| O | C | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.157/0.157 | 0.000 | 0.800 |
| O | C | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.158/0.158 | 0.000 | 0.810 |
| O | C | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0879/0.0879 | 0.000 | 1.000 |
| O | C | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0882/0.0882 | 0.000 | 1.000 |
| O | Z | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.127/0.127 | 0.000 | 0.000 |
| O | Z | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.132/0.132 | 0.000 | 0.000 |
| O | Z | Bounded pose | 1.000/1.000 | 1.000/1.000 | 987/1 | 1.000 | 0.000 |
| O | Z | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 973/1 | 1.000 | 0.000 |
| O | Z | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 1.18e+03/1 | 1.000 | 0.000 |
| O | Z | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.189/0.189 | 0.000 | 0.000 |
| O | Z | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.189/0.189 | 0.000 | 0.000 |
| O | Z | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.103/0.103 | 0.000 | 0.005 |
| O | Z | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.103/0.103 | 0.000 | 0.005 |
| P | C | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.114/0.114 | 0.000 | 1.000 |
| P | C | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.119/0.119 | 0.000 | 1.000 |
| P | C | Bounded pose | 1.000/1.000 | 1.000/1.000 | 1.04e+03/1 | 1.000 | 0.000 |
| P | C | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 1.03e+03/1 | 1.000 | 0.000 |
| P | C | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 1.24e+03/1 | 1.000 | 0.000 |
| P | C | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.159/0.159 | 0.000 | 0.660 |
| P | C | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.159/0.159 | 0.000 | 0.670 |
| P | C | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0887/0.0887 | 0.000 | 1.000 |
| P | C | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0891/0.0891 | 0.000 | 1.000 |
| P | Z | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.129/0.129 | 0.000 | 0.000 |
| P | Z | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.133/0.133 | 0.000 | 0.000 |
| P | Z | Bounded pose | 1.000/1.000 | 1.000/1.000 | 1.02e+03/1 | 1.000 | 0.000 |
| P | Z | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 1.01e+03/1 | 1.000 | 0.000 |
| P | Z | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 1.22e+03/1 | 1.000 | 0.000 |
| P | Z | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.19/0.19 | 0.000 | 0.000 |
| P | Z | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.19/0.19 | 0.000 | 0.000 |
| P | Z | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.103/0.103 | 0.000 | 0.000 |
| P | Z | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.104/0.104 | 0.000 | 0.000 |
## EMPIAR-10076

Completed 200/200; failures 0. Runtime 6217.8 seconds; peak resident 0.669 GB.

| Template | Target | Method | Same raw / selected coverage | Independent raw / selected coverage | Median raw / selected relative width | Fallback rate | Independent sign exclusion |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| O | C | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.106/0.106 | 0.000 | 1.000 |
| O | C | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.111/0.111 | 0.000 | 1.000 |
| O | C | Bounded pose | 1.000/1.000 | 1.000/1.000 | 651/1 | 1.000 | 0.000 |
| O | C | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 642/1 | 1.000 | 0.000 |
| O | C | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 777/1 | 1.000 | 0.000 |
| O | C | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.152/0.152 | 0.000 | 0.985 |
| O | C | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.152/0.152 | 0.000 | 0.990 |
| O | C | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0834/0.0834 | 0.000 | 1.000 |
| O | C | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0835/0.0835 | 0.000 | 1.000 |
| O | Z | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.126/0.126 | 0.000 | 0.000 |
| O | Z | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.131/0.131 | 0.000 | 0.000 |
| O | Z | Bounded pose | 1.000/1.000 | 1.000/1.000 | 679/1 | 1.000 | 0.000 |
| O | Z | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 668/1 | 1.000 | 0.000 |
| O | Z | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 809/1 | 1.000 | 0.000 |
| O | Z | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.185/0.185 | 0.000 | 0.000 |
| O | Z | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.185/0.185 | 0.000 | 0.000 |
| O | Z | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.101/0.101 | 0.000 | 0.000 |
| O | Z | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.101/0.101 | 0.000 | 0.000 |
| P | C | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.106/0.106 | 0.000 | 1.000 |
| P | C | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.111/0.111 | 0.000 | 1.000 |
| P | C | Bounded pose | 1.000/1.000 | 1.000/1.000 | 685/1 | 1.000 | 0.000 |
| P | C | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 676/1 | 1.000 | 0.000 |
| P | C | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 818/1 | 1.000 | 0.000 |
| P | C | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.152/0.152 | 0.000 | 1.000 |
| P | C | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.152/0.152 | 0.000 | 1.000 |
| P | C | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.0836/0.0836 | 0.000 | 1.000 |
| P | C | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.0837/0.0837 | 0.000 | 1.000 |
| P | Z | Fixed folded | 1.000/1.000 | 1.000/1.000 | 0.126/0.126 | 0.000 | 0.000 |
| P | Z | Fixed sum | 1.000/1.000 | 1.000/1.000 | 0.131/0.131 | 0.000 | 0.000 |
| P | Z | Bounded pose | 1.000/1.000 | 1.000/1.000 | 718/1 | 1.000 | 0.000 |
| P | Z | Mixed 0 | 1.000/1.000 | 1.000/1.000 | 707/1 | 1.000 | 0.000 |
| P | Z | Mixed .1 | 1.000/1.000 | 1.000/1.000 | 856/1 | 1.000 | 0.000 |
| P | Z | GP 2 fixed | 1.000/1.000 | 1.000/1.000 | 0.185/0.185 | 0.000 | 0.000 |
| P | Z | GP 2 pose | 1.000/1.000 | 1.000/1.000 | 0.185/0.185 | 0.000 | 0.000 |
| P | Z | GP 1 fixed | 1.000/1.000 | 1.000/1.000 | 0.101/0.101 | 0.000 | 0.000 |
| P | Z | GP 1 pose | 1.000/1.000 | 1.000/1.000 | 0.101/0.101 | 0.000 | 0.000 |

## Pose calibration and numerical convergence

| Stack | Template | Score bound | Median rotation RMS (degrees) | Median shift RMS (A) | Pose-ball inclusion (exact 95% MC interval) |
| --- | --- | ---: | ---: | ---: | --- |
| 10028 | O | 29.840 | 11.990 | 6.617 | 0.965 (0.929, 0.986) |
| 10028 | P | 31.150 | 15.568 | 8.051 | 0.995 (0.972, 1.000) |
| 10049 | O | 31.063 | 15.766 | 6.822 | 1.000 (0.982, 1.000) |
| 10049 | P | 31.580 | 16.013 | 7.419 | 1.000 (0.982, 1.000) |
| 10076 | O | 31.097 | 16.103 | 8.639 | 0.995 (0.972, 1.000) |
| 10076 | P | 31.849 | 18.481 | 9.282 | 0.995 (0.972, 1.000) |

10028: 4000 converged solves of 4000; every numerical result remains in the coverage evaluation.

10049: 4000 converged solves of 4000; every numerical result remains in the coverage evaluation.

10076: 4000 converged solves of 4000; every numerical result remains in the coverage evaluation.

## Nonlinear-remainder decomposition

Post-outcome reporting diagnostic of every completed fit, with no refitting or selection. Fractions below concern the raw mixed-common-zero width.

| Stack | Template | Target | Fits | Median remainder / raw width | Minimum | Maximum |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 10028 | oracle_reference | center | 200 | 0.994057 | 0.993908 | 0.994175 |
| 10028 | oracle_reference | contrast | 200 | 0.993334 | 0.993171 | 0.993528 |
| 10028 | independent_pilot | center | 200 | 0.994170 | 0.993950 | 0.994355 |
| 10028 | independent_pilot | contrast | 200 | 0.993390 | 0.993056 | 0.993590 |
| 10049 | oracle_reference | center | 200 | 0.994881 | 0.993835 | 0.995291 |
| 10049 | oracle_reference | contrast | 200 | 0.994564 | 0.993544 | 0.994991 |
| 10049 | independent_pilot | center | 200 | 0.994892 | 0.993859 | 0.995307 |
| 10049 | independent_pilot | contrast | 200 | 0.994582 | 0.993739 | 0.995066 |
| 10076 | oracle_reference | center | 200 | 0.994985 | 0.994800 | 0.995122 |
| 10076 | oracle_reference | contrast | 200 | 0.994236 | 0.993875 | 0.994458 |
| 10076 | independent_pilot | center | 200 | 0.995021 | 0.994804 | 0.995180 |
| 10076 | independent_pilot | contrast | 200 | 0.994306 | 0.994103 | 0.994531 |

A no-data fallback can provide coverage without using inference images or resolving a sign. The deterministic pose construction is simulation-assisted and its marginal tolerance guarantee averages over calibration datasets. Mixed intervals retain unverified conditional centering for estimated designs. These results do not calibrate experimental density coverage or repeat global ab initio reconstruction.
