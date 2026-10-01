# Saved event-probability decomposition

1 October 2026 UTC. This is an algebraic diagnosis of the existing Fisher-study test, not a new score, calibration or power run. All 60 score/cap combinations are retained in the CSV. The table shows κ=1.1, which the experimental metadata does not support as a measured cap.

| EMPIAR | Score | Observed 25% probability gap | Grouped envelope − held-out null | Outer mean increment | View-cap increment | Final bound − observed alternative |
|---|---|---:|---:|---:|---:|---:|
| 10028 | power_matched | 0.0147095 | 0.0208664 | 0.00459459 | 0.051173 | 0.0619245 |
| 10028 | power_fisher | 0.0150604 | 0.0203323 | 0.00445089 | 0.0490411 | 0.0587639 |
| 10028 | power_bispectrum_matched | 0.0126266 | 0.0226164 | 0.00444236 | 0.0492134 | 0.0636455 |
| 10028 | power_bispectrum_fisher | 0.0159912 | 0.0221343 | 0.0044736 | 0.0494626 | 0.0600793 |
| 10049 | power_matched | 0.0158081 | 0.0131259 | 0.00279408 | 0.0231532 | 0.0232651 |
| 10049 | power_fisher | 0.0159378 | 0.0126305 | 0.00276643 | 0.022671 | 0.0221301 |
| 10049 | power_bispectrum_matched | 0.0161896 | 0.0145521 | 0.0027939 | 0.0231307 | 0.0242871 |
| 10049 | power_bispectrum_fisher | 0.0173187 | 0.0143967 | 0.00276329 | 0.0225773 | 0.0224185 |
| 10076 | power_matched | 0.00402832 | 0.0105524 | 0.00222891 | 0.0114553 | 0.0202083 |
| 10076 | power_fisher | 0.00405121 | 0.0101051 | 0.00222975 | 0.0115103 | 0.0197939 |
| 10076 | power_bispectrum_matched | 0.00431061 | 0.0105243 | 0.00223267 | 0.0115408 | 0.0199871 |
| 10076 | power_bispectrum_fisher | 0.00432587 | 0.0111637 | 0.00223196 | 0.0115184 | 0.0205882 |

The three increments plus the held-out null rate equal the archived final bound to numerical precision. The last increment is the net result of the original minimum and clipping; it is not independently certified physical viewing uncertainty. The positive final gaps explain why the original test is insensitive at this alternative under the specified simulator. Sampling uncertainty in the held-out probabilities remains as reported in the original study.

**The second column of increments is not an estimate of finite-L bias.** It combines view-dependent amplitude, the cellwise event supremum, maximization after finite noise replication, and differences between calibration/test Monte Carlo draws. The old L=32 versus L=128 runs also used different views and noise, so their difference cannot isolate any one component. The complete CSV separately retains the global three-amplitude grid maximum, individual-noise envelope, mean radius, paired-product sampling/range terms and CVaR comparators. The three-amplitude maximum does not optimize a continuous amplitude law.

This completes the identifiable saved-data portion of the requested decomposition. A separate matched nested calculation would be needed to estimate inner-replication bias. That calculation cannot repair the demonstrated noise sensitivity, poor retained separation or unsupported experimental viewing model; the moment-test branch remains frozen.

Inputs and source hashes are in `results/uncertainty/development/matched-event-decomposition-v1/summary.json`.
