# Completed joint density/pose cubic design

30 September 2026. The declared 10049 region-1 fit completes in 2,320.89
seconds at 2.65 GB peak resident memory. The original thirteen-column span,
128 particles, 20 A Gaussian target, density radius 2 and one-degree/0.5 A
joint pose budget remain unchanged. This is one development case, not an
independent confirmation or an experiment establishing these physical budgets.

Two design evaluations select evaluation 2. The final-cut-rescored restricted
guide gap is .00088459; this does not certify full-space convergence. The first
master reports `optimal_inaccurate`, the second `optimal`; both remain recorded.
Design takes 571.95 seconds and the fresh final certificate 1,733.59 seconds.
Four probes and forty power steps use the predeclared fresh certificate seed.

The audited half-width is 2.130112, or .175192 of the same no-data half-width
12.158699. This is 2.55% smaller than the original .179785. The SD under the
supplied simulation model is .218692 and the bias upper bound is 1.553146.
At these same selected weights, the preceding joint bound gives 1.595052 and
the triangle bound 1.913367. These are bound comparisons, not separately
optimized estimators. The resolvent calculation does not fall back. Its
selected shift reaches the declared lower log-grid boundary; this is retained
and does not justify extrapolating an improvement beyond the search.

All six conditional reference checks are retained. Minimum sign power is
.014333 in the original frame and .978287 in the registered frame. No tested
conditional coverage failure occurs. These reference generators and supplied
noise cannot establish experimental coverage, pose calibration, or realistic
fine-detail recovery. Numerical guards are not rigorous floating-point
enclosures, and no full-space optimality claim follows from the restricted gap.

The result and exact arrays are in `joint-cubic-design-probe`; the integrity
record verifies 74 source snapshots and both improving checkpoint arrays.
The separately declared real-image application is reported in
[JOINT-ENRICHED-EXPERIMENTAL-RESULTS.md](JOINT-ENRICHED-EXPERIMENTAL-RESULTS.md).
