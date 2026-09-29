# Wider continuous pose budgets: development sensitivity check

29 September 2026. The frozen continuous-v1 study remains unchanged.

The development continuous audit and frozen study use rotation limits 0.1 and
0.5 degrees. These are sensitivity assumptions, not empirically calibrated
maximum errors of deposited poses. Recent orientation papers report errors
on larger angular scales, with their own reference-label and symmetry caveats.
We must not present the small-budget experiment as encompassing typical
unknown-pose inference.

Run the existing continuous nonlinear audit on the same saved development
weights for all three stacks, both width-0.07 features, and budgets 1, 2 and
5 degrees. Keep the shift radius and all other settings unchanged. Record raw
half-widths and the selected no-data fallback, including every vacuous result.
Check the exact continuous density-ball bias at the coherent allowed rotation.
These 18 post-hoc settings are a development sensitivity analysis, not new
confirmation observations or a calibrated model of experimental pose errors.

Command:

```sh
OPENBLAS_NUM_THREADS=4 .venv/bin/python scripts/audit_uq_continuous_pose.py \
  --angles 1,2,5 --output continuous-pose-large
```

The scientific question is whether the present post-audit remains informative
over a broader declared nuisance class. A loose upper bound does not prove the
feature is unidentifiable: compare its feasible lower bound and acknowledge
possible conservatism of the polynomial/spectral relaxation.

## Feasible nonlinear adversaries

The larger-radius sweep motivates a separately labeled adversary search. It
uses an order-16 Gauss quadrature surrogate and float32 MPS projected Adam to
propose poses, with three starts (zero, coherent rotation, random boundary),
both signed objectives, and 150 steps at learning rate 0.05. Test both broad
targets on all three development geometries at 0.5, 1 and 2 degrees. No frozen
continuous-v1 weights, cases, or decision rules change.

For every final proposed pose, recompute its pilot bias using exact constant-
cell Fourier projections and eliminate the entire continuous L2 density ball
with analytic sinc integrals in float64. Save poses, surrogate traces, exact
numerical biases, source hashes and independent direct-field checks. The
quadrature density approximation is only an optimization surrogate; a surrogate
objective is never reported as a feasible continuous lower bound.

A preliminary order-12, 20-step probe on 10028's central target at 0.5 degrees
passed independent Fourier-field checks (relative errors below 6e-6). Its
largest exact numerical feasible bias was 1.3325 versus the uniform upper
4.9854 and coherent-pose bias 1.2378. This probe motivated higher quadrature
order and more iterations; neither the probe nor the full search is a
preregistered confirmation study. Final numerical checks are not interval
arithmetic and do not certify the global nonlinear maximum.

```sh
OPENBLAS_NUM_THREADS=4 .venv/bin/python scripts/stress_uq_continuous_pose.py \
  --order 16
```
