# Development coverage summary

Source: `../coverage-benchmark.json`; 48 independent design settings, with paired
noise draws across methods within each signal. These are exact-dictionary,
fixed-pilot, known-bound simulations. They are not real-data calibration.

Coverage is conditional on each fixed signal; minimum means minimum over the
finite tested set, not a proof of the worst case. Random bounded and separately
constructed adversarial signals are both included. The bounds-violated stress
is excluded from this table. Widths count each independent design once.

The pooled Monte Carlo proportion is descriptive. Noise draws are paired across
methods, and multiple signals share a design, so treating all rows as independent
binomial trials would give invalid comparison error bars. Full Gaussian
posteriors separately passed the exactly matched prior-predictive check at 0.95;
fixed-parameter boundary failures do not contradict that result.

The unregularized baseline is present in the CSV (including its very large
widths); the figure focuses on the six main methods to remain legible. The
incomplete run is not included. This is an exploratory, not confirmatory, study.
