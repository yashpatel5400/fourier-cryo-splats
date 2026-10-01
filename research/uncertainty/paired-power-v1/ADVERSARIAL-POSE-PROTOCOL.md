# Adversarial continuous-pose falsification on the surviving stack

The shifted-covariance screen has completed all twelve cells. Only 10049
retains positive achieved growth at all declared shift SDs after repair on
10,000 fresh rotation/shift draws. A sampled maximum still is not a continuous
bound. Before any calibration or large noise trials, search directly for
counterexamples to each of these four saved directions.

For each shift SD 0,.5,1,2, use its saved union-repaired direction. Start at
the twelve largest-score points from the preceding 10,000-view check, plus
eight independent Haar rotations with uniform shifts on [-32,32)^2 (seed
261009). Optimize five coordinates: a local rotation vector and two common
translations, with the translations confined to the full periodic box.
The null here permits unknown continuous translation even in the zero-shift
alternative arm. L-BFGS-B gets 150 iterations and 2,000 function evaluations
per start. Retain all 80 statuses, initial values and selected poses.

The guide uses a padded 256-cube FFT of the known 64-cube constant-cell
density, cubic interpolation, the analytic cell sinc factor and half-cell
phase. It is only a fast objective for finding possible violations. Independently
evaluate every selected pose with the original continuous-cell NUFFT forward
operator, retaining both scores and their discrepancy. A positive exact-forward
score is a numerical counterexample; a failure to find one is never a
continuous certificate. Save initial points as well when a failed optimizer
moves to a worse guide score.

For each direction, subtract the largest *verified* normalized violation
times identity, plus the existing numerical margin, and recompute the
directional growth on its exact Gaussian-shift alternative second moment.
This remains a sampled repair, even after optimization. The Fourier grid is
reconstructible scratch data; saved poses and exact values are the scientific
outputs. No observation pixels, new random noise trials or selection based
on coverage enter this falsification screen.
