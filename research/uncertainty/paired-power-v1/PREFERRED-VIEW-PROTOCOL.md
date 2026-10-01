# Fresh preferred-view controls for the conditional-noise test

Post-outcome follow-up of the completed viewing-variance experiment. Its four
ranged-amplitude scores, thresholds and both candidate calibrations are frozen.
This does not refit a score or establish a three-stack positive result. The
10076 zero scores and all fixed-amplitude scores remain in the prior inventory.

For each of 10028 and 10049, simulate all nine viewing laws obtained by
conditioning Haar rotations on |R[2,j]| >= 1-1/kappa, for j=0,1,2 and
kappa=1.1,2,5. Here k=plane @ R, so row 2 is the beam direction in object
coordinates. Under Haar, its coordinate is uniform on [-1,1]; the acceptance
probability is 1/kappa and the conditional density relative to Haar is exactly
kappa on the cap union. These are specified model distributions, not inferred
experimental viewing laws. Rejection sampling retains complete rotations and
therefore retains in-plane-angle randomness. No search chooses a favorable axis.

Use 65,536 fresh views per law, the same first known CTF/noise profile, and both
true and removed physical-cell maps. At each view, use amplitudes .9,1,1.1 with
the same fresh proper Gaussian noise, and both power and power+bispectrum
contrasts. Every map/amplitude/score cell is retained. The two stacks and nine
laws use independent seeds 261014 + int(dataset) + 100000*law_index. Batches
contain 512 retained views; the rejection sampler discards surplus accepted
rotations without inspecting any score. It records proposed, accepted and used
counts. Draw noise only after sampling a complete batch of rotations.

Save each law's scalar scores, first rotations/means/noise and RNG states. Verify
the first view of each map by an independent direct physical-cell sum, and
compare every first cubic evaluation with direct moment features at all three
amplitudes. Report diagnostics for the law's support, conditional coordinate
second moment and acceptance fraction. Do not treat Monte Carlo agreement
with these moments as a proof of the sampling distribution.

The primary validation outcome is the correctly specified null, at each law's
matching kappa, for both candidate maps and all amplitudes. Compare the same
three procedures with their archived calibration bounds and delta=.001. Compute
single-image event frequencies with pointwise exact 95% intervals, and project
binomial rejection probabilities at n=1,000/10,000/100,000 as before. Additionally,
partition each law's 65,536 independent views into 64 disjoint groups of 1,024
particles, recording the actual count of rejected groups and its pointwise
binomial interval. Calibration is fixed: these 64 groups estimate conditional
rejection probability, not unconditional error over repeated calibrations.
All methods/maps/amplitudes share draws, so cross-cell outcomes are dependent.
There is no simultaneous-error claim or method selection using test outcomes.

The true-map data against the removed candidate remain the specified alternative
controls. Removed-map data against the true candidate are also retained, without
presuming that a one-sided score detects this opposite error. The law's kappa
is supplied exactly; this study cannot calibrate it from experimental images.
The score still uses the oracle full-removal alternative, no smaller changes or
multiple CTF profiles are included, and the numerical computations are ordinary
floating point rather than certified interval arithmetic. This checks an
assumption-respecting preferred-view simulation, not experimental validity.
