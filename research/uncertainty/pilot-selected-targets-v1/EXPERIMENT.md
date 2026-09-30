# Locked development experiment, 29 September 2026

This specification is written after pilot-only target selection and before any
reference outcome for the new targets. Selection lock SHA-256:
`57122bb2ee0d545584d854c1e8845adeced7620d4926a923fa9bf2c74af3e834`.
Keep all three selected regions and the matched-center control on all three
stacks. The 12 features form one multiplicity family for each declared pose
class. Use alpha_total = 0.05/12 and delta_spectral = 0.000001/12 per feature;
allocate their difference to Gaussian noise. Bonferroni then gives 95% family
coverage conditional on the entire independent pilot/design and the common
declared density, noise and pose assumptions. This does not calibrate those
assumptions. Different pose classes are sensitivity analyses, not 36 additional
independently chosen discoveries. Reusing the spectral event for deterministic
remainder refinements does not spend additional probability.

The fixed-pose baseline uses exact supplied poses and shifts. The three pose
audits use rotation radii 0, 1 and 2 degrees, all with 0.5 Å translation radius,
and the per-particle joint scaled five-dimensional unit ball. In particular,
the zero-degree audit is a translation-only class, not the fixed-pose baseline.
The global coordinate frame remains fixed, so coherent rotations remain in the
class. Width is Gaussian standard deviation 20 Å on each stack; Fourier radius
12 and 128 particles use seed 609315 and the existing inference_half0 split.
Use the archived radius-five prescribed noise scale, not an empirical estimate.
The class is the full unit-cube L2 ball, B=2 around the unit-L2 24-cell pilot.

Fit fixed-pose weights with order-80 quadrature, rank-1024 preconditioner,
relative sum-objective gap tolerance 0.005 and at most 100 outer iterations.
Retain incomplete optimization with its gap; it does not invalidate a correctly
audited feasible estimator. For poses use order 80, design order 12, four
independent Gaussian probes and 40 power steps. Certificate seeds are
620000 + 100*dataset_index + 10*feature_index + angle_index, with indices from
the lock and angle order [0,1,2]. Report the cross-term, sharp-cubic, known-pilot
bound and its deterministic per-particle minimum with the enclosing-cube
remainder. Do not choose between repeated spectral draws.

First time 10049/pilot_region_1 as a development compute pilot. Then execute all
12 fixed fits and 36 pose audits, preserving every outcome. Reference checks use
the unit-L2 deposited 64-cell map, nominal poses, coherent x rotation and random
joint-ball boundary poses; use the existing fixed seed 610281+dataset for the
random scenario across all targets and severities. Report bias, conditional
Gaussian coverage, sign power, width/no-data and fallback for every feature and
class. The reference is an approximate, processed comparator, not truth for
experimental density coverage. Report all failures rather than deleting them.

This is an exploratory utility test of a target-selection hypothesis. It is
not a frozen confirmatory study and does not resolve the reviewer's request
for calibrated experimental pose/noise/density inputs or a 10 Å application.
