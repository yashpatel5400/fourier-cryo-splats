# Coordinate follow-up: completed without improvement

30 September 2026. The original cubic fit's 0.940 relative surrogate gap
triggered the already declared coordinate follow-up. It starts from the same
fixed-pose weights and retains the same target, class, noise and one-degree /
half-angstrom nuisance budgets. This is post-outcome development, not a fresh
calibration study. The original result is not overwritten.

The run completed in 8,800.66 seconds: 7,220.44 seconds of optimization and
1,550.15 seconds for the independent final spectral audit. Peak RSS was
3,347,185,664 bytes. It reached the thirty-iteration limit, selected evaluation
34, and did not converge. All evaluations and retained checkpoints remain.
The guide objective decreased from 3.04493 to 2.81584, but the independently
audited surrogate upper is 2.87586, its dual lower is 0.0673149, and its relative
gap is 0.976593. These do not establish near-optimality.

The joint refined half-width is 2.71131 against a no-data half-width of
12.1587: relative width 0.222993. Minimum analytic reference sign power across
the same three nonlinear patterns is 7.73395e-7. All three prescribed
reference-coverage checks equal one, reflecting conservatism rather than
experimental calibration. The result is worse than the original fit's
relative width 0.179785 and minimum sign power 0.00653666. Under these compute
budgets, the coordinate change did not improve the scientific outcome.

The separately declared thirteen-column convex design uses the original fit,
not this follow-up, to define its span. Its queued startup initially encountered
the provenance guard while a new, unrelated likelihood module was uncommitted.
No subspace output directory or scientific fit had been created. That failed
startup log is archived. After committing the numerically checked module, the
identical subspace source, protocol, seeds, inputs and command were restarted.
This is a provenance-guard retry, not an outcome-selected scientific rerun.
