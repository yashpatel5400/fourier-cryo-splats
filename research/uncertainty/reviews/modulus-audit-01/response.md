# Response to the two-pose modulus audit

This focused report confirms the dual signs, density construction, common
scaling, testing threshold and phase rotation. It is not full ICML review 2.
The unconditional experimental-calibration and usefulness objections remain.

T1: negative padded residual squares now raise, with finite-input/output
checks. A zero exact residual uses zero amplitude rather than division by a
tiny number. A regression deliberately supplies an inconsistent residual and
checks failure. No such invalid square occurred in the complete v2 grid.

T2: an order-3 quadrature test compares residual norms, Gram actions and both
scaled/unscaled witness distances with the exact continuous sinc Gram. It
exercises nonzero integration pads. An explicit real rotation matrix tests
the wrapped target and preconditioner as well as the Fourier Gram.

T3: the assembled pilot mean is now compared with `pose_cell_forward` on a
nontrivial joint rotation/translation instance.

T4: best-upper and best-lower weights and iteration indices are saved
separately; the pilot-mean distance and weight norm are recorded; the JSON
hashes its NPZ companion. A fixed independent Hilbert example with deliberately
inexact inner solves has the best lower at iteration 1 and upper at iteration
6. The test re-evaluates both saved witnesses and matches their reported bounds.

T5: feature gaps now receive a downward floating-point magnitude pad. Mean
norms and dual objectives receive analogous upward pads. These are disclosed
heuristic guards, not validated interval arithmetic or special-function bounds.

T6: the JSON explicitly records the no-data half-width B||ell|| and the
unnormalized constructive half-width. This is the same denominator used by
`scripts/run_uq_pose_optimized_study.py` with B=2. Any comparison must also match
the geometry, target, noise and nuisance class; a future summary will check
these fields rather than infer equivalence from a plot label.

T7: both pipelines call the same row-vector `perturbed_geometry` convention
and use physical shift/field scaling. The lower-grid script now records and
asserts the maximum joint five-vector norm. The existing upper audit's joint
ball is explicit in THEORY.md and in `scaled_pose_radius`; its larger product
variant is distinct. The complete lower-grid snapshots retain those sources.

T8: unchanged limitation. The testing formula is conditional on the declared
white Gaussian law; it is not an experimental whitening guarantee.

T9: the failed v1 record has a new companion FAILURE-ANNOTATION.md explaining
that its postmortem error/log fields were added by the agent after the process
crashed. The record is not retroactively presented as a successfully saved fit.

Seven targeted tests pass in 0.74 seconds and the full suite passes 102 tests
in 7.78 seconds. The unchanged 30-case protocol is rerunning in v3 to save the
additional witnesses and numerical guards. The complete v2 outcomes remain
available. The reviewer could not execute the code; these execution results
are the project's independent checks, not actions attributed to the reviewer.

Completion addendum: the unchanged v3 grid subsequently completed all thirty
cases, with maximum fixed-pair relative gap 0.0048301131. The v2 records remain
unchanged. The v3 summary and every witness are archived in the development
results; this completion does not resolve T8 or imply a full-paper acceptance.
