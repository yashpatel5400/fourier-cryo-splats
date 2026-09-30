# Bounded CPU continuation of timed-out RELION fits

30 September 2026, declared after the first v2 initializer timed out and before
any continuation. This changes the compute allowance, not the scientific
settings or selection criterion. Original v2 records and checkpoints are
retained. The remaining v2 attempts continue under their original protocol.

For each of 10028, 10049 and 10076, wait for its v2 attempt to finish. Continue
only a wall-limited stage, from its highest numbered saved optimizer. A failed
software run or a successful full run is not silently retried by this rule.
Use a separate `relion-reconstruction-v3` output directory and record the full
parent record, checkpoint family and executable hashes before starting.

An interrupted initializer resumes to its original total of 100 VDAM batches,
with eight CPU threads and a four-hour additional wall limit. If successful,
perform the originally specified auto-refinement with three MPI ranks and four
threads per worker, fifty maximum iterations and a six-hour wall limit. A v2
refinement that itself hit its wall limit resumes its last saved optimizer to
the same fifty-iteration ceiling, with the same six-hour allowance. Retain the
original seed, images, half labels, angular/translation search, masking,
resolution limits, regularization and symmetry settings. Continuation inherits
the settings stored in RELION's optimizer file. Thread-dependent reductions
and checkpoint serialization mean this is not promised to be bitwise identical
to an uninterrupted run.

Process datasets sequentially, at most eight refinement-worker threads at once.
The maximum additional allowance is thirty CPU wall-clock hours for the full
batch, with early termination on a completed, failed or timed-out stage. This
uses the authorized Mac; it does not rent compute or depend on GPU availability.
Any time/iteration limit remains an unconverged outcome. No force-convergence
flag, reference-based checkpoint selection or altered restart is allowed.

Use the original evaluation rule: pilot-only frame/hand alignment, verified
exposure half labels, native half FSC plus approximate-reference/cross-method
agreement with all numerical limits explicit. The last available paired maps
may be reported as unconverged, never as a successful converged baseline.
Conventional reconstruction agreement and RELION's angular-accuracy estimates
are not experimental density-coverage or uniform pose-confidence guarantees.

Primary implementation checked: RELION tag 5.0.1, `MlOptimiser::read` and
`parseContinue` in `src/ml_optimiser.cpp`. The `--continue` path restores the
optimizer state; `--iter` denotes the total ceiling, not extra iterations.
https://github.com/3dem/relion/blob/5.0.1/src/ml_optimiser.cpp
