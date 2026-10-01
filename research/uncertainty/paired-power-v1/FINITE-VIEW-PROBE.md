# Declared finite-view oracle gate

1 October 2026 UTC, after saved-array R13/R14 diagnosis and before this probe.
The old pose bounds are abandoned as a route to practical widths. This gate
asks whether the distinct paired-power candidate has enough discrimination to
justify continuous-orientation work. It is neither a new coverage study nor
an experiment on independent acquired exposure pairs.

Reuse the exact Fourier means from the earlier discrete mixture screen on all
three stacks: 64 declared views, radius 12, the true map, a half-removed local
region, a fully removed region, and zero signal. Use CTF profiles at positions
0,16,...,112 of each fixed 128-particle cohort. All eight profiles are retained.
The alternative signal is the uniform-view mean power of the true map, supplied
as an oracle. No observed image or simulated noise is used to learn weights.
Real-coordinate covariance upper bounds are v in {1,2,4} times the generating
unit variance. The probe therefore has 3 x 8 x 4 x 3 = 288 optimization cases.

For each case maximize the expected log factor

    s' x / v + sum log(1-x_j^2)

subject to every sampled candidate view's power vector p satisfying
`p' [x/(1-x)] <= 0`. The dimensional weights are t=x/v. The constraint is convex
because powers are nonnegative. A finite-view constraint is weaker than the
continuous null. Oracle choice of weights and a sampled null are optimistic;
a successful gate is not a valid continuous test.

The primal numerical solve initially uses |x|<=.95, CLARABEL, 300 iterations,
1e-8 feasibility/gap tolerances. Preserve all statuses and failures. Repair
small numerical violations by moving the transformed coefficient uniformly
downward, and retain both raw and repaired values. Check a separable Lagrange
dual **over the full interval -1<x<1**, using nonnegative solver multipliers.
This distinguishes a numerical constrained lower value from an upper value
for the complete finite-view relaxation. Floating-point dual evaluation is
diagnostic, not validated arithmetic. A large gap triggers investigation,
not a claim that discrimination is impossible.

Negative controls: the true-map cone must have zero optimal expected growth;
global amplitude and complex conjugation leave that cone unchanged. The zero
map has an independent unconstrained scalar solution. Sum the eight expected
log factors with equal multiplicity 16 only as a declared 128-particle design
diagnostic. Exceeding log(20) in expectation is **not rejection power**. Report
each profile and each noise bound, including failures; make no Type-I estimate
without a later independent simulation and continuous-orbit certificate.

Stop this screen at 30 minutes overall, recording incomplete cases. No GPU or
paid compute is needed. Any next stage depends on whether local removal
retains useful expected growth with finite dual gaps, rather than only the
easy zero-signal control.

## Numerical runner retry

The first invocation at source `50c0beb` failed before saving any case, because
`variance_upper` was supplied twice while formatting the first solved record.
Its incomplete summary and traceback are retained in `paired-power-finite-view-v1`.
The first iterate was computed but was not written before that formatting
exception; no statistical result is inferred from it. The corrected retry uses
a new `paired-power-finite-view-v2` directory, removes only the duplicate key
and adds partial-array preservation. Scientific settings and stopping rules
are unchanged.

The v2 retry exposed a second serialization issue: a NumPy Boolean in the
solver metadata was not JSON serializable, including in the exception handler.
Its untouched initial summary, log and saved partial weight array remain. An
external failure manifest records this explicitly. The v3 retry converts that
flag to a Python Boolean and adds a JSON-serialization check to the numerical
test. No mathematical or experimental parameter changes.

The v3 run retained 139 completed cases before CLARABEL failed on 10049,
profile 48, full region removal, variance bound 2. The source mistakenly
terminated the entire sweep on a numerical solver exception, despite the
protocol requiring failed cases to remain. The v4 runner preserves a failed
case in its denominator and continues. It repeats the same 288-case design
and settings, retains every old record, and marks an aggregate upper bound
unavailable if any of its eight profiles fails; zero betting supplies the
feasible lower value for a failed profile. A convergence failure will not be
interpreted as evidence for or against discrimination.
