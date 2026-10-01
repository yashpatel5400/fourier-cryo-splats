# Larger-orientation cone diagnostic

1 October 2026 UTC. Post-outcome follow-up to the complete 288-case finite-view
gate. The half/full local-removal candidates have weak oracle expected growth
at 128 particles. Before studying continuous orientation certificates, test
whether that separation disappears when the candidate includes more views.

For each of the same three original maps and the same half/full local removal,
append nested batches of 256, 1,024 and 4,096 independent Haar rotations to the
original 64-view catalog. The seed is 261003 plus dataset ID. Keep the original
true map's uniform-64-view mean power as the alternative, so the scientific
question changes only by enlarging the null. Retain the true-map control.
Verify recomputed original Fourier templates against their saved arrays before
using a new orientation. Known CTF diagonal factors are applied only when
translating a cone approximation into a bound on the earlier eight profiles.

Retain a zero-added-view baseline as well. For each of 36 dataset/candidate/catalog cells, solve the linear feasibility
problem P' lambda = s, lambda >= 0, where rows of P are candidate projection
powers and s is the alternative mean power. Normalize equations by max(s_j,
1e-6 max(s)) for conditioning. Use HiGHS with 1e-9 primal/dual feasibility
tolerances. If equality is numerically infeasible, minimize the maximum
scaled equation residual with two linear inequalities per frequency. Retain
both solve statuses, nonnegative mixture weights, residuals and input arrays.

A nonnegative vector need not solve the equality to be useful: it is a valid
Lagrange multiplier for the power-cone constraint. Evaluate its separable
full-domain dual for every earlier CTF profile and v in {1,2,4}. This upper
value concerns the oracle expected log factor, not rejection probability.
The dual is evaluated numerically, not with validated interval arithmetic.
Any error remains visible; do not replace a failed cell with a zero bound.
Stop at 30 minutes, retaining partial arrays and all completed cells.

## Elementary interpretation

If s = sum_l lambda_l p_l with lambda_l >= 0, any admissible t obeys
sum_j t_j s_j <= sum_j [t_j/(1-t_j v_j)] s_j <= 0, since
t/(1-tv)-t=t^2 v/(1-tv)>=0. The logarithmic noise normalizer is nonpositive.
Hence no member of this linear cross-power betting family has positive
expected log growth for this alternative. The null may use probabilities
lambda/sum(lambda) and amplitude sqrt(sum(lambda)) to match its mean power.
This does not make the full image distributions identical, rule out other
statistics, or prove failure of every reconstruction/validation method.

A feasible finite mixture is also a mixture of continuous orientations, so
such a witness can diagnose a limitation without bounding all rotations.
An approximate numerical witness is reported with its residual and dual
upper value rather than asserted to be an exact algebraic equality.
