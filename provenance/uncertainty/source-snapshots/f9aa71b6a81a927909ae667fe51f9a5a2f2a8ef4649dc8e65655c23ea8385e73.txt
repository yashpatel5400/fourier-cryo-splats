# Exact zero-rotation transform reduction

Declared 2026-09-30 before any shift-only audit in this locked feature family.
The early one-degree audits are complete and unchanged. This amendment changes
only implementation cost, not targets, particles, weights, nuisance radii,
quadrature, scales, seed allocations, spectral dimension or stopping rules.

When all rotation radii are zero, the first/second pose polynomials depend on
translation alone and are constant in the spatial coordinates. Only the
constant channel among the operator's ten NUFFT channels has a nonzero
coefficient. The new subclass computes that channel and returns exact zeros
for the unused channels. In the adjoint, the other spatial moments are unused
because their multiplying coefficients are identically zero. For every
nonzero rotation radius the original transforms run unchanged. The original
operator and all frozen study implementations remain unmodified.

Four direct/NUFFT tests compare fields, adjoints, weight derivatives and paired
inner products with nonuniform scales, at both zero and nonzero rotation. A
single actual 128-particle/radius-12/order-80 operation check on the first 10049
feature gives zero measured output difference in both directions. Combined
single-call time decreases from 11.85 to 1.40 seconds (8.47-fold) on one thread
under concurrent load. This is not a pipeline speed or steady-state benchmark.
The exact record is
`results/uncertainty/development/audit-regressions/shift-only-operation-equivalence.json`.

The same future audit runner dispatches to this subclass; its source snapshots
include the amendment and implementation. No completed result or random
certificate is recomputed, selected or replaced. No inference pixels are read
by the numerical equivalence check.
