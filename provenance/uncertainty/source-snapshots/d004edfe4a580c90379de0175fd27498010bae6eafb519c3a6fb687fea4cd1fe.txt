# Envelope-mixture-guided refinement: development protocol

30 September 2026. Declared after all continuous-mixture-v1, anchor and local
curvature outcomes, before running this follow-up. This is a computational
development experiment on the same observations, not fresh validation.

Reuse each of the three parent simulations, including its 8,224-box covering
partition, exact center kernels and feasible mixture witness. Do not use the
saved generating orientations, oracle support likelihood, signal, or local
curvature optimization outcomes as optimizer inputs. Known noise and zero
translation remain explicit restrictions. No predictive numerator or rejection
frequency is computed in this study.

Fit mixture weights to the cell envelopes, using matrix EM for at most 100
iterations per cover and checking the original-log-kernel scaled Lindsay gap
every 25 iterations against 0.01. Matrix EM may floor tiny terms at 1e-250;
certificate evaluation uses the original finite log kernels. Its convergence
is not required for the resulting real-arithmetic upper bound.

Refine batches of 64 leaves. For fitted weights w and row mixture z, prioritize
`w_c sum_i (U_ic/z_i)(1-p_ic(center)/U_ic)` plus
`max(sum_i U_ic/z_i-n,0)/number_of_cells` and a 1e-12 tie/exploration term
times the dual score. This is a heuristic computational choice, with no
optimal-complexity claim. Split each selected box along its widest Euler
coordinate; allocate half its guide weight to each child. Intersect every
child kernel envelope with its parent's envelope. This is valid because the
child is a subset of the parent, and makes the saved upper matrices monotone
under subdivision. Keep every leaf, including those with negligible weights.

Use the earlier first-order/coordinate enclosure for boxes with a sum of Euler
half-widths greater than five degrees. At or below that threshold, also use the
published curvature enclosure. This avoids spending its extra cost on large
boxes where the completed diagnostic found no benefit. No approximation of
Gaussian tails is introduced.

Every four batches, and at the final split budget, improve the finite feasible
lower bound using the union of the four highest evaluated center kernels per
image (at most 512 support points), 200 EM iterations and tolerance 0.01.
Retain the best previous feasible witness even when that restricted support
does worse. This support restriction never supplies a continuous upper bound.

Budgets are 32,768 additional splits, 900 seconds per geometry (checked between
completed batches), or continuous upper-minus-feasible gap at most one log
unit, whichever comes first. A batch and final serialization may extend wall
time beyond the nominal cap. Run all three geometries regardless of earlier
outcomes. Preserve all cover histories, bisections, evaluated centers, original
and best feasible witnesses, upper anchors, source/input hashes and failures.
Report runner and solver time separately. No outcome-dependent extension or
replacement of this family is permitted; further work requires a new protocol.

## Validity of the refinement

For each box C, let U_ic bound p_i(R) over C. Any shared orientation measure G
assigns masses w_C to the partition. Then integral p_i dG <= sum_C w_C U_ic.
Consequently the optimal finite envelope-mixture likelihood upper-bounds the
continuous likelihood. A scaled positive-anchor dual further upper-bounds
that finite optimum. We use the established Lindsay (1983) mixture dual, not
a new general likelihood theorem. The envelope-mixture primal is only a lower
bound on this relaxation's optimum, never on the continuous optimum.

Bisection preserves the initial complete cover. Taking the minimum of valid
parent and child enclosures remains valid on the child. The pointwise guide
and its floors only determine which cells to split. Retaining every cell and
re-evaluating the dual in original log units preserves the real-arithmetic
upper-bound argument irrespective of priority, EM convergence or stopping.
The finite center mixture independently gives a feasible continuous lower.
These facts do not certify floating-point rounding; no interval-arithmetic or
experimental calibration claim is made. Unit tests compare replayed witnesses,
partition volume/membership, inherited bounds, a dense independent orientation
fit and an exact zero-map control before any recorded geometry calculation.
