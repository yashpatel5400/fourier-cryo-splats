# Separate anchor diagnostic for the first continuous cover

30 September 2026. Declared after observing large intermediate upper bounds
in continuous-mixture-v1, while that unchanged run continues. This is explicit
post-outcome algorithm development. Do not alter its source, observations,
partition, or recorded outcomes.

After each stack completes, reuse its exact best and final Euler covering
partitions. Each column of their kernel-envelope matrix is already a valid
continuous cell upper bound under the specified Gaussian model. Fit the
finite envelope mixture with at most 500 EM iterations, gap target 0.01.
The feasible objective in this *envelope* problem is not a lower bound on
the true continuous likelihood. Label it only as an envelope-relaxation
primal value.

Compute Lindsay's scaled-anchor upper bound from each fitted mixture's primal
and upper-anchor weights, and from the per-image row maxima. The latter
anchor is arbitrary but admissible, and guarantees an upper no larger than
the sum of row maxima. Also retain that independent-image upper and the
absolute Gaussian maximum (normalization only). The minimum of valid upper
bounds remains valid pointwise; this is not statistical model selection.
Report the gap to the original finite-center feasible likelihood, and to
the separately recorded oracle-support feasible likelihood, without calling
either the true continuous maximum. Record the envelope fit's own gap.

Wait at most one hour per stack for its existing parent computation. Run all
three stacks and both cover choices, preserving every intermediate
and final outcome. Preserve anchors and weights with hashes of their input
cover arrays. Check that every proposed bound contains both feasible lower
values. No new observations, e-values or rejection counts are generated.
The purpose is to distinguish anchor looseness from cell-envelope looseness
before deciding whether further continuous refinement is worth its cost.
This diagnostic cannot establish a tight bound if its envelope primal itself
remains high above every available continuous feasible likelihood.
