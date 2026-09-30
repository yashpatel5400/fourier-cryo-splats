# Local curvature-bound diagnostic on the unchanged continuous simulations

30 September 2026. Declared after the large first continuous likelihood gap,
before these local curvature outcomes. This is post-outcome development.
The parent global computation and anchor study continue unchanged.

Reuse all three completed continuous-mixture-v1 simulations and their exact
maps/noise/CTFs. On each stack use particles at array positions 0, 64 and 127.
Center a local Euler box at that particle's true simulated orientation. Set
each of its three half-widths to one third of an angular budget in
{0.1, 0.5, 1, 2, 5, 10} degrees, so the sum bounds its geodesic displacement.
This uses the truth only to diagnose bounds; it is not a pose estimator or
a partition of all SO(3).

For every one of the 54 cases, compare the original envelope with the
second-order log-kernel envelope and its retained pointwise minimum. Run a
local likelihood maximization from the center and all eight corners, using
analytic gradients, L-BFGS-B, at most 100 iterations per start, ftol=1e-12 and
gtol=1e-9. Keep every start's result and status. These searches give feasible
lower values, not a guaranteed global box maximum. Report upper-to-feasible
gaps and the third-order remainder, not unqualified true envelope errors.

Save every bound, center likelihood, half-width, local optimum candidate,
input/source hash, and stopping status. Allow thirty minutes overall for
case calculations, checked between cases. Wait at most one hour per parent
stack before computation. Failed/incomplete runs remain recorded. Do not
change observations, drop cases or seed a new global optimizer from these
oracle-local results without a separate declared protocol.

An improvement here is only a prerequisite for a new global calculation.
It does not establish structural-test power, a learned predictive numerator,
or calibration of an experimental uncertainty statement.
