# One proposed repair of the orientation-integration gate

1 October 2026 UTC. This protocol is written before catalogue-proposal outcomes.
Execution remains contingent on the focused methods consultation; the code is
prepared so that any recommendation is concrete. This is the one substantive
revision permitted after the [first failure](ADAPTIVE-POSE-INTEGRATION-RESULTS.md).
It is classical numerical integration, not a new uncertainty method.

## Measured obstacle and prediction

The first proposal used eight optimized modes and 5% Haar defense. It failed
the fixed tolerance on 10049/10076, with a median ESS below 256 on 10076.
Post hoc inspection of saved weights shows that large weights frequently come
from the defensive component far from the eight modes. This suggests inadequate
mode coverage; it does not prove that modes are the only error source.

The prediction is that broadening the proposal, without more samples or local
optimization, will pass the unchanged gate on all three stacks: at least 90%
of 128 images per stack must have two-bank log-ratio discrepancy at most .01,
and all four median ESS values must be at least 256. A failure ends this
numerical branch. No third proposal, relaxed threshold, extra draws, selected
images, new noise seeds or favorable subset is permitted under this protocol.

## Fixed proposal and computation

Reuse every observed image and all eight local centers, covariances, weights
and optimizer outcomes from attempt one, checked by hashes. Do not refit.
For the same 32,768 Haar training orientations, evaluate midpoint-map scores
using the existing exact training means. Give every catalogue orientation
weight proportional to exp(score/2). Temperature two broadens support relative
to the untempered midpoint likelihood; it is fixed before outcomes. Underflowed
weights may be zero, but the Haar component retains full support. There is no
top-k selection, per-image scale tuning or empirical truncation of the density.

At each catalogue orientation use an isotropic ACG quaternion kernel with
Sigma = diag(s²,s²,s²,1), s = (8 degrees in radians)/2. Eight degrees is a fixed
broad proposal scale, not a pose confidence radius. Retain the original local
ACG kernels without their embedded Haar component. The complete proposal is

    q(R|y) = .05 Haar + .45 local-ACG-mixture + .50 catalogue-ACG-mixture.

For unit center c and unit quaternion u, the isotropic component density
relative to normalized Haar is

    s^(-3) [ (1 - (cᵀu)²)/s² + (cᵀu)² ]^(-2).

Evaluate the full mixture density over all 32,768 centers in memory chunks.
The catalogue is used only to propose continuous rotations. It does not
replace the continuous likelihood with a discrete-support model. Fresh draws
are generated conditionally on the frozen image-specific proposal. The maps,
noise, CTF, frequencies, 64 source indices, two states, two 8,192-draw banks,
4,096-draw prefixes and seed formula are unchanged from attempt one. The new
sampling algorithm changes the realized proposals but does not select seeds.
Evaluate physical likelihoods with the same constant-cell NUFFT.

Retain every observation, quaternion, component/family label, catalogue weight,
log proposal, both physical log kernels, timing and all integration summaries.
Run all three stacks. Record the first attempt's fitting failures as inherited
facts; do not mislabel reused modes as newly converged optimizations.

## Verification and interpretation

Targeted tests compare densities with independently transformed full 4D
covariances, antipodal invariance, the Haar limit and inverse-density moments.
Replay every stored integral/ESS/gate. Check proposal densities independently
using full covariance matrices, and selected physical kernels with direct sums
over all 64³ cells. State exactly which samples each independent check covers.

Passing is only a necessary computational-feasibility result. Two banks can
miss the same region, and a per-image .01 tolerance does not control a product
likelihood over a large stack. No finite-sample coverage, known experimental
noise, uniform true viewing law, genuine conformation population, numerical
error certificate or new theorem follows. A synthetic timing check (uniform
random catalogue, no experimental outcomes) took about one second for one
8,192-point full-density evaluation on this Mac. No GPU rental is needed.
