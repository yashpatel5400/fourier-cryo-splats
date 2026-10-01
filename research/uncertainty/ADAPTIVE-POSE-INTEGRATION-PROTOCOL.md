# Numerical gate: image-specific pose importance sampling

1 October 2026 UTC, before any adaptive integration outcomes. The matched
Haar comparison shows median ESS 1.23 / 7.75 / 28.7 at 32,768 points and
material between-bank log-ratio variation. Address this numerical obstacle
before choosing a statistical uncertainty method. This is a classical
importance-sampling implementation, not a new UQ contribution or another
moment score. Likelihood integration and importance sampling are established
cryo-EM methods; see [Brubaker et al., CVPR 2015](https://www.cs.toronto.edu/~mbrubake/projects/CVPR15.pdf)
and the [BioEM analysis in cryo-BIFE](https://pmc.ncbi.nlm.nih.gov/articles/PMC8249403/).
The Brubaker paper's importance-sampling experiment was retrieved; the complete
paper and author implementation have not yet been read or executed for this gate.

Use the unchanged three candidate/region cell densities, frequencies, transfer,
unit independent white noise and Haar law. Select 64 indices by integer
linspace from 0 to 8191 in the existing matched replay. For each index evaluate
both its unchanged and 25%-deleted image: 128 observed images per stack.
Generating poses may be used only for reporting, never proposal construction.

For each observed image fit a proposal to the midpoint map (12.5% deletion).
Start from eight high-likelihood orientations from the first 32,768 saved
training draws, greedily separated by at least ten degrees. Optimize each
start using a 256-padded cubic Fourier interpolant, BFGS, at most 120 iterations.
Retain all convergence flags and improvements. A worse terminal point reverts
to its start and remains labelled; this is not a global alignment certificate.
Use a finite-difference 3D rotation Hessian with step .001 radian for proposal
shape. Clip its inverse-curvature rotation standard deviations to [.25,30]
degrees and record the original eigenvalues. These are proposal scales, not
pose uncertainty estimates. Keep duplicate converged modes instead of silently
discarding a failed or redundant start. Assign mode weights proportional to
the local Laplace mass (likelihood times square root covariance determinant).

The proposal is 5% Haar plus 95% of the eight rotated angular central Gaussian
(ACG) components on unit quaternions. For local covariance
Sigma=diag(C/4,1), with C the clipped rotation covariance, its density relative
to normalized Haar measure is |Sigma|^(-1/2)(q^T Sigma^(-1)q)^(-2).
This follows directly by radially integrating a centered four-dimensional
Gaussian. Antipodal symmetry makes it a density on SO(3). Gaussian sampling
followed by normalization exactly samples this proposal in real arithmetic.
The defensive Haar component guarantees support but does not guarantee useful
importance variance. Source: Tyler, Biometrika 74(3), 579–589 (1987),
DOI 10.1093/biomet/74.3.579; publisher access failed, so no claim of reading
that full article is made. Independently check the quaternion density against
the full transformed covariance and check isotropic and inverse-density
integrals in numerical tests.

Freeze each image's fitted proposal before its two independent importance
banks. Use 8,192 draws per bank with nested prefix 4,096. Seeds are
261020 + dataset + image_position*100 + state*10 + bank. Evaluate both map
likelihoods with the physical constant-cell NUFFT, not the proposal interpolant.
Save every sampled quaternion, proposal log density, two physical log kernels,
mode/Hessian/optimizer records and resulting log integrals, ratio, ESS and
maximum normalized weight. Model-dependent normalizing constants cancel only
within the same image's ratio. Check selected samples with direct cell sums
and replay importance weights independently.

**Prediction/gate:** on each of the three stacks, at least 90% of the 128
observed images have absolute between-bank log-ratio difference at most .01
at 8,192 draws, and the median ESS for each of the two models in each bank is
at least 256. Report all four ESS medians, density-integral disagreement,
prefix changes and every failed fit. Passing only a favorable stack is not a
pass. A failure permits one substantive revision directed at the measured
failure, not repeated seeds or a particle-count search. Two missed predictions
end this numerical branch under the round-4 stopping rule.

Passing would support computational feasibility under the declared simulator.
It would not certify quadrature error, prove correct orientation recovery,
establish novelty, justify Haar/noise assumptions for real particles, or give
coverage/size for a future statistical procedure. No external GPU spend.
