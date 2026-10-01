# An exact unbiased bounded score cannot remove an arbitrary amplitude law

1 October 2026. A limitation of the proposed feature-response construction,
not a selected method or novelty claim. It strengthens the affine observation
in the population algebra note and applies even with the pose known.

**Proposition.** Fix two mean vectors m_A,m_B in R^d and known nonsingular
Gaussian noise. Suppose the amplitude may take every value in a nonempty
open interval for each state, with otherwise unrestricted state-specific
amplitude laws. There is no bounded measurable image score psi such that

    E[psi(Y) | state=A, amplitude=a] = 1,
    E[psi(Y) | state=B, amplitude=a] = 0

for all allowed amplitudes. The two intervals may differ, and neither needs
to contain zero. Means and noise can first be whitened without changing the
claim. The result is conditional on a fixed allowed pose, so observing that
pose cannot by itself avoid it.

**Proof.** For z=A,B define

    F_z(a) = integral psi(y) phi_d(y-a m_z) dy,

where phi_d is the standard Gaussian density. For bounded psi this function
extends to an entire function of complex a. Indeed for a=u+iv the absolute
value of the complex Gaussian integrand is

    |psi(y)| phi_d(y-u m_z) exp(v^2 ||m_z||^2 / 2).

On compact subsets of complex a it and all required derivatives have
integrable dominating Gaussian tails. Differentiation under the integral
therefore proves analyticity. If F_A=1 on an open real interval, the identity
theorem implies F_A=1 everywhere; similarly F_B=0 everywhere. At a=0 both
integrals equal E[psi(epsilon)], a contradiction. The value zero is used only
for analytic continuation, not asserted to be an allowed physical amplitude.

**Consequence for bounded features.** If h has finitely many bounded
coordinates and the fitted affine coefficients are finite, then
psi=lambda_0+lambda^T h is bounded. Exact zero response bias simultaneously
over both continuous amplitude intervals is impossible. A finite catalogue
with zero calibration residual does not overcome this fact: the residual
must appear somewhere between catalogue points unless assumptions change.
The proposed method already allows nonzero beta; this proposition explains
why that bias cannot simply be set to zero from catalogue feasibility.

**What this does not prove.** It does not give a quantitative lower bound on
the smallest approximate bias at a specified variance, rule out nonlinear
multi-image inference, prove non-identifiability, or forbid useful population
intervals. Disjoint noiseless mean supports can still identify the population
through Gaussian deconvolution with infinite data. Exact unbiasedness,
statistical identifiability and stable finite-sample estimation are different
properties. It does not establish a root-n impossibility theorem without a
separate regularity/tangent-space argument. Such a theorem is not claimed.

The general linear-functional deconvolution literature is relevant, including
[Pensky (2017)](https://sciences.ucf.edu/math/mpensky/wp-content/uploads/sites/9/2018/02/AOS2017-1498Published.pdf).
Only its primary abstract/introduction was screened here; its detailed rate
theorems have not been audited for this multivariate support-constrained model.
The elementary analytic proof above stands independently. No priority claim
is inferred from that targeted search, and no new numerical experiment or
post-failure region selection is authorized by this note.
