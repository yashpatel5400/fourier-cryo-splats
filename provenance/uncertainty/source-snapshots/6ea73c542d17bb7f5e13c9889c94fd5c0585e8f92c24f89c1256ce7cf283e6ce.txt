# Prospective bounded cubic-pose diagnostic

Declared 2026-09-30 after the remainder-only two-case probe. This is a
post-development diagnostic, not confirmation or an independently selected
benchmark. Its sole empirical case is the already completed fixed-weight
EMPIAR-10049 pilot_region_1 (128 particles, radius 12, 20 A Gaussian target)
at one degree and 0.5 A in the per-particle joint pose ball. The fixed weights,
pilot, target, density class B=2/P=1 and marginal error allocation are unchanged.
No inference images or reference outcomes choose the cubic operator or scales.

Retain every Taylor term through degree three, including rotation curvature
and the cubic term of the rotation exponential. For the five scaled pose
coordinates u, use normalized symmetric lifts
z_alpha=sqrt(m!/alpha!) u^alpha, |alpha|=m, m=1,2,3.
Each lift has squared norm ||u||^(2m). Thus any positive particle/degree
scales d give a combined lifted radius sqrt(sum d) for joint unit balls.
The field has 55 columns per particle and 20 spatial monomials. Compute scales
on order-12 quadrature only as a design choice, then audit order-80 quadrature
with a degree-six analytic product-kernel error bound. Use four fresh Gaussian
probes, 40 power steps, seed 640001 and failure probability 1e-6/12. Reusing a
successful random event after selecting weights or changing scales is forbidden.

For Phi(t)=Phi0+t Phi1+t^2 Phi2+t^3 Phi3+..., the retained exponential terms are
E1=i Phi1, E2=i Phi2-Phi1^2/2, E3=i Phi3-Phi1 Phi2-i Phi1^3/6.
Phi3=-(angle^2/6)||u_rot||^2 Phi1_rot follows from the skew-matrix identity.
Use the degree-three remainder from the Bell-polynomial/Sobolev diagnostic,
including its planar-embedding residual pad. Pair the pilot and target with
all twenty analytic continuous moments. Retain the nominal residual/pose cross
term using the existing joint-bias inequality; include its quadrature pad.
Ordinary floating point and NUFFT checks do not constitute interval arithmetic.

Before any empirical run: test the operator against an independent scalar
Taylor recurrence of the exact nonlinear phase, adjoint consistency, the old
quadratic operator when cubic lift coordinates vanish, NUFFT/direct agreement,
continuous cell and Gaussian moments against independent quadrature, and the
bound against sampled nonlinear perturbations on a toy continuous problem.
Commit source and this declaration before computing the empirical outcome.

Persist completion/failure, all scales, provenance, spectral history, cross
vectors, remainder components, runtime, peak memory, and every existing
reference scenario (nominal/coherent_x/random_boundary). Preserve the source
quadratic result and compare both at the same weights and allocations. Do not
choose the better random certificate without paying its failure probability:
the diagnostic cubic interval stands alone. Do not run additional empirical
cases or change iterations/scales based on this result under this protocol.
If utility remains poor, report that fact and diagnose the separate polynomial,
density, remainder and noise contributions before proposing more experiments.

## Conditional bias proposition and proof

Let F map the concatenated scaled lifts to the real cubic Fourier-adjoint
field on the continuous unit cube. Let d_im>0 and divide each degree-m,
particle-i column by sqrt(d_im). The multinomial theorem gives
sum_{|alpha|=m}(m!/alpha!)u^(2 alpha)=||u||^(2m)<=1.
Consequently the corresponding scaled lift v has ||v||<=L=sqrt(sum_im d_im).
An upper bound s^2 on the continuous spectral norm squared gives
sup_v ||Fv||<=f=L s. The degree-six product-kernel pad bounds the difference
between the quadrature Gram and the continuous Gram by epsilon||a||^2,
where a contains absolute polynomial-Fourier column coefficient sums after
scaling. Thus a valid quadrature spectral upper plus this pad supplies s^2.
Only the randomized spectral event consumes its declared delta; degree and
scales are fixed before those probes. This is a classical relaxation, not a
new concentration theorem.

Write h=ell-A0* w, ||h||<=H, ||F* h||<=c and ||F* rho0||<=p. For any allowable
lift, ||h-Fv||^2 <= H^2+2 L c+f^2 and
|<rho0,Fv>| <= min(P f,L p). The degree-three Taylor remainder e_u obeys
||e_u||<=r from the fourth-derivative Bell/Sobolev construction. For every
rho=rho0+g with ||g||<=B and ||rho0||<=P, the absolute noiseless estimation
error is therefore at most

    B sqrt(H^2+2 L c+f^2) + min(P f,L p) + (B+P) r.

The ordinary triangle expression B(H+f)+min(P f,L p)+(B+P)r is another bound
on the same event, so their deterministic minimum is allowed. The usual
bias-aware Gaussian interval then has marginal coverage at least
1-alpha_noise-delta under the stated independent Gaussian-noise model.
No union with a separately selected quadratic random event is claimed here.
No theorem asserts that the actual experimental nuisance radii are correct.
