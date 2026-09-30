# Continuous support and outside-density sensitivity

This is a post-review development experiment for R2, not an empirically
calibrated support estimate or a new frozen confirmation. It keeps the unknown
density continuous. The following coordinate identity and triangle bounds are
elementary consequences of the existing dual audit; no novelty is claimed for
the change of variables itself.

Let S_l=[-l/2,l/2]^3, 0<l<=1, inside the original unit cube D. For a density
supported on S_l define f(u)=l^(3/2) rho(l u). This is an L2 isometry. Its
Fourier observations are obtained from the original implementation with

    k' = l k, q' = l q, C' = l^(3/2) C, shift' = shift/l.

For a normalized Gaussian target of width sigma and center c, the transformed
target has width sigma/l, center c/l and coefficient l^(-3/2). Hence both the
physical target and every nonlinear rotated/translated observation are
unchanged. Rotation commutes with the isotropic coordinate rescaling. CTF
values are evaluated at the original physical frequencies before rescaling.
The supplied noise level is unchanged.

For rho = p_S + d_S + t, suppose p_S and d_S are supported on S_l,
||d_S||_2 <= B, t is supported on D\S_l, and ||t||_2 <= epsilon. The estimator
is <ell,p_S> + w^T(y-A_0 p_S). An inside-support bias certificate applies to
p_S+d_S. The remaining error is <A_u^*w-ell,t>, bounded by epsilon times the
full-cube perturbed residual norm. A convenient conservative upper bound is

    epsilon [ ||A_0^*w-ell||_2 + sum_iq |c_iq| d_iq ],
    d_iq = min(2, 2*pi*sqrt(theta^2*||k_iq||^2/12
                                 + shift^2*||q_iq||^2)).

To prove it, use |exp(ia)-exp(ib)|<=min(2,|a-b|). For uniform integration on
D, the linear spatial phase has zero mean and variance ||Delta k||^2/12,
and the constant shift phase has squared magnitude at most
shift^2*||q||^2. A rotation of size at most theta has
||Delta k||<=theta||k||. Apply the L2 triangle inequality to the Fourier sum,
then restrict to D\S_l. This loses cancellation and uses the full residual
for the outside part; it is an upper bound, not a sharp support complement
calculation. These statements concern real arithmetic; the inherited numerical
quadrature and floating-point qualifications still apply.

The zero-weight alternative has half-width at most
B||ell restricted to S_l|| + epsilon||ell||. Its center must also revert to
<ell,p_S>; retaining the fitted center while replacing only the width is invalid.

## Declared development probe

`scripts/probe_uq_continuous_support.py` starts with 128 existing development
particles from EMPIAR-10049, radius 12, a centered Gaussian target with physical
sigma 10 Angstrom, 1 degree rotation and 0.5 Angstrom translation. Cube side
fractions are 0.5, 0.75 and 1.0. The original radius-5 simulation noise scale
stays fixed across these designs. Estimator weights minimize the fixed-pose
continuous objective, then receive a matrix-free pose audit with fresh probes;
this is not the in-progress pose-aware weight optimizer. The outside density
radii are 0, 0.01, 0.05 and 0.1. These are assumptions, not estimated error bars.

Both the original deposited-map generator and its exactly cropped version are
evaluated, without renormalizing after cropping. The former is labeled outside
the class whenever its omitted L2 norm exceeds the supplied allowance. The
latter measures a controlled support ablation and is not the original molecule.
Reference values are never used to choose weights, support or interval widths.
The inside class has B=2 and the pilot's actual restricted norm; triangle
inequality verifies inside membership for these unit-norm reference generators.

Before any fitting outcome was available, the original map generators showed
outside L2 norms of 0.285/0.079 (10028), 0.320/0.145 (10049), and 0.429/0.287
(10076) for side fractions 0.5/0.75. Thus the strongest hard-mask assumptions
exclude the original maps; their apparent precision cannot be reported as
validated molecular uncertainty. The observed tail norms are diagnostics of
deposited, processed approximate references, not externally calibrated bounds
for the unknown experimental density.

The implementation tests compare exact cell signals, target pairings and
nonlinear poses before and after the isometry. A separate direct Fourier-field
calculation checks the tail inequality for nonlinear perturbations. Two tests
passed; the development fitting study is still running at this checkpoint.

## Completed three-size sweep

All three 10049 side fractions completed. The full-cube control required
6,324.04 seconds and reached a nominal sum-objective gap of 0.001867, but
every tested tail allowance fell back to no data after pose auditing. Its
integrated cubic bias alone was 115.36. This is a negative conditional result,
not an impossibility theorem and not evidence that the outside-density
allowances are empirically calibrated. Exact completed weights are included
in the separate v0.4 development release.
