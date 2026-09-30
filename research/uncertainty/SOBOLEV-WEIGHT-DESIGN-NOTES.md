# Convex remainder calculations for weight design

Implementation groundwork, 2026-09-30; no new empirical fit or outcome is claimed.
The bounded cubic audit remains a fixed-weight diagnostic. These calculations
prepare a potential pose-aware optimizer if its separate polynomial, remainder
and noise terms justify further work. A new empirical protocol would be required.

On a symmetric enclosing ball or cube, write the known real field as
Re sum_q t_q (w_Rq+i w_Iq) exp(2 pi i k_q.x). Let K(v) be the domain's Fourier
kernel. Its order-j derivative-tensor squared norm has no real/imaginary cross
block. The two blocks are

    G_j^+ [q,r] = (2 pi)^(2j)/2 (k_q.k_r)^j t_q t_r
                  [K(k_q-k_r)+(-1)^j K(k_q+k_r)],
    G_j^- [q,r] = (2 pi)^(2j)/2 (k_q.k_r)^j t_q t_r
                  [K(k_q-k_r)-(-1)^j K(k_q+k_r)].

This follows by realifying the existing complex Fourier tensor identity;
odd integrands vanish on a symmetric domain. Both are Gram matrices of real
functions, hence positive semidefinite in exact arithmetic. The norm is
sqrt(w_R^T G_j^+ w_R+w_I^T G_j^- w_I). A positive Bell-coefficient sum of these
norms, plus the approximate-embedding amplitude pad, is convex and homogeneous
in the estimator weights. Its subgradient is the summed Gram action divided
by each nonzero norm, with zero chosen at zero. The implementation's diagonal
pads protect numerical convexity heuristically; they are not interval bounds.

Tests compare this independent real-block calculation with the existing complex
pair identity, verify the gradient by directional finite differences, check
Euler homogeneity and convex supporting inequalities, and retain both planar
and distorted geometries. They do not show that weight optimization will improve
usefulness. No original estimator or calibration model has been modified.

The separate `DifferentiableCubicPoseFieldOperator` also supplies the weight
adjoint of v^T F(w)u and of the known-pilot moment pairing. A future convex
surrogate can combine z||w|| + B||ell-A0*w|| + B L||F(w)|| +
L||F(w)*rho0|| + (B+P)r(w), with fixed positive scales and a fresh final
spectral audit. This triangle-based surrogate is distinct from the tighter
joint cross-term post-audit; no convexity is claimed for that post-audit formula.
Tests verify directional derivatives and the homogeneous bilinear identity.
No convergence certificate, fitted weights or performance improvement has yet
been produced by this proposed optimizer.
