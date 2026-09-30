# Cubic weight-design surrogate and its lower certificate

Development derivation, 30 September 2026. This specializes standard norm
majorization and weak duality; it is not a new general coverage theorem. The
protocol and successful prerequisite checks were committed in 2123d33 before
the empirical run. The empirical result is pending.

## Setting

Work in the real Hilbert space H=L2([-1/2,1/2]^3). A0 is the nominal realified
Fourier observation operator with declared noise whitening. The known pilot is
rho0, ||rho0||<=P, and rho=rho0+h with ||h||<=B. The fixed target is ell in H.
For every allowed per-particle joint pose, suppose the adjoint difference is
F_cont(w)v+r(w,v), where ||v||<=L and ||r(w,v)||<=r_PSD(w). The estimator is
<ell,rho0>+w^T(y-A0 rho0). Its bias magnitude is bounded by

  B||ell-A0*w|| + B L||F_cont(w)||op
    + L||F_cont(w)*rho0|| + (B+P)r_PSD(w).

This follows by expanding the error into terms paired with h and rho0. The
known-pilot polynomial term is evaluated directly, not replaced by its norm
bound. The remainder still pairs with both h and rho0. A numerical upper of the
polynomial norm may use fixed positive scales and quadrature; its probability
budget must be counted separately from observation noise.

## Proposition: a convex, continuous-valid design majorant

Let G=A0 A0*, let G_Q be the quadrature Gram, and a=A0 ell. Suppose its kernel
error E gives |w^T(G_Q-G)w|<=E(sum_j |t_j| b_j(w))^2, with b_j the Euclidean norm
of each real/imaginary weight pair. Set eta=E sum_j t_j^2. Then
G_Q+eta I-G is positive semidefinite, by Cauchy--Schwarz and sum b_j^2=||w||^2.
Consequently

  H(w)^2 = ||ell||^2 - 2 a^T w + w^T(G_Q+eta I)w

is a squared affine Hilbert norm and H(w)>=||ell-A0*w||. In particular it is
convex; its derivative at H>0 is (G_Q w+eta w-a)/H.

For the cubic field, suppose the degree-six quadrature error gives
||F_cont(w)||op^2 <= ||F_Q(w)||op^2 + E6||M b(w)||^2,
where M is the nonnegative matrix of absolute scaled column-coefficient sums.
The function

  J(w) = B L sqrt(||F_Q(w)||op^2+E6||M b(w)||^2)
          + L||F_cont(w)*rho0|| + (B+P)r_PSD(w)

is convex and positively homogeneous. To see this, each b_j is convex and
nonnegative, M preserves the componentwise ordering, and Euclidean norm is
convex and componentwise nondecreasing on the nonnegative orthant. The spectral
norm is convex; the square root is the Euclidean norm of the nonnegative pair
(||F_Q||op,sqrt(E6)||M b||). The pilot pairing is a linear function of w
followed by a Euclidean norm. Each remainder term is a positive multiple of a
PSD quadratic-form norm or complex-pair amplitude. Their sum has the same
properties. Thus z||w||+B H(w)+J(w) is a convex majorant of the sum-width
objective. It is not the exact folded-normal width or the cross-term refinement.

The derivative-Gram implementation adds reported nonnegative diagonal pads
before evaluating its norms. These preserve convexity of the computed PSD
objective. Their size is a heuristic floating-point safeguard, not a validated
operator-error envelope; all real-arithmetic statements above are conditional
on valid integration bounds and numerical evaluation.

## Proposition: arbitrary approximate modes still give a valid dual support

Let v_j be any unit spatial vectors. Put u_j=F_Q(w)*v_j/||F_Q(w)*v_j|| for
nonzero adjoints, and zero otherwise. Then the gradient g_j of
v_j^T F_Q(w)u_j satisfies g_j^T w' <= ||F_Q(w')||op for all w'. A convex mixture
g_F has the same support property, even if the directions are inaccurate Ritz
vectors. Let s be their log-sum-exp norm surrogate, p=sqrt(E6) M b(w), and
D=sqrt(s^2+||p||^2). Combine (s/D)g_F with the chain-rule amplitude support from
(p/D)^T sqrt(E6) M b(w'). Cauchy--Schwarz gives a support of the unsmoothed
joint spectral/integration norm because (s/D)^2+||p/D||^2=1. At D=0 use the zero
support. Add norm supports of the pilot and remainder terms to obtain g_J with

  g_J^T w' <= J(w') for every w'.

For positive smoothing, entropy makes g_J^T w generally smaller than the
reported smoothed value. Euler equality is not asserted. Approximate directions
and their smoothed value guide optimization only; neither is an upper bound on
the full spectral norm. Partial eigensolver convergence is retained explicitly.

## Proposition: a continuous weak-duality lower bound

At the chosen candidate w, take an upper Hbar>=||ell-A0*w||, set
f=B/Hbar and v=f(ell-A0*w). Its norm is at most B. Suppose e_G bounds
||(G_Q-G)w||. Let

  d = ||f(a-G_Q w)-g_J|| + f e_G,
  c = min(1,z/d),
  lower = max(0,c f (||ell||^2-a^T w)).

With the natural d=0 convention, ||A0(c v)-c g_J||<=z. For every w', Hilbert
norm duality and the support property imply

  z||w'|| + B||ell-A0*w'|| + J(w') >= <c v,ell>.

Replacing the continuous residual by its larger H(w') preserves this lower
bound. The zero lower bound is also valid. Notice the absence of eta*w in the
approximation a-G_Q w: the dual uses the actual continuous residual, not a
fictitious observation operator that includes the majorant's padding.

An independently drawn spectral upper gives the primal upper at the SAME pose
quadrature order. Their relative difference bounds optimization error for the
unsmoothed sum surrogate, subject to the stated numerical qualifications.
Replacing the primal by a narrower folded-normal or joint-cross interval would
change the objective and cannot inherit that reported gap.

## Statistical scope and independent numerical checks

Weights, targets, scales and fallback rules must be fixed independently of the
inference noise under the conditional model. The final Gaussian spectral probes
are drawn after selecting weights and have their own failure budget. Neither
condition proves experimental pose calibration, a density-class membership,
common noise covariance, homogeneity or independence of deposited poses.

The four-coordinate conic check gives optimum 2.071569215466, primal upper
2.071569215459 and lower 2.071569195636, within ordinary floating-point accuracy.
Separate finite differences, norm supports, a provable zero-weight optimum and
partial-Ritz tests are recorded. These establish implementation agreement at
the tested geometries and do not establish experimental usefulness.
