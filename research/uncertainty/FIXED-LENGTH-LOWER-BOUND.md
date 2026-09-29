# A classical two-point lower bound for fixed-length intervals

This specializes standard Gaussian testing/optimal-recovery ideas, not a new
general lower-bound theorem. It supplies a quantitative comparison for the
continuous fixed-pose experiments without fitting to simulation outcomes.

Let y ~ N(A rho, I), ||rho-rho0||_H <= B, and t=<ell,rho>. Suppose a procedure
returns intervals of deterministic length 2q and has marginal coverage >=1-alpha
for every density in this class, with 0<alpha<1/2. Its center may be nonlinear
and data dependent. For any h in H and any lambda with

    lambda ||h|| <= B, lambda ||A h|| < z_(1-alpha),

consider rho+/-=rho0 +/- lambda h. Their target separation is
2 lambda |<ell,h>|. Their observation laws have total variation
2 Phi(lambda ||A h||)-1 < 1-2alpha. If the deterministic interval length is
smaller than their target separation, the two containment events are disjoint.
Coverage >=1-alpha under each law would then require total variation
>=1-2alpha, a contradiction. Taking limits at the strict boundary yields

    q >= |<ell,h>| min(B/||h||, z_(1-alpha)/||A h||).

When Ah=0, the second constraint is inactive. This lower bound applies to every
fixed-length interval, not just affine centers. It does not bound the length
of a randomly varying interval at every observation. For expected length of
arbitrary intervals, the separate weaker TV expression in the paper applies.

The existing continuous solver supplies h=ell-A* w. Its norm and target pairing
are analytic; Ah=A ell-Gw. Existing Gauss remainders bound integration errors.
Upper bounds on both norms in the denominators preserve a conservative lower
bound in real arithmetic. Ordinary roundoff is not formally enclosed. Also
evaluate h=ell, retaining the larger valid lower bound. This post-hoc analysis
does not alter a frozen estimator or study outcome.

Fixed poses form a submodel of the bounded-pose class, so the lower bound remains
necessary there. A large gap to a pose-robust upper bound does not prove that
the gap is unavoidable under pose uncertainty. The comparison separates
fixed-pose efficiency from additional pose conservatism.
