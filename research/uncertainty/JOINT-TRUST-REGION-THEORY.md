# Joint density/pose design and a residual-controlled trust-region audit

30 September 2026. This development follows the triangle-objective cubic
studies and the experimental noise-calibration failures. No empirical outcome
is asserted here. General robust least squares, trust-region duality, Krylov
spaces and cutting planes are classical, not claimed as new theorems.

Primary sources consulted: El Ghaoui and Lebret (1997), *Robust Solutions to
Least-Squares Problems with Uncertain Data*, SIAM J. Matrix Anal. Appl. 18(4),
1035–1064, DOI 10.1137/S0895479896298130, author's full-text PDF
https://people.eecs.berkeley.edu/~elghaoui/Pubs/rob-ls.pdf (abstract/robust
formulation inspected); Adachi, Iwata, Nakatsukasa and Takeda (2017), *Solving
the Trust-Region Subproblem By a Generalized Eigenvalue Problem*, SIAM J.
Optim., https://doi.org/10.1137/16M1058200 (publisher abstract inspected,
including its handling of the hard case). No full-paper-reading claim is made.

Let the target be ell in L2 of the cube; h(w)=ell-A* w; F(w) be the scaled
cubic pose-field operator, linear in w; and L be the existing enlarged pose
radius. The unchanged density uncertainty ball has radius B. The joint
surrogate uses J(w)=B sup_{||v||<=L} ||h(w)-F(w)v||, followed by the same known-
pilot, quartic-remainder and noise terms. For each fixed v, the expression
inside the supremum is the norm of an affine function of w, so J is convex.
This retains the shared residual/pose geometry and is at most the former
triangle term B(||h(w)||+L||F(w)||). The global lifted ball is still a
relaxation of the physical polynomial-pose set. No pose set is reduced.

## Residual-controlled upper proposition (real arithmetic)

For a fixed selected w, write G=F_Q(w)*F_Q(w) for the quadrature pose Gram,
and b for the implemented continuous-target/quadrature-nominal cross vector.
Assume the existing enclosures give

- ||h|| <= H;
- ||G_cont-G|| <= epsilon_G;
- ||F_cont* h-b|| <= epsilon_b;
- lambda_max(G) <= U on the separately drawn spectral event.

For any lambda>U, any x in the pose-coordinate space, and r=b-(lambda I-G)x,

    sup_{||v||<=L} ||h-F_cont v||^2
    <= H^2 + L^2 epsilon_G + 2L epsilon_b
       + lambda L^2 + 2 b'x - x'(lambda I-G)x
       + ||r||^2/(lambda-U).

Proof: expanding the norm and using the two integration errors reduces the
pose-dependent expression to v'Gv-2b'v. Put D=lambda I-G positive definite.
Completing the square gives v'Gv-2b'v <= lambda L^2+b'D^-1 b. For arbitrary
x, b'D^-1 b=2b'x-x'Dx+r'D^-1 r, and D >= (lambda-U)I. These are identities
and inequalities, not convergence claims about the solver. Add the original
pilot term and (B+P) times the quartic field remainder after multiplying the
square root by B. An H bound is not substituted for the cross-vector error.

All choices of lambda and x are valid on the SAME spectral event for fixed w.
Thus minimizing their computed upper expressions and the old cross/triangle
upper bounds requires no additional probability allowance in real arithmetic.
This does not justify choosing new weights after seeing the audit probes.
Existing quadrature error bounds and normal-noise/class premises still apply;
ordinary floating-point and NUFFT errors are guarded heuristically, not by
validated arithmetic. The implementation states that limitation explicitly.

A cached fully reorthogonalized Krylov basis Q stores explicitly computed GQ.
For each lambda, solve the small Galerkin system and evaluate r from Q and GQ;
no exact Lanczos recurrence or exact orthogonality is assumed in the identity.
A truncated or inaccurate reduced solution only changes the residual term.
The spectral event covers the Gram, not errors in its floating-point actions.

## Design and independent checks

For a fixed pose witness v, ||ell-(A*+F(.)v)Q_w c|| <= t is a second-order
cone constraint in reduced weights c. A cutting-plane master accumulates
these constraints; a small Krylov trust-region problem supplies new feasible
witnesses. A separate random start supplements the cross-vector start to
help see a leading eigenvector orthogonal to b (the classical hard case).
Pointwise target quadrature, PSD factor guards, approximate separation and
conic stopping gaps guide design only. They do not certify a continuous
optimum. The final audit recomputes the continuous density/cross bounds and
uses fresh spectral probes after fixing the chosen weights.

Independent tests compare (1) analytic hard/near-hard solutions, (2) the
resolvent expression with direct linear solves and a deliberately nonorthogonal
cache, (3) a separately assembled semidefinite trust-region dual, and (4) a
robust-norm S-lemma SDP against the joint cutting-plane design. A first near-hard
witness test exposed root-tolerance sensitivity; the feasible maximizing
PSD direction is now extended to the ball boundary, and both attempts are
retained. Passing finite checks does not validate experimental premises.
