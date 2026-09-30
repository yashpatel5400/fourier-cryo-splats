# Spectral support exchange for pose-aware weight fitting

This is exploratory numerical-method development, not a new statistical
coverage theorem. It addresses poor convergence of the smoothed Ritz/L-BFGS
implementation. The final continuous density class and pose model are unchanged.

Write the previously defined objective as

  J(w) = z||w|| + B||ell-A*w|| + p(w) + r(w),

where p(w)=R sqrt(sum d) sqrt(||F(w)||op²+E||M|w|pair||²) is the shared-field
pose upper bound and r(w) is the nonnegative weighted sum of complex-pair
amplitudes for the cubic remainder. Block scales d are fixed from the original
fixed-pose fit. R=B+||pilot|| and F is linear in estimator weights.

## Restricted conic problem and spectral cuts

Restrict weights to w=V beta without restricting the density class. Every pair
of unit vectors (u,v) gives a full-weight linear support

  g^T w = v^T F(w) u <= ||F(w)||op.

For a finite set of such supports, replace the spectral norm by t>=0 and
constraints t>=g_j^T V beta. The noise norm, continuous density residual and
cubic group norms are second-order-cone representable. Omitting the nonnegative
integration pad from this inner problem is an explicit optimization relaxation;
it remains included in the independent final upper audit. A restricted conic
objective is neither a full-space lower bound nor a coverage certificate.

After each solve, the leading matrix-free Ritz directions supply new supports.
Only their unit norms and the adjoint identity are needed for support validity;
a Ritz value is never used as an upper bound. The prototype starts with shell
components of the original and saved weights. Every third iteration adds a
preconditioned full-weight stationarity direction, orthogonalized twice against
existing columns. This enrichment is a heuristic; no convergence rate is claimed.

## Full-space dual support and final audit

If mu_j are conic multipliers for t>=g_j^T w, the t>=0 constraint implies
mu_j>=0 and sum mu_j<=R sqrt(sum d). Numerically project multipliers onto these
constraints before exporting gp=sum mu_j g_j. For EVERY full-space weight x,

  gp^T x <= R sqrt(sum d) ||F(x)||op <= p(x).

Thus gp is a valid full-space pose dual support even when the weight subproblem
is restricted. The cubic norm's supporting gradient gr is similarly feasible.
Choose v=(B/||h||)h for h=ell-A*w, using a padded residual norm in computation.
The remaining noise dual constraint is ||A v-gp-gr||<=z. Scaling v,gp,gr by
min(1,z/defect) makes the full dual feasible. The known continuous target pairing
then yields a full-space lower bound on J. The implementation includes the
quadrature Gram-action error in this defect, as in the earlier solver.

The candidate is selected using optimization diagnostics only. Fresh Gaussian
power probes then supply the final probabilistic spectral upper bound, with
the original integration pad, nonlinear remainder and alpha allocation. The
full objective gap compares this upper to the feasible full-space dual value.
It is not inferred from a conic solver status or a restricted model value.

## Independent check and current probe

An independent small full conic problem uses the exact spectral matrix norm.
The cut implementation brackets its optimum and converges to it on a nonzero
solution; exported dual supports are tested on independent arbitrary weights.
A first random fixture happened to have a zero optimum and was strengthened
before accepting the test. Both the original zero-fixture check and a NumPy
scalar-shape fixture failure are retained in logs; neither is presented as a
completed scientific experiment.

The initial bounded archive-geometry probe uses the completed 10049 center
case at two degrees and 0.5 Angstrom, whose previous sum-objective gap is 99.27%.
It runs 30 conic/support rounds and an independent final audit. Success or
improvement is not assumed in advance; checkpoints and failures are preserved.

The follow-up implementation exports the cubic cone's dual vectors directly.
At nearly zero complex weight pairs, normalizing the tiny primal pair gives
an unstable boundary gradient and can make the full-space dual residual
unnecessarily large. The cone dual can instead choose a vector inside the
known per-pair dual ball. Each exported vector is projected onto that ball,
so it remains a valid global support even with solver tolerances. An additional
independent test has an exactly zero active pair with a strictly interior
optimal subgradient and checks full stationarity. The best feasible full-space
lower bound is retained across iterations, independently of which primal
candidate is best. The original normalized-gradient probe remains unchanged.

## Full observation weights with an inner density projection

A subsequent bounded trial removes the estimator-weight restriction entirely.
Let phi_j be an orthonormal continuous Hilbert basis spanning the target ell
(first) and selected observation representers A_i*. Exact Gaussian/sinc kernel
identities implement pivoted Hilbert Gram-Schmidt. Let F_ij=<A_i*,phi_j> and
b_j=<ell,phi_j>. Then F b=A ell, ||b||=||ell||, and

  ||b-F^T w|| <= ||ell-A*w||     for EVERY full observation-weight vector w.

The conic inner solve uses this smaller density residual, all observation
weights, cubic pair cones and spectral cuts. It is an explicit lower relaxation
of the original objective. No new density-model assumption is made. Final
candidate selection and the final interval still use the full continuous
residual and independent pose audit; poor projected candidates are retained
as diagnostics but cannot silently replace a better ambient candidate.

The conic density dual u lifts to v=sum_j u_j phi_j. Its norm is ||u|| and
A v=F u. After projecting ||u||<=B, the spectral and cubic dual supports as
before, a common scale enforces ||F u-gp-gr||<=z. The target pairing b^T u
is therefore a feasible FULL-space lower bound, not merely a restricted-weight
optimum. In exact arithmetic the projection Gram residual is PSD; tests compare
it with the independent analytic sinc Gram and verify full conic bracketing.
Numerical orthogonality remains subject to the disclosed floating-point caveat.

An independent test caught a wrong sign when exporting the density cone dual;
it was corrected before any archive-geometry trial. The failing test log is
retained. The initial full-weight trial uses rank 256 and six exchange rounds,
with a fresh final spectral seed. Neither completion nor improvement is assumed.

The six-round trial completed in 201.50 seconds. Its best primal candidate
remained the input weights: projected candidates had full Ritz objectives
31.1--32.8, much larger than the input's roughly 11.2. The fresh audit gives
relative width 0.69024 and a 0.76228 gap; the small difference from the input
width is the new random upper certificate, not a changed estimator.

An adaptive follow-up retains all observation weights and adds the omitted
representer after every inner solve. If F,b denote the current continuous
projection and w a candidate, put c=G w-F(F^T w), s²=w^T c. For s²>0,
append c/s as a column of F and zero to b. This is Hilbert Gram--Schmidt on
(I-Pi)A*w. It makes the projected residual exact for this candidate while
preserving a lower relaxation for every other candidate. The bounded-size
probe uses the independently assembled analytic sinc Gram for these products,
not a silently unpadded approximate Gram. Tests check PSD of G-FF^T and exact
representation of all previously added candidates. Significant negative
residuals abort and are preserved; floating-point orthogonality is still a
numerical caveat. The initial follow-up has 30 rounds and a fresh audit seed.
