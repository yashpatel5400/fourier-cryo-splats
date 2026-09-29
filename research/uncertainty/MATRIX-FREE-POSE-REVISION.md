# Matrix-free shared pose fields and numerical probability budgets

Development after the first independent review. This is not a revision of the
frozen v1/v2 outcomes, nor a completed scale or usefulness result.

## Field operator

The first and half-scaled symmetric second pose derivatives have ten spatial
monomials: 1, x, y, z, x², y², z², xy, xz, yz. Their Fourier coefficients depend
on the particle, frequency and derivative direction. Summing coefficient
combinations before transforming evaluates the full field with ten simultaneous
type-3 NUFFTs, rather than materializing twenty fields for every particle.
The transpose transforms ten monomial-weighted quadrature fields and contracts
the resulting frequency coefficients. It is the same quadrature linear operator
as the original dense field matrix, not a lower-resolution approximation.
Storage is linear in the number of particle frequencies and quadrature nodes.

The implementation accepts particle-specific rotation and translation radii.
The five-dimensional unit ball remains joint, not five independent boxes.
The coordinate frame is externally declared: coherent global rotations are
included. Removing those modes would require a gauge-invariant estimand or an
explicit quotient class and is not silently done to obtain smaller widths.

## Spectral upper bounds need more than a Ritz value

For a fixed positive semidefinite matrix C=F F*, independent g_j~N(0,I), and
m probes, put c=sqrt(2) erf^{-1}(delta^{1/m}). Fix a unit eigenvector v for the
largest eigenvalue before drawing the probes. Then

    P(max_j |<v,g_j>| < c) = delta.

Outside this event, for every positive integer p simultaneously,

    lambda_max(C) <= [max_j ||C^p g_j|| / c]^(1/p).

Proof: at least one projection has magnitude at least c, and its component
in C^p g_j has magnitude lambda_max(C)^p |<v,g_j>|. Taking the norm and p-th
root yields the bound. Stable normalization tracks the product of norms in
log space. Taking the minimum across powers is valid on the same event. Any
deterministic trace upper bound may also cap the result. A Rayleigh or Lanczos
Ritz value is only a lower bound and is never used as the uncertainty upper.

This is classical Gaussian randomized norm estimation, not a new general
theorem. Relevant prior art includes Kuczynski and Wozniakowski, *Estimating the
Largest Eigenvalue by the Power and Lanczos Algorithms with a Random Start*,
SIAM J. Matrix Anal. Appl. 13(4), 1094–1122 (1992), DOI 10.1137/0613066.
The publisher record and the introduction of the authors' March 1989 Columbia
technical report CUCS-465-89 have been checked. The latter is a scanned primary
source; its first five pages were rendered and OCR-read. It studies random-start
Krylov error/failure bounds and explicitly discusses the difficulty of stopping
criteria. This is related prior art, not an assertion that our displayed
Gaussian-projection formula is copied from a numbered theorem in that report.
The short argument above is independently stated in full.

For the grouped field inequality with positive scales d, its continuous norm
upper becomes

    sqrt(sum(d) * [lambda_upper(C_quadrature) + integration_pad]).

The existing polynomial quadrature remainder bounds the omitted integral. It
does not bound floating-point/NUFFT errors; those remain numerical checks and
must not be relabeled validated arithmetic. Streaming particle-block norms
constructs d without a global dense Gram. Zero fields contribute zero.

For the 10,000-particle run, a cheaper alternative chooses each d as the
Euclidean norm of its polynomial coefficient envelopes. Every spatial monomial
has absolute value at most one on the unit cube. Consequently the squared
Frobenius norm of the quadrature matrix is bounded by the sum of squared scaled
coefficient envelopes. The block scales need only be positive, not exact field
norms, for the shared-field inequality. These scales affect tightness, not the
matvec's resolution. A block may be omitted only if it remains identically zero;
an initially zero but subsequently optimized block must receive a positive scale.

## Convex estimator objective and a feasible dual lower bound

Fix all positive block scales before optimizing w. Let A be the whitened
nominal continuous-density observation operator, ell the target, B the density
radius about the independent pilot, and R=B+||pilot||. Write F(w) for the scaled
quadrature pose-field matrix; it is linear in real estimator weights. Let
a_j(w)=sqrt(w_Re,j²+w_Im,j²), M(w)=C a(w) for a nonnegative coefficient-envelope
matrix C, E the polynomial integration error constant, and r>=0 the integrated
cubic coefficients. The sum-objective upper bound is

    J(w) = z ||w|| + B ||ell-A* w|| + p(w) + r' a(w),
    p(w) = R sqrt(sum d) sqrt(||F(w)||op² + E ||M(w)||²).

Each a_j is a norm. Because C is nonnegative, ||C a(w)|| is convex: the
componentwise triangle inequality followed by monotonicity of the Euclidean
norm on the nonnegative orthant proves this directly. Composition with the
nonnegative, coordinatewise increasing two-dimensional Euclidean norm proves
convexity of p. All four summands are nonnegative convex functions; the last two
are also positively homogeneous. No differentiability or unique leading
singular vector is assumed. J is an upper objective, not the exact minimax
nonlinear-pose risk or the exact folded-normal interval length.

The optimizer uses approximate leading singular modes. A log-sum-exp over k
singular values smooths the maximum in exploratory optimization. Exact leading
values would give an additive singular-value excess at most mu log(k), but Ritz
values do not supply a spectral upper certificate. No guarantee of final
optimization error is derived from a small gradient or a solver-success flag.
The final reported primal bound uses fresh independent power probes and the
analytic density/pose integration pads.

The following feasible-dual construction is useful even when the optimization
has not converged. Select v in the continuous density Hilbert space with
||v||<=B, and linear functionals g_p, g_r satisfying, for every x,

    g_p' x <= p(x),   g_r' x <= r' a(x),
    ||A v - g_p - g_r|| <= z.

Then J(x)>=<v,ell> for every x: apply the support inequality to each norm and
cancel the linear terms. This is weak convex duality, not a new generic theorem.

Our implementation uses v=B(ell-A* w)/r_plus, where r_plus bounds the exact
continuous residual norm from above. For any unit spatial modes v_j, set
u_j=F(w)'v_j/||F(w)'v_j|| when nonzero. Any convex mixture of v_j u_j' has
nuclear norm at most one; it therefore yields a supporting linear functional
for the spectral norm, even if the modes are not leading eigenvectors. For the
coefficient term, take t=M(w)/||M(w)||>=0. Combining t with the unit complex-pair
directions gives a supporting functional for ||M(x)||. Weight the spectral and
coefficient supporting functionals by nonnegative beta_1,beta_2 with
beta_1²+beta_2²<=1, and multiply by R sqrt(sum d), to obtain g_p. The implemented
beta values use the smoothed spectral value and sqrt(E)||M(w)||; the support
argument does not require that spectral value to be an upper bound.

The cubic g_r is the weighted complex-pair norm subgradient. A padded upper
bound D on ||A v-g_p-g_r|| includes the continuous Gram-action quadrature error.
Scale v, g_p and g_r together by tau=min(1,z/D). The scaled supports remain
feasible because their original support sets contain zero. The reported lower
bound is max(0,tau <v,ell>). An independent small CLARABEL conic optimization
checks this bracketing against a separately assembled exact density Gram and
row-compressed pose matrices. This check does not prove real-data usefulness or
replace the support argument.

## Interaction with inference and selection

Weights, scale choices and optimization iterations may use design information
and separate random draws. The final numerical certificate requires F fixed
BEFORE its fresh Gaussian probes. Reusing probes to choose F invalidates the
simple fixed-matrix probability statement. A collection of candidate matrices
may instead allocate delta by a union bound before selection.

Allocate alpha_noise + delta <= alpha_total. Conditional on design, on the
numerical event the bias-aware interval has Gaussian-noise coverage at least
1-alpha_noise; its overall failure probability is at most alpha_noise+delta.
All numerical randomness must be independent of inference noise. This is an
explicit change from a deterministic spectral calculation, not a claim that
a randomized approximate eigenvalue has the old deterministic guarantee.

The implementation checks dense field agreement, adjoint agreement,
particle-specific radii, estimator-weight gradients, direct matrix powers and
the sufficient Gaussian-projection event on small matrices. Completed 128- and
1,024-particle audits take 50.3 and 353.3 seconds under concurrent Mac load, with
peak resident memory 390 and 434 MB. The latter's relative half-width is 0.740:
scalability does not imply a useful interval. The independent tiny conic check
brackets the optimum with a 0.026% relative primal/dual gap. The first larger
20-iteration pose-aware solve has a 97.8% gap and is not converged. Its width
improvement is exploratory. Higher-band 10,000-particle timing and the longer
optimization remain in progress. No reviewer objection is resolved by this
derivation alone.
