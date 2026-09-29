# Working theory and proof audit

These are development statements. General confidence-set and convex-duality
principles are classical; the intended contribution is their physically explicit,
computable cryo-EM specialization and validation, not priority for these principles.

## 1. A finite-sample certificate with separate error sources

Condition on an independent pilot and on the inference design. Suppose

    y = A c0 + A delta + J u + r + epsilon,
    epsilon ~ Normal(0,I),
    ||delta|| <= B, ||u_i|| <= eta_i, ||r_i|| <= gamma_i.

J is block diagonal over particles. Let l be a specified density functional and
w be any weights determined without the inference noise. Define

    t_hat = l'c0 + w'(y-Ac0),
    s = ||w||,
    b = B||l-A'w|| + sum_i eta_i||J_i'w_i|| + sum_i gamma_i||w_i||.

Then an interval centered at t_hat with half-width q(s,b,alpha) has conditional
coverage at least 1-alpha for l'(c0+delta), uniformly over the displayed class.
Here q solves Phi((q-b)/s)-Phi((-q-b)/s)=1-alpha when s>0, and q=b when s=0.

Proof: the estimation error is w'epsilon + (A'w-l)'delta + sum_i w_i'J_i u_i
+ sum_i w_i'r_i. Cauchy-Schwarz bounds the deterministic part by b. The Gaussian
part has variance s^2. For a symmetric interval, normal coverage is an even,
nonincreasing function of the magnitude of its mean; its minimum on [-b,b]
occurs at an endpoint. The critical-value equation gives the result. No claim
about a posterior law is needed. The looser width z_(1-alpha/2)s+b is also valid.

For a finite family selected independently of inference noise, applying this
result with alpha/K and the union bound yields simultaneous coverage. Gaussian
maxima can give a tighter simultaneous noise term if the covariance of the fixed
weights is included. Searching over targets on the inference residuals without
such a correction is invalid.

## 2. Convex optimization and a computable dual gap

Minimize F(w)=z||w||+B||l-A'w||+sum eta_i||J_i'w_i||+sum gamma_i||w_i||.
Its Fenchel dual is

    maximize l'v
    subject to ||v||<=B, ||t_i||<=eta_i, ||r_i||<=gamma_i,
               ||A v - J t - r|| <= z.

Proof: represent each norm as the support function of its dual Euclidean ball.
The infimum over w is finite precisely when the sum of the dual linear terms
vanishes. Eliminating the dual variable for z||w|| gives the last constraint.
The primal is finite and continuous everywhere, and z>0 makes it coercive;
standard finite-dimensional Fenchel duality gives equality of the optima.
Every feasible dual tuple is a lower bound even before convergence.

The identity sqrt(||x||^2+e^2)=min_(t>0) [(||x||^2+e^2)/t+t]/2 gives a quadratic
majorization iteration. Holding t fixed yields

    w = D^(-1) A (A'D^(-1)A + beta^(-1)I)^(-1) l,

where D_i=d_i I+k_i J_iJ_i'. The nuisance rank is small. Woodbury permits exact
block inversion, while conjugate gradients avoid storing the full covariance.
The unsmoothed primal objective and an explicitly feasible scaled dual candidate
are recorded. The gap certifies accuracy for F, not optimality of the sharper
folded-normal interval. The inference result holds for every fixed iterate.

## 3. Nonlinear nuisance remainder

For each particle let A_i(theta) be twice differentiable along the segment
theta0+t*u, t in [0,1]. Suppose ||u||<=eta and ||delta||<=B. Let

    J_i u = D A_i(theta0)[u] c0.

If L_i bounds ||D A_i(theta)[v]||_op and H_i bounds
||D^2 A_i(theta)[v,v] c0|| over that segment for every unit v, then

    A_i(theta0+u)(c0+delta) = A_i(theta0)(c0+delta) + J_i u + r_i,
    ||r_i|| <= eta B L_i + (eta^2/2) H_i.

Proof: decompose the remainder into [A_i(theta0+u)-A_i(theta0)]delta and the
second-order Taylor remainder on c0. Integrate the first and second derivatives
along the segment. A derivative evaluated only at theta0 is insufficient.
Use a Lie-algebra rotation path R exp(t*[u]_x), rather than subtracting Euler
angles; include chart scale factors for joint rotation, translation, and CTF.

## 4. Noise scaling and the cost of nuisance protection

Let P project onto the orthogonal complement of the block nuisance columns.
On a chosen coefficient subspace with its gauge fixed, assume A'PA >= kappa*n*I.
Choose w=PA(A'PA)^(-1)l. Then A'w=l, J'w=0 and

    ||w|| <= ||l||/sqrt(kappa*n).

If every remainder bound is gamma, Cauchy-Schwarz gives
sum_i gamma||w_i|| <= gamma*sqrt(n)||w||. Consequently the sum-width is at most

    ||l||/sqrt(kappa) * [z/sqrt(n) + gamma].

With the previous remainder result, this separates a vanishing noise term from
an O(eta*B+eta^2) local-model term. It does NOT promise shrinking uncertainty for
fixed pose error, nor does it apply to an unidentified gauge direction. Since
the optimized certificate minimizes F over all w, its objective is no larger
than that of this feasible nuisance-orthogonal choice.

### Structured interaction refinement

Let u be the normalized joint pose vector, ||u_i||<=1, and define
D_ia = partial A_i / partial u_a at zero. Keep the bilinear term rather than
charging all of it to an isotropic remainder:

    A_i(u_i)(c0+delta) = A_i(0)(c0+delta)
        + sum_a u_ia D_ia c0 + sum_a u_ia D_ia delta + e_i.

If H_i bounds ||D^2 A_i(tu)[v,v] c0|| and K_i bounds
||D^2 A_i(tu)[v,v]||_op over the entire unit ball and unit directions,
Taylor's integral formula gives ||e_i|| <= (H_i+B K_i)/2.
The radii in physical units are already included in u's parameterization.

For T_i(w_i)=[D_i1' w_i,...,D_iq' w_i],

    |sum_a u_ia w_i' D_ia delta|
       = |delta' T_i(w_i) u_i|
       <= B ||T_i(w_i)||_op <= B ||T_i(w_i)||_F.

Consequently the previous coverage theorem applies after adding
B sum_i ||T_i(w_i)||_F to the bias and replacing the coarse remainder
by (H_i+B K_i)/2. This is a relaxation: delta is shared across all particles,
and delta*u_i' is rank one, neither of which is exploited by independent
Frobenius balls. An optimality statement is only for this relaxed certificate.

Each extra Frobenius norm is a Euclidean norm after vectorization. The same
primal-dual and quadratic-majorization derivations hold for multiple groups;
D_i becomes d_i I + sum_g k_ig J_ig J_ig'. Invert in the smaller of the
combined nuisance-coordinate dimension and image-pixel dimension. Independent
conic checks cover both implementations.

The full-Gaussian implementation bounds radial derivatives on the entire
rotation ball. For a Gaussian of width sigma and distance r, use
g0=exp(-r^2/(2 sigma^2)), g1=r*g0/sigma^2, and the conservative Hessian bound
g2=(1+r^2/sigma^2)*g0/sigma^2. Their maxima on a distance interval occur at
its smallest endpoint for g0 and at the point nearest sigma for g1 and g2.
For frequency k and rotation radius a, the displacement is at most
2||k|| sin(a/2), speed a||k||, and acceleration a^2||k||. Translation of radius
s has phase speed at most 2*pi*||q||*s/box. The product rule then bounds the
second derivative by speed^2*g2 + acceleration*g1 +
2*speed*phase_speed*g1 + phase_speed^2*g0, separately for each conjugate pair.
Summing squared coefficient-column bounds yields a Frobenius upper bound on K_i.
This theorem assumes fixed CTFs; unmodeled CTF error is not covered.

### Two further computable interaction bounds

The shared delta also gives

    |sum_i delta' T_i(w_i) u_i|
        <= B sum_i ||T_i(w_i)||_op,
    |sum_i delta' T_i(w_i) u_i|
        <= B sqrt(n) ||[T_1(w_1),...,T_n(w_n)]||_op.

The second bound uses the same density perturbation across all particles and
||concat(u_i)||<=sqrt(n). Either can be smaller than the independent Frobenius
relaxation. For fixed weights, the minimum of valid deterministic upper bounds
is still valid. However, the minimum of convex bounds need not be convex, so
the existing optimization gap does not certify optimality for that minimum.
These are development refinements to evaluate, not implemented final claims.

### Density-energy sensitivity class

For the Gaussian feature dictionary, the real coefficient Gram matrix G is
block diagonal between real and imaginary coefficients. With

    P_jk=(2*pi*sigma_j^2*sigma_k^2/(sigma_j^2+sigma_k^2))^(3/2)
         * exp(-||mu_j-mu_k||^2/(2*(sigma_j^2+sigma_k^2))),

and Q_jk replacing mu_j-mu_k by mu_j+mu_k, its blocks are 2(P+Q)
and 2(P-Q). This follows by integrating products of Gaussians; Parseval then
gives ||rho_delta||_L2^2=delta'G delta. On a declared non-null dictionary
subspace, choose T with T'G T=I and set delta=T b. Transform A to A T,
the target l to T'l, and derivative matrices to D_ia T. The coefficient bound
||b||<=B now means a whole-space integrated density-energy bound. This does
not supply B from data or cover density components outside the dictionary.
Excluding small numerical eigenmodes narrows the declared model subspace and
must be recorded rather than treated as harmless uncertainty removal.

## 5. A necessary uncertainty floor

If h lies in ker(A), the two parameters c0 +/- B*h/||h|| induce the same data
law under zero nuisance and remainder. Their target values differ by
2B*|l'h|/||h||. Any interval with uniform coverage 1-alpha must contain both
target values with probability at least 1-2alpha under their common law; hence
its expected length is at least (1-2alpha) times that separation. Maximize over
h to obtain the projection onto ker(A). For the proposed deterministic-width
sum certificate, ||l-A'w|| >= ||Proj_ker(A) l|| directly yields its bias floor.
An arbitrarily narrow image-prediction interval cannot remove this obstruction.

## Assumptions still requiring practical work

- Bounds B, eta, gamma and the noise law must be justified, or the result must be
  explicitly reported as conditional sensitivity analysis.
- Learned dictionaries, c0, target selection, noise scaling and regularization
  chosen on the inference data need additional analysis; independence is the
  initial implementation's deliberate design.
- Supplied consensus poses estimated using all particles break a literal claim
  of independent end-to-end reconstruction. Real-data results are conditional.
- Micrograph correlations and windowing violate independent scalar Gaussian
  noise. Whitening, block-aware diagnostics and violation experiments are needed.
- Homogeneous density uncertainty does not validate a conformational population.

## Numerical checks completed so far

2026-09-29: five tests pass: folded-normal coverage equation, scalar analytic
optimum/unidentified case, independent CVXPY/CLARABEL primal optimum and dual gap,
100,000 draws at a simultaneously adversarial bias/nuisance/remainder boundary,
and agreement of sparse versus dense implementations. These are implementation
checks and do not yet establish cryo-EM performance.

## 6. Audit the density space, not just the fitting dictionary

Let H be a declared finite voxel space (possibly restricted by a specified
support mask), A:H->R^m its forward operator, and ell in H a fixed functional.
Suppose ||rho-rho0||_H<=B. For arbitrary fixed image weights w the correct
reconstruction-bias support function is

    b_H(w)=B ||ell-A* w||_H.

If S:R^p->H has orthonormal columns, a dictionary-restricted computation instead
uses b_S(w)=B||S*(ell-A*w)||. It can be arbitrarily smaller: by Pythagoras,

    b_H(w)^2 = b_S(w)^2 + B^2 ||(I-SS*)(ell-A*w)||_H^2.

Thus increasing the forward interpolation accuracy within S cannot establish
confidence over H. Evaluating the full adjoint residual does account for errors
outside S, without separately estimating their direction from image residuals.
It still requires the total radius B and the ambient grid/support assumptions.

**Coverage.** Use the same center ell'rho0+w'(y-A rho0) and replace the restricted
bias by b_H. The theorem in Section 1 applies directly; weights may have been
computed using any Gaussian or neural surrogate. Centering must use the ambient
operator too, or its deterministic approximation error must also be bounded.

**Reduced-space lower bound.** The optimum of
F_S(w)=z||w||+B||S*(ell-A*w)|| is at most the full optimum F_H. More usefully,
a feasible restricted dual v has ||v||<=B and ||A S v||<=z; its lift S v is
feasible for the full dual and has objective ell'S v. Therefore a full adjoint
residual supplies a primal upper bound and a restricted dual supplies a lower
bound on the *same ambient* optimization problem. Their gap is a legitimate
stopping criterion. This is a specialization of reduced-space convex duality,
not a claim of a new general duality theorem.

**Enrichment.** Append the normalized component of ell-A*w orthogonal to S.
The subspaces are nested, so their exact restricted optima are nondecreasing and
bounded by the ambient optimum. Store the best ambient primal iterate rather
than assuming every enrichment decreases its objective. If the outside residual
is zero and the restricted primal-dual gap is zero, the current w solves the
ambient problem. A finite gap certifies the corresponding tolerance. No rate
for this greedy scheme is asserted. All choices depend on the design and
functional, not inference noise, so enrichment preserves the conditional result.

**Matrix-free comparator.** With fixed poses, the majorization subproblem is
(A*A+lambda I)u=ell and w=A u. For any approximate u, set
v=u min(B/||u||,z/||A u||). This is a feasible ambient dual point independently
of CG convergence. The implemented full-grid solver is a comparator to reduced
Gaussian enrichment, not evidence that the latter is invariably faster.

FINUFFT forward/adjoint evaluations are checked against direct Fourier sums
and adjoint identities at tight tolerances. These are floating-point numerical
certificates, not interval-arithmetic bounds on roundoff. A rigorous operator
error bound epsilon_A, if available, can be charged by adding B epsilon_A to
the evaluated adjoint residual norm and separately bounding pilot-centering
error. Tests do not by themselves prove a uniform NUFFT error bound.

## 7. Near-indistinguishability gives an uncertainty floor

Let two allowed densities/nuisances have Gaussian means mu0,mu1 and target
values t0,t1, with identity noise covariance. Their total variation distance is
TV=2 Phi(||mu1-mu0||/2)-1. If an interval has coverage at least 1-alpha at both
parameters, then under the first law it contains both targets with probability
at least max(0,1-2alpha-TV): transfer the second coverage event to the first law
and apply the union bound. Consequently its expected length under the first law
is at least |t1-t0| max(0,1-2alpha-TV). Exact null directions recover Section 5.
This familiar two-point argument shows why confidence cannot be rescued merely
by accurate image prediction or small half-map disagreement. Searching a finite
set of pairs gives a lower bound, not the exact modulus over all densities.

## Solver provenance

The additional implementation uses classical Chambolle--Pock primal-dual hybrid
gradient, with blockwise dual step sizes. For K_g of norm at most L_g, choosing
tau=.99/sum L_g and sigma_g=.99/L_g gives
||sqrt(Sigma) K sqrt(tau)||^2 <= tau sum sigma_g L_g^2=.99^2<1.
Zero operators receive arbitrary finite positive steps. The noise and per-image
remainder norms have nested groups, so their proximal map first shrinks each
image block and then the entire vector. Independent conic solves validate both
this implementation and the original majorization solver, including rescaling.
See Chambolle and Pock (2011), doi:10.1007/s10851-010-0251-1.

## 8. Combining ambient auditing with nonlinear pose bounds

For a finite supported voxel model, its complex Fourier column at voxel x is
C_i(q) exp(-2 pi i k(u)'x / box) exp(-2 pi i q't(u)/box), divided by the known
real/imaginary noise standard deviation. This gives full-space derivatives,
without assuming a Gaussian dictionary spans the true density.

Let a,s be the physical rotation/shift radii of the normalized joint unit ball.
Along any segment tu, frequency rotation preserves ||k||. Define per voxel and
frequency

    L(q,x)=2 pi/box * (a ||k|| ||x|| + s ||q||),
    H(q,x)=2 pi/box * a^2 ||k|| ||x||.

The phase derivative magnitude is at most L and its second derivative at most
H throughout the segment. Since |exp(i phase)|=1, the complex column second
 derivative is bounded by |C|/noise * (L^2+H). Summing these positive column
bounds against |rho0_x| gives the pilot bound H_i after taking the frequency
Euclidean norm. Their frequency/voxel Frobenius norm gives K_i. Taylor's formula
therefore gives gamma_i=(H_i+B K_i)/2 for any ||delta||<=B in the support.
The realification norm equals the complex Euclidean norm, so no duplicate
Friedel-pair or real/imaginary factor is added.

The exact first derivative of each column is i times its phase derivative times
the column. Its pilot contraction supplies J_i. Flattening its density/pose
indices gives D_i; a factor of D_i D_i' supplies the support function
B||D_i' w_i|| without storing the full tensor for every particle. A small
positive Cholesky padding enlarges this bound. The resulting joint confidence
certificate uses the full ambient adjoint residual, pilot-pose term, full-density
interaction term and uniform remainder in a single objective.

The general majorization solver now accepts a forward/adjoint LinearOperator.
For the Fourier voxel operator with independent real/imaginary samples,
diag(A'A) is exactly sum_{i,q} |C_i(q)/noise|^2 for every supported voxel:
cos^2+sin^2=1. This supplies its CG preconditioner without forming the design.
Particle-block nuisance inverses are small. The noise, support and local-pose
assumptions still require validation; integration does not make them automatic.

## 9. Independent calibration of a scalar noise upper bound

This is a standard chi-square confidence construction, not a new distributional
result. Suppose calibration coordinates z_j=mu_j+sigma e_j, j=1,...,nu, have
independent standard Gaussian e_j and arbitrary fixed means. Their noise is
independent of the inference noise. Set

    sigma_upper = ||z|| / sqrt(chi2_quantile(beta, nu)).

Then P(sigma_upper >= sigma) >= 1-beta. At mu=0 this is the exact chi-square
pivot. For nonzero mu, ||z/sigma||^2 is noncentral chi-square; its lower-tail
probability at a fixed threshold is at most that of the central distribution.
Thus residual signal increases conservatism; it need not be declared absent.
Subtracting a fitted mean without adjusting the argument is not this procedure.

Conditional on the independent calibration, weights may be chosen using its
upper scale. On the scale-coverage event, the inference noise SD of a linear
estimator is at most s_upper=sigma_upper ||w||. For alpha<1/2 the folded-normal
critical width q(s_upper,b,alpha) is at least b. Its coverage for any |bias|<=b
is nonincreasing in the actual SD over [0,s_upper], so the conditional coverage
is at least 1-alpha. Integrating gives at least (1-alpha)(1-beta), and therefore
at least 1-alpha-beta. Multiple upper scales use an explicitly allocated beta
budget and the corresponding diagonal weighted norm. This does not estimate
arbitrary correlations, justify solvent homogeneity, or make windowed Fourier
coordinates independent. Those assumptions still need separate diagnostics.

## 10. Retaining pose curvature instead of an isotropic quadratic remainder

The first-order certificate may pay heavily for a direction-independent bound
on the entire second-order change. A higher-order specialization keeps that
curvature inside weight-dependent support functions. This uses Taylor's theorem
and Euclidean relaxations; no novelty is claimed for Taylor bounds themselves.

Let E_iab be the second pose derivative of A_i at zero, and define
Q_i(w)_ab=w_i' E_iab rho0 and U_i(w) with columns E_iab' w_i for all ordered
pairs a,b (symmetric off-diagonal terms appear twice). For ||u_i||<=1,
||u_i tensor u_i||_F=||u_i||^2<=1, so the quadratic pilot and density terms are
bounded by ||Q_i(w)||_F/2 and B||U_i(w)||_F/2. Add these to the previous density,
pilot-linear and density-linear-interaction terms. Replace the isotropic
second-order remainder by a uniform cubic remainder.

For the finite Fourier voxel column, retain L,H from Section 8 and set
T(q,x)=2*pi/box * a^3 ||k|| ||x||. Throughout each allowed pose segment,

    |column'''| <= |C|/noise * (L^3 + 3 L H + T).

Indeed differentiating exp(i phi) three times gives terms with magnitudes
|phi'|^3, 3|phi'||phi''| and |phi'''|; rotation preserves frequency length,
and translation is linear. Contract the positive column bound with |rho0| for
a pilot norm bound H3_i, and use its Frobenius norm for K3_i. Taylor's integral
remainder is at most gamma3_i=(H3_i+B K3_i)/6. Thus the general coverage theorem
applies with all four structured nuisance terms and the cubic remainder.
It remains a relaxation: the repeated density perturbation and the rank-one
quadratic pose tensor are not jointly optimized over their exact set.

For implementation, phi_a uses x cross k for rotations and q for translations.
Its only nonzero second derivatives are the rotation block:

    phi_ab = (-2*pi/box)*a^2 *
             [ (k_a x_b+k_b x_a)/2 - 1(a=b) k'x ].

The column Hessian is exp(i phi)*(i phi_ab - phi_a phi_b), multiplied by C/noise.
Pixel-space Gram factors compress the density/Hessian interaction exactly up
to the same conservative positive numerical padding as the first derivative.
Both pilot and interaction norms are Euclidean groups, so the existing convex
solver and feasible dual diagnostics apply without a new optimization claim.
Finite differences check diagonal/mixed rotation and translation entries;
independent nonlinear projections check the uniform cubic remainder. Whether
this tighter expansion actually improves useful intervals is an experiment,
not an assumption from its higher order.

## 11. Auditing pose terms without building full-grid derivative tensors

For fixed image weights, the first and second pose-derivative adjoints can be
computed from twenty Fourier moment fields. Let f=(k_x,k_y,k_z,q_x,q_y),
c_q=(w_real+i*w_imag) C_q/noise, and define
F_r(x)=sum_q c_q f_qr exp(+2*pi*i*k_q'x/box),
F_rs(x)=sum_q c_q f_qr f_qs exp(+2*pi*i*k_q'x/box).
There are five first moments and fifteen distinct symmetric second moments.

The positive-phase derivative is linear in f: psi_a(x,q)=M_ar(x) f_qr.
M contains the cross-product coefficients of x with k for rotation, and scaled
q components for translation. The rotation phase Hessian similarly has the
linear form psi_ab=N_abr(x) f_qr. Then

    D_a' w(x) = Re[ i sum_r M_ar F_r ],
    E_ab' w(x) = Re[ i sum_r N_abr F_r - sum_rs M_ar M_bs F_rs ].

Pilot inner products and supported-voxel Frobenius norms give exactly the same
structured bounds as the full derivative tensors. This avoids forming their
large pixel-space Grams when auditing fixed coarse-space weights on a finer
space. The fields can be evaluated by batched type-1 NUFFTs or direct Fourier
sums; the direct path is currently selected for at most 128 frequencies. That
threshold is an implementation heuristic, not an optimal complexity theorem.
Both paths are tested against explicit column derivatives, including mixed
rotation/translation curvature. For a new ambient grid, centering and the
target also use that grid, and physical voxel-volume scaling must be preserved.
These are larger finite-space guarantees, not assertions about all continuous
unresolved density.

## 12. Shared-density spectral post-audit and feasible nonlinear adversaries

The earlier sum B*sum_i ||T_i||_F allows a different worst density direction
per particle. Retain the shared density vector in the linear/quadratic Taylor
polynomial instead. Let h=A'w-ell, and form blocks H_0=h, H_i=T_i, and
H_(n+i)=U_i/2. Their coefficient vectors are 1, u_i, and u_i tensor u_i,
respectively. Every block coefficient has norm at most one. For arbitrary
positive deterministic block scalings d_j, Cauchy--Schwarz and the operator norm
give the uniform inequality

    ||sum_j H_j v_j|| <= sqrt(sum_j d_j) ||[H_j/sqrt(d_j)]||_op.

Multiplying by B bounds the combined density residual/interaction, retaining
its shared direction. This is an elementary classical spectral inequality.
The pilot-linear, pilot-quadratic and uniform cubic remainder bounds are added
unchanged. No coverage assumption is removed. The coefficient-one and quadratic
rank-one constraints are relaxed, so the bound can still be conservative.

The implementation compares the block triangle bound, separate linear and
quadratic spectral bounds, an affine residual-plus-linear spectral bound, and
one joint bound. Scales d_j=||H_j||_F are deterministic for fixed weights; zero
blocks are dropped. The minimum remains a uniform upper bound. All candidates
are independent of inference noise. Dense symmetric float64 eigensolves have
a small roundoff pad; no interval-arithmetic claim is made. The original convex
solver's objective gap does not imply optimality of this post-audited minimum.

A complementary nonlinear stress test eliminates the density ball exactly.
For each feasible pose u, put h(u)=A(u)'w-ell. The largest signed bias over that
ball is

    F_sign(u)=sign * rho0'[A(u)'w-A(0)'w] + B||h(u)||.

It is attained by delta=sign*B*h(u)/||h(u)|| when h(u) is nonzero. Projected
Adam on each pose unit ball therefore finds feasible lower bounds on worst
bias. Four starts for each sign are used in current development tests, with
150 iterations per start. This is a local search, not a global optimum or proof
of worst-case coverage. MPS float32 accelerates the search; every saved pose
is reprojected onto its ball and recomputed by independent SciPy rotations and
direct complex Fourier sums in float64. Autograd checks and explicit first
pose-derivative comparisons validate the search objective and gradient.
