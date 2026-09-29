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
