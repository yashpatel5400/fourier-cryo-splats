# Viewing-law misspecification and mixture optimization

30 September 2026. Three primary sources, with targeted reading scopes below;
none is an author-code replication. Local PDFs are hashed but not redistributed.

**Xu, Zhang and Singer (2025), nonuniform orbit MLE.**
[Primary preprint](https://arxiv.org/abs/2509.22945v1).
Read model, sections 3.2--3.3, simulations and discussion; not all proofs.
The projected-model bias theorem requires an analytic parametrization and a
nonzero gradient witness; it is an existence statement about a misspecifying
viewing law, not bias for every nonuniform law. The cryo-EM witness is described
as numerically verifiable, with details omitted. Low-noise bounds require the
true law's density relative to Haar to be bounded above and away from zero,
with further geometric conditions; the initial bound is one-sided. Joint MLE
consistency assumes an identifiable model, bounded structure coefficients and
an in-plane-invariant viewing density of fixed Wigner bandlimit. It is not a
finite-sample confidence theorem or proof of an optimizer's global convergence.
Experiments use continuous MRA, not raw cryo-EM reconstruction. The general
model fixes a known projection operator, with image-specific CTFs noted as
an additional practical complication. Implication: unknown viewing laws are
scientifically consequential, but a blanket claim that uniform priors always
fail would overstate this evidence.

**Lindsay (1983), geometry of mixture likelihoods.**
[Published article](https://doi.org/10.1214/aos/1176346059);
[primary-paper mirror](https://people.csail.mit.edu/jrennie/trg/papers/lindsay-mixture-83.pdf).
Visually read printed pages 87--90, including Theorems 3.1, 4.1 and 5.1 and
Corollary 5.2. For a compact likelihood curve, the maximizing fitted likelihood
vector is unique and admits a mixing measure supported on at most the number
of distinct observations. The mixing measure itself need not be unique.
Independent nonidentical kernels have the analogous n-dimensional statement.
The directional derivative and supporting hyperplane yield a global
optimality condition and dual. If the largest directional derivative is
delta, the likelihood gap is at most n log(1+delta/n), improving the unscaled
tangent bound. These are direct precedents for our finite and continuous
mixture certificates. A finite-support existence theorem does not justify
restricting unknown orientations to a predetermined finite grid.

**Zhang, Cui, Sen and Toh (2024), scalable NPMLE.**
[JMLR article and code link](https://jmlr.org/papers/v25/22-1120.html).
Read fixed-support formulation (6), section 2's primal/dual derivation and
Algorithm 1, plus the algorithmic overview; not all convergence proofs.
The work solves the same finite simplex likelihood problem used in our
prototype. Its augmented Lagrangian method uses semismooth Newton subproblems
and exploits sparsity in the Hessian structure. It compares with established
EM, interior-point and SQP approaches. Neither EM weight fitting nor the
finite convex formulation is novel here. Its main examples are astronomical
empirical-Bayes denoising, not cryo-EM orientation envelopes. Its reported
large-scale performance motivates a solver comparison if our larger finite
mixtures become the bottleneck; we have not measured that author solver on
our likelihood matrices. Continuous-support approximation is distinct from
a finite solver's convergence certificate.

**Development implication.** The prospective contribution must live in useful
continuous imaging-model bounds and validated scientific inference, rather
than in renaming established mixture duality. The current finite-view results
are feasibility screens. The new continuous Gaussian-orbit derivation is in
`CONTINUOUS-MIXTURE-THEORY.md`; computational usefulness is still unestablished.
