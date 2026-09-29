# Candidate method: uncertainty separated by source

Development notes, 2026-09-29. Not final claims or results.

## Scientific target

Quantify uncertainty of a specified linear density functional (a local average,
contrast between regions, or Fourier component), conditional on an explicit
class of density and imaging-parameter errors. Distinguish this from predictive
coverage of held-out noisy images, conformational variability, and confidence in
map enhancement relative to a processed reference.

The central failure to test is narrow but wrong intervals under regularization
and pose error, including cases where independent half reconstructions agree.
This failure and basic Bayesian uncertainty are established prior work; the
candidate contribution is a computable, nuisance-aware confidence procedure with
explicit separate error terms and an adversarial calibration benchmark.

## Local experiment and a confidence certificate

Use an independent pilot to select a Fourier Gaussian dictionary, pilot density
c0, noise whitening, functionals l, and tuning constants. The inference images
obey the realified, whitened local model

    y_i - A_i c0 = A_i delta + J_i u_i + r_i + epsilon_i,
    ||L^(1/2) delta|| <= B,  ||u_i|| <= eta_i,  ||r_i|| <= gamma_i,
    epsilon ~ N(0, I).

Here delta=c*-c0; J_i is the pose/CTF/scale derivative of the pilot forward
image. The first-order nuisance term is distinct from the remainder. For any
weights w fixed independently of the inference noise, the estimator and bound are

    t_hat = l'c0 + sum_i w_i'(y_i-A_i c0),
    bias_bound = B ||L^(-1/2)(l-A'w)||
                 + sum_i eta_i ||J_i'w_i|| + sum_i gamma_i ||w_i||,
    half_width = z_(1-alpha/2) ||w|| + bias_bound.

This is a direct support-function construction in the optimal-recovery/honest-CI
tradition (Donoho 1994; Armstrong and Kolesar 2018), not a claim of a new general
inference theorem. Prove coverage with explicit conditioning. Also derive a
simultaneous Gaussian/Bonferroni version for a predeclared family of functionals.
No calibration on noisy images can establish this density guarantee by itself.

Important: B, eta and gamma are assumptions, not automatically valid estimates.
Report sensitivity surfaces, deliberate violations and achievable interval width.
A favorable coverage rate with oracle bounds is only a theorem check; it is not
a practical experiment. Investigate defensible pilot-only selection and physical
units for these bounds before the methodology is finalized.

## Nuisance-aware optimization

Minimize the displayed half-width over w. This is a convex sum-of-norms problem.
Pure nuisance orthogonalization J_i'w_i=0 is a limiting baseline, not automatically
optimal: it can discard substantial structural information. Joint Hessian/pose
sensitivity is already studied by Rangan et al. (2024).

An iteratively reweighted quadratic algorithm uses
||x|| = inf_(t>0) (||x||^2/t + t)/2. At each step its subproblem is

    min_w w'Dw + beta ||l-A'w||^2,
    w = D^(-1) A (A'D^(-1)A + beta^(-1)I)^(-1) l.

D is block diagonal over particles, with an identity term plus a small-rank
J_i J_i' contribution. Apply D^(-1) with Woodbury and solve the coefficient-space
system without forming a dense covariance. Compare small problems against an
independent conic solver, record objective convergence, and never equate an
approximate minimizer with an invalid confidence bound: any fixed w has the same
coverage certificate when its bound is evaluated correctly.

## Fourier Gaussian-specific issues

Analytic derivatives of conjugate-paired Fourier kernels supply J. The old
implementation uses a hard spatial cutoff in Fourier distance. That kernel is
discontinuous at its cutoff and is unsuitable for an unqualified Taylor remainder
theorem. Either use full Gaussians in the theory and charge truncation to r, or
implement a twice continuously differentiable compact taper and derive its
derivative bounds. Do not silently ignore truncation discontinuities.

For a smooth forward map A_i(theta), a local density-error ball and a nuisance
ball imply a remainder of order

    ||r_i|| <= eta_i B sup ||D_theta A_i L^(-1/2)||
               + (eta_i^2/2) sup ||D_theta^2 A_i c0||.

Suprema must cover the entire nuisance ball. A Hessian evaluated only at the
pilot is not a certified supremum. Gauge ambiguity (global rotation/translation,
amplitude) must be fixed or reflected in interval widths.

## Experiments required before a paper claim

- Correct-model repeated-sampling coverage, uniform over difficult directions,
  comparing full Gaussian posterior, diagonal variational covariance, sampling
  covariance, bootstrap, unregularized inference where possible, and the proposed
  certificate; equal estimands and transparent prior/bound choices.
- Pose perturbations (random and coherent/adversarial), CTF error, missing views,
  correlated or heavy-tailed noise, representation mismatch, heterogeneity.
- Width-versus-coverage, sensitivity to assumed bounds, power to certify signs,
  runtime/memory, finite numerical tolerance, multiple seeds and Monte Carlo CIs.
- Real-geometry simulations from independent molecular maps; any truth-informed
  hyperparameters explicitly labeled oracle and excluded from main comparisons.
- Real-particle held-out prediction and reconstruction stability with source-group
  splits, without asserting coverage of an unavailable true density.
- Retain negative results. If bounds are too wide to be useful or closest prior
  work already subsumes the contribution, revise the method rather than relabeling
  the result as publication-ready.
