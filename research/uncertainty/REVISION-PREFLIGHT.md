# Preflight while the first independent review runs

This note records an implementation constraint, not a reviewer finding.
The scientific packet for round 1 is frozen at `a4117ab`; no scientific
result or manuscript in that packet has changed while review runs.

## Closest published probabilistic reconstruction

Ullrich et al. (PMLR 115, arXiv:1906.07582v2) equation 14 includes both
expected negative image log likelihood and KL divergence to a standard
Gaussian prior. Equation 15 scales the prior KL by minibatch size / total
particles. The main experiments condition on poses, use diagonal Gaussian
Fourier coefficients with Hermitian restrictions, and do not model CTFs.
Their original paper already discusses model bias; a representation-bias
experiment must not be presented as refuting an unclaimed coverage theorem.

The released repository at 013149a0927183b5969095cb774d43b387644771 is a
tutorial/operator release. Its displayed training loop optimizes only sampled
negative likelihood, with log scale initially -100; it is not the full stated
variational objective. Dependencies target Python 3.6 / PyTorch 1.1. Running
this loop unchanged would not produce a faithful full uncertainty baseline.
The authors identify this as an observation-model tutorial, not a maintained
complete benchmarking pipeline.

A credible expanded comparison would explicitly implement equation 14 with
Hermitian parameterization, likelihood/CTF/noise unit checks, and the same
linear density targets. For fixed poses and a linear interpolation operator,
the exact diagonal-Gaussian variational optimum has mean equal to the ridge
posterior mean and coordinate variances equal to reciprocal precision
diagonals. An analytic solution provides an independent ELBO/stationarity
check of stochastic training. The full posterior uses inverse precision,
which is different from reciprocal precision diagonals. Prior-predictive
coverage and fixed-signal coverage must both be reported without conflation.
Adding CTFs and our data preprocessing is a documented extension, not original
code reproduction. A baseline should retain a prior/regularization sweep
chosen without test outcomes and reference the original code's GPL license
if code is reused.

This preflight does not decide the revision scope; the independent review
will identify priorities. No extra baseline fit is claimed here.
