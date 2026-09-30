# Fixed-pose Fourier variational baseline

Primary source: Ullrich et al., *Differentiable Probabilistic Models of Scientific
Imaging with the Fourier Slice Theorem*, PMLR 115, 399--411 (2020),
https://proceedings.mlr.press/v115/ullrich20a.html. Sections 4.3 and 5.1 and
equations 14--15 were checked in the primary text. The repository tutorial
limitation is recorded separately in `REVISION-PREFLIGHT.md`.

The new implementation solves the paper's Gaussian variational objective
analytically in a specified linear, known-pose acquisition model. It is an
independent reimplementation of that objective, not a run of the original
software or reproduction of all its experiments. It does not learn poses,
noise or a flexible variational family. No third-party source code was copied.

## Parameterization and extensions

An odd 33^3 Fourier grid represents a real volume through independent Hermitian
coordinates. They are coefficients of the orthonormal real Fourier functions
1, sqrt(2) cos(2*pi*f*x), sqrt(2) sin(2*pi*f*x) on the unit cube. Corresponding
complex grid values are (x_cos-i*x_sin)/sqrt(2) and their conjugates. The odd
grid avoids self-conjugate Nyquist ambiguity. Observations interpolate the
complex grid trilinearly, then multiply by the published CTF and divide by a
declared common noise standard deviation. This interpolation is an approximation
to the continuous Fourier transform; its error is measured explicitly.

CTFs, an independent Gaussian-reconstruction pilot as prior mean, the unit-L2
normalization, linear Gaussian density targets and the explicit prior-scale
sweep are project extensions. The original paper's main experiments omit CTFs
and describe a standard Gaussian prior. Neither the absolute numeric prior
scale nor the target normalization is silently attributed to that work.

Let y=A x+e, e~N(0,I), and x~N(m0,tau^2 I). The posterior precision is
Q=A^T A+tau^-2 I. Minimizing the full Gaussian negative ELBO over a diagonal
Gaussian q gives

    mean_q = m0 + Q^-1 A^T(y-A m0),  var_q,j = 1/Q_jj.

The stationary variance identity follows by differentiating
(tr(Q diag(v))-sum_j log(v_j))/2. The mean solves the ordinary normal equation.
This analytic optimum includes both likelihood and KL terms; it does not repeat
the released tutorial's likelihood-only, near-zero-variance training loop.
For a target ell, solve Q v=ell and set w=A v. The full posterior target
variance is ell^T v; diagonal VI uses sum_j ell_j^2/Q_jj. Their means agree.
The former is not obtained by reciprocating only precision diagonals.

The independent check

    ||w||^2 + tau^2 ||A^T w-ell||^2 = ell^T Q^-1 ell

verifies prior-predictive calibration of the full Gaussian interval under its
matched model. Neither posterior promises uniform fixed-signal coverage over
the density/pose class of the robust audit. A diagonal approximation can be
wider or narrower for a particular target; it must not be described as
universally overconfident.

## Completed experiment

All three geometries use 1,024 distinct development particles and 220 independent
Fourier pairs per particle (radius 12). Noise is the supplied radius-5 simulation
scale held fixed across designs. Two physical targets have Gaussian sigma
10 Angstrom: a central average and the difference of averages centered at
z=+10 and -10 Angstrom. Prior expected squared deviation norms are 0.25, 1 and
4, specified as a sweep before outcomes. No prior is selected using coverage.

For each target and prior, controls include a signal generated in the same
interpolated Fourier model, the independent continuous 64-cell deposited-map
generator at nominal poses, and coherent/random boundary poses at 0.5, 1 and
2 degrees with 0.5 Angstrom translation budgets. Reported fixed-signal coverage
is analytic over the known simulated measurement noise. It is not observed
coverage on an unknown experimental molecule. Consensus-pose dependence and
experimental noise calibration remain unresolved.

All 18 solves completed in approximately 3.3--3.9 seconds per dataset for this
sparse conditional baseline. Each sparse acquisition matrix occupies 45.0 MB.
The continuous-reference forward discrepancy from trilinear interpolation and
finite representation is 16.8%, 16.2% and 25.3% across the stacks. This is a
material baseline-model limitation, not something to suppress in comparisons.
The variance identity residuals are at floating-point scale. Full posterior
prior-predictive coverage is 0.95. All tested diagonal-VI feature intervals are
wider than their full-posterior counterparts (approximately 1.07--1.42 times).
Poor fixed-map coverage therefore cannot be blamed generally on diagonal
covariance: prior shrinkage and representation/acquisition discrepancy matter.
The raw records retain favorable and unfavorable powers and coverage values.

Tests independently compare sparse realification with complex regular-grid
interpolation, Hermitian conjugacy and Parseval norm, and verify ELBO
stationarity and the full posterior solution against dense linear algebra.

## Holding the estimator fixed for a continuous audit

`audit_uq_fourier_variational.py` separately audits all 18 fitted weight vectors
in the original continuous cube at B=2. It explicitly includes the known
difference between the interpolated prior center and exact pilot center. If
that offset is kappa, the original estimator's fixed-pose bias upper bound is
|kappa|+B||ell-A* w||. A recentered comparison subtracts kappa from the estimator
and omits |kappa| from its bias bound; neither uses the reference outcome.

All audits completed. Half-widths are 0.257--0.815 of no data and 56--352 times
the full posterior's credible half-width. The fixed continuous-map cases are
covered conservatively, with numerically zero sign power throughout. Correcting
the known pilot offset does not restore power. These results expose the cost of
the supplied broad class while keeping the estimator fixed; they do not prove
that every physical prior, target or estimator must be uninformative.
