# Experimental calibration attempt after round 1

This development experiment addresses part of R1/R8; it does not establish the
missing end-to-end experimental guarantee. The frozen original studies remain
unchanged. No new particle bytes are required: both archive cohorts are already
downloaded and hashed.

## Conservative noise calculation with nonzero signal

Suppose independent calibration vectors obey Y_i~N(mu_i,Sigma), with a common
unknown covariance and arbitrary fixed means. Let T=tr(Sigma), S=sum ||Y_i||².
For lambda>0 the noncentral Gaussian Laplace transform gives

    E exp(-lambda S)
      <= det(I+2 lambda Sigma)^(-n/2)
      <= (1+2 lambda T)^(-n/2).

The second inequality uses the product expansion over nonnegative eigenvalues.
Applying the exponential Markov inequality and optimizing lambda gives

    P(S <= n t T) <= exp(n*(1-t+log(t))/2),  0<t<1.

Choose t so that this upper bound equals beta. Then T_plus=S/(n t) is an upper
confidence bound for T and Sigma<=T_plus I on that event. This is an elementary
Gaussian concentration calculation, not new general statistical theory. It
permits signal contamination but can lose a factor comparable to dimension
relative to a known white-noise variance. Unknown group-specific covariances or
dependent calibration samples are not covered by the statement.

For inference vectors independent of calibration and of the fixed design,
whitening by sqrt(T_plus) makes the stacked noise covariance at most I if the
particle noises are independent. The usual norm bound on the scalar estimator
noise then holds. A bias-aware interval at alpha_noise and calibration beta has
failure probability at most alpha_noise+beta. The invariance of S under separate
orthogonal transformations does not establish independence for data-fitted
poses; that remains a separate limitation.

## Data use and declared sensitivities

`probe_uq_experimental_noise.py` fits the scalar amplitude of an existing
unit-norm pilot on old pilot-pool particles. It selects one particle from each
fresh exposure group and splits groups into calibration and diagnostic halves.
Calibration uses unmasked Fourier coordinates, without subtracting a fitted
signal, using alpha_noise=0.045 and beta=0.005. The revised grouped attempt also
selects only one inference particle per old inference exposure group.

The initial development attempt had 128 inference particles, including repeated
source groups. It is retained under `experimental-noise-attempt`; it additionally
requires the unverified within-exposure independence assumption. The grouped
attempt avoids that particular replication, but cannot prove independence
between groups or a common noise covariance from the same archive.

The targets are central Gaussian density averages with sigma=0.07 of the field
and sigma=10 Angstrom. Neither is selected based on the observed inference
values. Supplied density radius B=2, rotations 0/1/2 degrees, and nonzero-pose
translation radius 0.5 Angstrom define a sensitivity table. They are not learned
confidence bounds. Both the center and width switch to the pilot-only interval
if the audited interval exceeds the no-data bound.

The deposited map is an approximate comparison after an amplitude fit on the
pilot pool. Its feature value is not a coverage label. Consensus poses were fit
using broader datasets, original extraction/preprocessing is not reversed, the
density energy radius is not physically calibrated, and noise-model stability
remains empirical. These unresolved conditions must accompany any resulting
interval. A table of numbers alone cannot resolve R1.

## Explicit structural assumption clarification

The inference theorem also assumes a shared density in the declared class.
This is not established for the experimental stacks. The primary archive
identifies [10076](https://empiar.ipr.pdbj.org/en/entry/10076/) as a heterogeneous
mixture of L17-depleted 50S assembly intermediates with preferred orientations.
The arbitrary means allowed in CALIBRATION do not extend the INFERENCE model
to arbitrary particle conformations. An ensemble-mean interpretation needs a
justified relation between state, viewing direction and the sampling operator;
one deposited class map cannot establish it. The simulation violation controls
already demonstrate failures for view-coupled mixtures. Original numerical
records are retained unchanged, with this additional scope clarification.
