# Classical constrained Fisher direction: derivation and limits

Let d be the finite-training alternative-minus-null moment mean and U the
matrix of nonzero, normalized null-mean blocks. With power and bispectrum
blocks these columns have disjoint coordinates. Let S be positive definite.
For h satisfying U^T h=0, consider J(h)=(h^T d)^2/(h^T S h).

Define

    r = S^{-1}d - S^{-1}U (U^T S^{-1}U)^{-1} U^T S^{-1}d.

If U is empty, the second term is absent. This expression requires independent
constraint columns; zero columns are omitted. Then U^T r=0. Further,
S r=d-U c for a vector c, so for every feasible h,

    h^T d = h^T S r.

Cauchy-Schwarz in the S inner product yields J(h) <= r^T S r, with equality
for h proportional to r when r is nonzero. If r=0, all feasible mean gaps
are zero. Choosing h=r/||r|| orients toward the alternative because
r^T d=r^T S r>0. This also supplies an independently checkable solution:
if B spans ker(U^T), r=B(B^T S B)^{-1}B^T d.

In the experiment, S=.9*Cov_hat + .1*trace(Cov_hat)/p*I, using the unbiased
covariance of 65,536 noisy null-candidate features. All eigenvalues of Cov_hat
are nonnegative in exact arithmetic, and positive trace makes S positive
definite. The code solves systems by Cholesky rather than forming an inverse.
The numerical verifier uses the independent null-space expression above.

The two mean blocks scale as a^2 and a^3 because the power coordinates are
noise debiased and the nonredundant Gaussian bispectrum has zero additive
noise mean. Their separate annihilation therefore removes the finite-training
null mean at every scalar amplitude. It does not eliminate finite-view
training error, viewing-law shifts, amplitude effects on variance/tails,
or unknown noise/CTF errors.

This is an elementary equality-constrained Fisher/Rayleigh result, not a new
optimization theorem or proof of uniformly most powerful testing. The alternative
need not have the null covariance, and moments are not assumed Gaussian in the
calibration. The Rayleigh criterion does not optimize the actual rejection
probability or a worst-viewing-law criterion. The midpoint threshold is fixed
by training and may itself be inefficient.

The independent calibration guarantee is unaffected by the choice of training
procedure: condition on the training sample, score and threshold; apply the
existing probability-bound theorem to new calibration views and new test
images. Integrating this conditional inequality over training preserves the
stated size bound under the same image, amplitude, viewing and independence
assumptions. This is conditional composition, not experimental verification of
those assumptions. Selecting among calibrated variants using test outcomes
would require a separate multiplicity or selection argument; none is used.
