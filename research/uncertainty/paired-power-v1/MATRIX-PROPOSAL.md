# Matrix paired-exposure extension: development gate

The enlarged-view diagnostic defeats the diagonal cross-power statistic on
10076. Retaining correlations between Fourier coefficients is a different
hypothesis. It retains phase information but loses automatic common-translation
invariance. The initial screen fixes translations exactly; no unknown-pose
experimental guarantee or novel general Gaussian theorem is claimed.

## Conditional construction

After independently fixed whitening and an orthonormal projection, let X,Z be
real Gaussian vectors with a common mean m and independent errors. Their
separate covariance matrices satisfy S_X,S_Z <= I in Loewner order. Choose a
symmetric matrix T with all eigenvalues strictly between -1 and 1. Write

    U = T(I-T)^(-1),
    log E = X' T Z + (1/2) log det(I-T^2).

If m'U m <= 0 for every mean in the null, E has expectation at most one.
For a candidate map with arbitrary nonnegative particle amplitude, it suffices
to impose m(R,s)'U m(R,s) <= 0 at unit amplitude for every allowed rotation
and translation. Any arbitrary viewing distribution then satisfies the null.
A finite view catalog supplies only a weaker, discrete null.

At S_X=S_Z=I, diagonalizing T reduces the Gaussian integral to independent
scalar bilinear integrals. Its value is

    E exp(X'TZ) = det(I-T^2)^(-1/2) exp(m'U m).

For smaller separate covariances, the central determinant increases with each
covariance because T S_Z T and T S_X T are positive semidefinite. The
noncentral coefficient for the joint quadratic form H with off-diagonal
blocks T/2 is K(S)=(I-2HS)^(-1)H. Along any positive-semidefinite covariance
increment its derivative is 2 K(dS)K, positive semidefinite. The moment domain
holds on the path to the identity: the relevant singular values of
S_X^(1/2) T S_Z^(1/2) are below one. Thus the bounding identity covariance
upper-bounds both factors. This proves the stated conditional expectation
inequality. Independent products and Markov's inequality supply a level-alpha
test at 1/alpha. These are standard Gaussian and e-value arguments.

The two means must agree. Independent exposures do not by themselves guarantee
this in the presence of dose effects, relative motion or processing dependence.
The covariance upper bound is a matrix inequality, not a marginal noise SD.
Compression must be fixed without the inference exposures. An oracle mean
used to choose a test is legitimate only for this development power screen,
not an experimental test learned on the same particles.

## Cone diagnostic and its limit

Let S=E[m_alt m_alt'] denote the alternative's *uncentered* second moment.
The expected log factor equals tr(TS)+(1/2)log det(I-T^2). If
S=sum_l lambda_l m_l m_l' for lambda_l>=0 and null means m_l, then

    tr(TS) <= tr(US) <= 0,

because U-T=T^2(I-T)^(-1) is positive semidefinite. The log determinant is
nonpositive. Hence no member of this matrix bilinear family has positive
expected log growth for that alternative. Matching this moment does not
match the complete distribution and is not a general impossibility theorem.

Conversely, a symmetric separator H with m_l'H m_l<=0 for every catalog
member and tr(HS)>0 yields positive expected growth for sufficiently small
positive c, using U=cH and T=U(I+U)^(-1). Its derivative at c=0 is tr(HS),
while the log determinant has zero first derivative. The moment domain
requires lambda_min(cH)>-1/2. This is a finite-dimensional cone-separation
argument, not a new general identifiability result.

## Translation tradeoff

Suppose nonredundant complex frequencies are distinct up to sign and exclude
zero. A symmetric real bilinear statistic invariant to every common continuous
translation must commute with all block phase rotations. In complex notation,
off-diagonal conjugated terms carry q_j-q_k and unconjugated terms carry
q_j+q_k. Their coefficients must vanish for invariance. Each remaining real
2-by-2 block must be a scalar multiple of identity. Thus the invariant
statistic reduces to diagonal cross power. Repeated or opposite frequencies
would introduce additional blocks, so the nonredundancy premise matters.

This rules out claiming that the matrix extension both retains general
off-diagonal information and automatically eliminates arbitrary common shifts.
Continuous bounded-shift certification or separately defensible centering
would be an additional methodological requirement.

## Declared first screen

Reuse all 4,160 saved orientations for the alternative as well as the null;
this avoids exploiting the previous alternative's 64-view rank. Use each
stack's first saved CTF/noise profile. Its known CTF signs can be removed by
an orthogonal diagonal change of coordinates, so the saved squared transfer
suffices. The three stacks are 10028, 10049 and 10076. Compare the true-map
control and full-removal candidate at ranks 8,16,32: eighteen cells total.

For each stack, choose a shared orthonormal basis from the top 32 eigenvectors
of the average true/full-removal real mean second moments. This is an oracle
compression; report retained signal energy. Nest ranks using this same basis.
For each cell, solve a nonnegative covariance-cone approximation using a
coordinate-scaled infinity residual. Off-diagonal symmetric coordinates use
sqrt(2) scaling, preserving the Frobenius inner product. Normalize each null
outer product to trace one and each moment equation by max(abs(target),
1e-3 max(abs(target))). A 60-second HiGHS solve limit applies to each cell.
All timeouts/failures are retained, not treated as nonmembership.

The LP dual supplies a separator. Repair its maximum finite-view violation
by subtracting that value times identity, with a 1e-12 numerical margin.
Along this one direction, evaluate zero and 181 log-spaced positive scales
within the moment domain, then locally refine the best adjacent interval.
Report feasible expected-growth values at noise covariance bounds 1,2,4;
these are lower values for this sampled-view/oracle experiment, not a global
optimum or continuous-pose certificate. Save weights, basis, moments, mixture
coefficients, source hashes, solver metadata and all eighteen outcomes.
No new observations, rejection-power trials or experimental claims occur.

The initial unit-test run exposed an advanced-indexing shape mismatch in the
feature matrix (four tests passed, two failed before reaching the LP). The
indexing was corrected before any molecular screen. Subsequent test output
is retained separately. Tests check implementation, not physical calibration.

Closest prior art remains the moment-comparison and SubspaceMoM sources in
`../MOMENTS-AND-IDENTIFIABILITY-READING.md`. This proposal has not established
a novel practical contribution beyond those methods.
