# Covariance-guided design with independent final calibration

The original grouped experimental estimator used a scalar trace envelope to
choose weights. The independent directional audit then held those weights
fixed. It improved noise precision but could not undo their conservative fit.
This new exploratory study separates a useful design metric from a valid final
noise bound. It is not a newly frozen confirmation or a simultaneous claim.

Using only the original first fresh-exposure half, define the positive metric
C = 0.8*(sum_j Y_j Y_j^T/n) + 0.2*tr(sum_j Y_j Y_j^T/n)/m * I.
Means are not subtracted. This empirical second moment may contain signal; it
is merely a fitting guide. The shrinkage fraction 0.2 is chosen in the script
before running the new fit and is not tuned on the final calibration pool.

If O_i is the orthogonal realification of particle i's supplied centering phase,
let W be block diagonal with blocks (O_i C O_i^T)^(-1/2). Write physical weights
as w=W x. Minimizing z||x||+B||ell-A*W x|| is then the original continuous
fixed-pose convex design problem with target vector W A ell and Gram W G W.
The numerical preconditioner transforms the existing low-rank Gram factors;
it does not discretize the unknown density. Quadrature residual errors are
unchanged after evaluating at W x. The Gram-action error is multiplied by an
upper bound ||W||=lambda_min(C)^(-1/2) for the transformed dual diagnostic.

The proxy objective's noise norm and optimization gap refer ONLY to C. They
are not noise confidence bounds. After selecting w, pull it back into the raw
Fourier frame and apply the independent second-half directional calibration
lemma. Report the resulting actual calibrated noise upper together with the
full continuous density/pose sensitivity audit. Preserve no-data fallbacks and
all approximate-reference comparisons. A better proxy objective need not
produce a better final interval, and an estimated covariance alone cannot
justify this replacement of the final audit.

Dense independent tests check the real centering transforms, matrix powers,
colored-noise variance identity, transformed Gram, Woodbury inverse and
quadrature-error transformation. Common Gaussian covariance, independence
between exposures and supplied poses, and the density/pose class remain
unverified on experimental data.
