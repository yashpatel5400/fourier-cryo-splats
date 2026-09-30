# Independent calibration of a feature estimator's noise variance

This is a post-review development refinement of the elementary trace bound,
not a new general concentration theorem. It changes only the noise audit of
existing fixed estimators. Supplied pose/density bounds, Gaussianity, exposure
independence, and covariance stability remain assumptions.

## Conditional statement and proof

Let the uncentered realified Fourier noise in every inference and calibration
particle have one common covariance Sigma. Particle noises are independent;
calibration vectors Y_j may have arbitrary deterministic, nonzero means.
Condition on training data and the design. Suppose the feature estimator's
inference noise is sum_i v_i^T epsilon_i, with fixed raw-frame weight rows v_i.
Set V to the matrix with these rows. Its variance is

    s² = sum_i v_i^T Sigma v_i = tr(V Sigma V^T).

The projected calibration vectors V Y_j have common covariance V Sigma V^T
and arbitrary means. The earlier trace lemma therefore applies, giving

    s² <= S_V / (n t),    S_V = sum_j ||V Y_j||²,

with probability at least 1-beta, where
n*(t-1-log(t))/2 = log(1/beta). Correlations between the components of V Y_j
are allowed. The rank or number of inference particles need not be smaller
than the number of calibration samples. This result does not imply a matrix
upper confidence bound for arbitrary adaptively chosen weights.

For the SAME calibration pool and beta, S_V <= ||V||_F² sum_j ||Y_j||².
Thus the direct feature-variance bound cannot exceed the original trace bound
times the estimator's squared weight norm. The experimental comparison uses
different independent pools for weight selection and this audit, so that
pointwise sample inequality alone does not promise improvement in every case.

Let |b| <= b_plus be an independently valid deterministic bias bound. On the
calibration event, Gaussian estimator noise has standard deviation at most
s_plus. Using the folded-normal critical half-width at alpha_noise with
(b_plus,s_plus) is valid: its width is at least b_plus, and for |b|<=b_plus
its outside probability increases with noise scale. Independence of inference
noise and calibration then gives failure at most beta+alpha_noise. This is a
pointwise feature interval. Simultaneous statements need a multiple-testing
budget; selecting targets after inspecting inference values is not covered.

## Raw-frame centering and independence

If the centered observation is p_iq Y_iq/sigma_old, |p_iq|=1, and the centered
realified weight is w_R+i*w_I, the corresponding raw coefficient is
conj(p_iq)*(w_R+i*w_I)/sigma_old. The implementation checks the complex/real
inner-product identity. Applying the same weights directly to uncentered
calibration vectors would generally be wrong for colored noise.

The original grouped experiment used the first half of fresh exposure groups
to select a conservative scalar whitening factor and fit weights. This new
audit holds those weights exactly fixed and uses ONLY the second, independent
half to calibrate their variance. The second half previously supplied aggregate
energy diagnostics; this is transparent exploratory reuse, not newly frozen
confirmation. No signal subtraction or mask enters the calibration statistic.
Pilot amplitudes are selected using the separate old pilot pool. Any dependence
introduced by published full-data consensus poses remains unresolved.

## Verification and scope

The tests independently verify the stacked-covariance variance identity,
raw/centered Fourier estimator equivalence, domination for a common calibration
sample, and the exact noncentral chi-square lower tail for correlated rank-one
noise. The experiment reconstructs the original calibration statistic and
original interval widths before making the new comparison. Source weights,
source-result hashes, exposure indices, all outcomes, and failures are retained.
Neither narrower intervals nor agreement with an approximate deposited map
establishes the unverified conditions or true density coverage.
