# Classical covariance-aware score comparison, version 1

Declared after inspecting all candidate-score-v1 outcomes, before producing
any new covariance-training, calibration or held-out outcome. This is a
classical comparison, not a claim of a new Fisher discriminant or optimal
robust testing procedure. No previous calibration/test draws are reused.

Use the same three fitted Gaussian candidates, normalization, physical-cell
operator, inherited frequency/transfer fixture and candidate-only 20-Angstrom
region rule. Generate 65,536 fresh Haar training views per stack. At each view,
retain the noiseless candidate and 25%-deleted Fourier means and one fresh
proper unit Gaussian noise sample for the unchanged candidate. Compute exact
conditional moment means and the unbiased sample covariance of the noisy
candidate features, which includes both viewing and measurement variation.
Accumulate centered second moments in batches; save both mean arrays and the
noise array to independently replay this calculation.

For each of power-only and power-plus-512-bispectrum features, compare:
(1) the Euclidean matched contrast projected off the two candidate-mean blocks;
(2) a constrained Fisher direction using S=.9*sample_covariance +
.1*trace(sample_covariance)/dimension*I. There is no shrinkage selection.
Both use the same 65,536-view mean gap and nuisance constraints. For each,
orient toward the 25%-deletion alternative, normalize to Euclidean norm one,
and use the midpoint of training null/alternative score means as threshold.
A direction of norm at most 1e-12 is zero and retained. The two candidate mean
blocks enforce finite-training-law amplitude cancellation; no continuous-law
zero identity is claimed. Alternative features need not be Gaussian or share
null covariance: the Rayleigh criterion is only a design criterion.

Calibration uses 32,768 new Haar views, two independent 32-noise groups,
16 amplitude cells on [.9,1.1], center .5 and confidence failure .001.
Testing uses 131,072 new Haar views and all deletions 0/.1/.25/.5/1 and
amplitudes .9/1/1.1. Retain all seven existing probability-bound methods,
kappas 1/1.01/1.1/2/5 and particle counts 1,000/10,000/100,000, with test
level .049. Also retain actual outcomes of 128 disjoint groups of 1,024.
Do not combine methods or score variants by an unadjusted minimum.

Seeds are 261017+int(dataset)+stage*1,000,000 for training/calibration/testing
stages 0/1/2. Batch sizes are 512/32/512. Both directions, feature families,
alternatives and amplitudes share draws within a stage. Save all RNG states,
input/source hashes, directions, thresholds, calibration envelopes, scores,
counts, numerical timings and direct-cell first-view checks. Retain every
attempt. There are 18,900 scalar projections and 6,300 group-outcome cells.

Correct-null controls and every stack are mandatory. Independently verify the
covariance and the constrained optimizer, rather than checking only its own
formula. Improvement would diagnose a weak original score, not establish
experimental calibration, unconditional calibration-repeat error, regional
attribution or novelty. No new full-review acceptance claim follows.

Classical positioning: Friedman's Regularized Discriminant Analysis technical
report, Stanford SWI NSF 17 (October 1987; journal version 1989):
https://statistics.stanford.edu/technical-reports/regularized-discriminant-analysis
The primary software documentation describes covariance-aware LDA and
shrinkage: https://scikit-learn.org/stable/modules/lda_qda.html . Metadata and
the documentation's mathematical/shrinkage sections were consulted; the full
Friedman report has not been read in this protocol. The block constraints are
a direct equality-constrained Rayleigh optimization, with proof recorded in
the study, not a novelty claim attributed to those sources.
