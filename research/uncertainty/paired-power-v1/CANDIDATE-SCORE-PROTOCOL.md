# Candidate-only local-perturbation validation gate

Post-outcome method development after v0.7.4-dev, motivated by the unresolved
oracle-score criticism. No new positive three-stack outcome is presumed. Use
the already fitted `results/final/{10028,10049,10076}/gaussian-mean.mrc` maps,
not an EMDB reference, for score/region design. These maps were fitted to
8,192-particle stacks using supplied consensus poses. Their history does not
provide an experimental uncertainty guarantee. All new observations are fresh
simulations independent of those fits.

## Fixed design from the candidate alone

Interpret each saved 64-cubed map as physical cells on the unit cube, normalize
its Euclidean coefficient norm to one, and use its saved pixel size to set the
field in Angstroms. The signed map is used without positivity clipping. Locate
the largest Gaussian-smoothed density value within the central half-field cube;
the smoothing standard deviation is 20 Angstroms, zero extension and a four-
standard-deviation discrete filter. Use that peak as the center of a Gaussian
regional mask of the same physical width. This prespecifies one region per
candidate; it does not search the new test images or an external reference.

Draw 4,096 fresh continuous Haar training views. Use the same 220 frequencies,
512 frequency-selected triads and first known transfer/noise profile as the
preceding studies. Evaluate both the candidate and its masked region with the
physical-cell Fourier operator. At each moment family (power or combined),
form the mean feature difference between a 25% regional deletion and the
unchanged candidate. Compare the normalized difference and a scale-orthogonal
version. The latter removes its projection on the training mean power vector
and, separately, the mean bispectrum vector. These blocks scale as a^2 and a^3;
their removal makes the candidate's training-law mean contrast zero at every
global amplitude. This identity is for the finite training mean, not a
continuous-view validity assertion. The contrast is a classical matched
direction, not a novel score-training principle or a discovery of unknown
errors. Below norm 1e-12, preserve a zero score as uninformative.

Set each threshold halfway between its training mean under the candidate and
the 25%-deletion mean. These four scores and thresholds then remain fixed. The
score never receives the true held-out deletion amount, image noise or an
external reference. The hand-specified deletion family remains a strong
inductive assumption; a test for one nominated local error is not an omnibus
test of map accuracy.

## Independent calibration and held-out controls

For each of all three candidates, use 32,768 fresh Haar calibration views,
two independent groups of 32 noises per view, 16 amplitude cells on [.9,1.1],
and center b=.5 fixed without viewing the calibration. Calibrate all four
scores with the same three original procedures (individual ratio, grouped
ratio, paired viewing variance), alpha=.05 and delta=.001. Physical amplitudes
are independent of noise conditional on view. Kappa values are 1/1.01/1.1/2/5.
Each method is a separate test: no uncorrected selection among the four scores
or three calibration procedures is used to make an overall decision.

Generate 131,072 independent Haar held-out views. At each view retain deletion
amounts 0/.1/.25/.5/1 and amplitudes .9/1/1.1, using shared fresh proper Gaussian
noise across every method/map/amplitude. The unchanged candidate at all three
amplitudes is the correct-null control. Full deletion is no longer the sole
alternative. Deletion means are computed by linearity from the candidate and
region projections, with a direct physical-cell check at the first view of
every phase. Known noise and transfer scaling are held fixed; neither the
candidate normalization nor the simulator is claimed experimentally calibrated.

Use seeds 261016+int(dataset)+stage*1,000,000 for training/calibration/held-out
stages 0/1/2. Batches have 512 training/held-out views and 32 calibration views.
Preserve every score array, input/source hash, RNG state, center/mask,
training direction and threshold. Retain one projection row for each
score/kappa/procedure/n, storing all 15 map/amplitude outcomes in that row.
There are 540 rows total over three stacks and 8,100 scalar rejection outcomes.
Also partition held-out views into 128 disjoint groups of
1,024 particles for actual conditional rejection counts. Calibration is fixed;
these do not estimate unconditional error over calibration repetitions.

Report all three stacks, both feature families, both scale variants, every
deletion size and all correct-null controls. A favorable result would establish
only a candidate-driven simulator test under supplied noise/view/amplitude
assumptions, not experimental validity, three successful reconstructions,
omnibus detection or ICML-level novelty. Calibration repetitions, extremal
view laws and real-noise controls remain separate unresolved experiments.
