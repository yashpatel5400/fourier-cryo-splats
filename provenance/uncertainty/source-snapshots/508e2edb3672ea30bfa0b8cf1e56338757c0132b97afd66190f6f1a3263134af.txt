# Conditional application to experimental particles

Prospective amendment, before reading the experimental feature estimates for
the locked targets: apply all 12 fixed-weight estimators and their completed
pose audits to the actual cached particle images. Retain all four pose classes
(fixed, shift-only, one degree, two degrees), all features and all failures.
The weights and selected targets are already specified by the pilot-only lock;
do not select them using the calibration or inference outcomes.

Use the scalar pilot amplitude from the earlier grouped experimental attempt,
which was fitted on the separate old pilot pool. Calibrate the variance on the
second half of that attempt's fresh exposure representatives, now at radius
12. The first half was used in other exploratory method development. The second
half has also been inspected in previous development experiments, so this is
transparent exploratory reuse, not new confirmation. Verify its exposure groups
are disjoint from every old pilot/inference group before accessing observations.

The new 128-particle estimators can include multiple particles per exposure.
Allow arbitrary jointly Gaussian dependence within each inference exposure;
retain independence across exposures and one common marginal raw-coordinate
noise covariance Sigma across inference and calibration. If exposure g has
n_g particles and scalar noise contributions X_i = v_i^T epsilon_i, then

    Var(sum_{i in g} X_i) <= (sum_{i in g} sqrt(Var X_i))^2
                           <= n_g sum_{i in g} v_i^T Sigma v_i.

Therefore the total estimator variance is at most tr(V_g Sigma V_g^T), where
row i of V_g is sqrt(n_g)*v_i. Apply the existing noncentral-Gaussian trace
lemma to V_g times each independent calibration observation. This elementary
group-size inflation removes the extra within-exposure independence assumption;
it does not establish common covariance, Gaussianity or between-group independence.
Without the inflation the comparison would be unjustified for these estimators.

Pull weights back to uncentered Fourier coordinates, including the supplied
centering phase and prescribed whitening factor. Calibration observations are
raw image Fourier samples divided by the independent pilot amplitude. Allocate
beta_calibration=0.005/12, delta_spectral=0.000001/12, and
alpha_noise=(0.045-0.000001)/12 to each feature. The total is 0.05/12, so
Bonferroni gives conditional simultaneous coverage across the 12 features for
each declared pose class. This extra calibration budget differs explicitly
from the known-noise development calculation. The same calibration event and
spectral event can be reused for deterministic remainder refinements.

Report observed centers, calibrated standard-deviation bounds, bias terms,
interval widths, fallback, zero exclusion, and the amplitude-aligned deposited
map's feature and inclusion indicator. Also report that map's continuous L2
distance from the pilot to check the declared class for this comparator.
The map comparison is not a coverage measurement. Uniform pose radii, density
radius, a single shared density and independence of supplied consensus poses
remain unverified. EMPIAR-10076 is explicitly heterogeneous; arbitrary
calibration means do not make the inference theorem valid for its heterogeneous
particle population. Do not present this attempt as resolving those limitations.
