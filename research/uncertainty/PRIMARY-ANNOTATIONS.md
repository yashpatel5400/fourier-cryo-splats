# Primary-source annotations (living; targeted reading, not a full-screening count)

These notes record which claims and sections were checked, rather than treating
PDF retrieval as completed review. The source/hash manifest identifies the local
versions. General survey synthesis is in SURVEY.md.

## Xu, Balanov, Singer and Bendory: Bayesian orientation estimation

Source: arXiv:2412.03723v3 (23 February 2026), subsequently Acta Crystallographica D,
doi:10.1107/S2059798326001415. Checked introduction/contribution list, propositions,
heterogeneity experiment, discussion and code availability.

The estimator minimizes posterior expected chordal rotation loss; averaging
rotation matrices and projecting onto SO(3) is distinct from choosing a posterior
mode. The paper establishes a high-SNR connection to maximum likelihood and
studies downstream heterogeneous reconstruction with RECOVAR. Its Figure 6
compares ground-truth, MMSE and MLE pose inputs; this is a pose-estimation and
heterogeneity-quality experiment, not a coverage experiment for density CIs.
The discussion explicitly identifies individual-rotation confidence regions as
future work. Therefore its existence does not supply empirically validated pose
balls for the present conditional theorem.

Author code: https://github.com/AmnonBa/bayesian-orientation-estimation.
The currently documented toolkit is MATLAB, with Image Processing Toolbox and
3D alignment/MMSE demos. Reproducing a 3D demo is not equivalent to reproducing
its full single-particle 2D-projection experiment. Inspect the actual released
components before describing any run as that baseline.

## Rangan et al.: joint likelihood soft modes

Source: arXiv:2411.13263v2 (24 November 2024). Checked formulation, Hessian blocks,
soft-mode examples, discussion, and appendix contents including translations and
noise marginalization.

This is a direct nearest neighbor for pose--volume ambiguity and for asking
which extra orientations/defocus values could resolve structural uncertainty.
The likelihood representations invoke a low-temperature approximation or a
noise-marginalized regime. Local soft eigenvectors characterize sensitivity,
not a uniform fixed-parameter confidence guarantee under arbitrary nonlinear
pose errors. Several illustrative calculations are low resolution, and the
heterogeneity example describes a uniform view distribution and CTF equal to one.
Do not infer the scope of its demonstrations merely from the broad formulation.
The discussion explicitly anticipates ill-posed or spurious multi-particle
heterogeneity even when a related single-particle problem is well posed.

## Ullrich et al.: probabilistic Fourier-slice reconstruction

Source: PMLR 115, UAI proceedings published 2020; local arXiv:1906.07582v2.
Checked abstract/contributions, probabilistic model, uncertainty/variance and
model-bias discussion, and pose-inference limitations.

This work already connects differentiable Fourier reconstruction to explicit
density uncertainty, missing-view diagnostics and held-out model fit. Its joint
pose inference can encounter symmetry-related local optima. A diagonal Gaussian
mathematical control in this project reproduces a covariance approximation under
a matched linear model, not the paper's entire learning pipeline. The original
paper recognizes bias as a separate issue, so presenting it as claiming that
posterior variance automatically covers representation bias would be a straw man.

## Armstrong and Kolesar: inference over convex classes

Source: Econometrica 86(2), 655–683 (2018), doi:10.3982/ECTA14434.
Checked setup, fixed-length intervals, folded-normal critical values, and the
scope-for-adaptation discussion.

The central bias-aware Gaussian construction used here is prior art. Uniform
coverage over a broad class limits data-adaptive tightening toward smoother
subclasses. This is directly relevant to proposals to estimate a convenient
small density radius from the same noisy residuals. A numerical Gaussian basis,
a new name for its confidence bound, or an independent conic solver does not
create a new statistical theorem. Our possible contribution must be the cryo-EM
physics, practical audit, and convincing validation, with explicit scope.

## LocScale-2.0: confidence for enhancement

Source: Nature Communications (2026), doi:10.1038/s41467-026-75327-8.
Checked model/uncertainty calibration and validation descriptions.

MC dropout variances are empirically recalibrated and used for voxel scores
relative to a phase-conservative baseline map. Its target is the map-processing
comparison defined in the paper. The proxy should not be renamed unknown true
density, and the method should not be penalized for failing a different estimand
it did not promise. The distinction suggests a useful future matched experiment:
can ambient particle evidence validate or reject enhancement-induced features?
That experiment has not yet been implemented here.

## Beckers, Palmer and Sachse: significance maps

Primary full text: https://journals.iucr.org/d/issues/2020/04/00/rr5194/.
Checked signal detection, multiplicity, simulations and limitations.
The method tests voxel intensity against estimated background, with dependence-
robust multiple-testing options. Its FWER experiment trades detection sensitivity
for fewer false water-like peaks. The authors explicitly discuss misalignment,
preferred views and reconstruction artifacts that can be significant without
representing correct structure. This is therefore neither an unacknowledged
failure discovered here nor a method promising density-error coverage. A useful
comparison should separate background detection from robustness to a declared
forward-model error. Noise-region selection and homogeneity are material
assumptions; smoothing the displayed confidence surface is a visualization step.

## Sorzano et al.: bias, consensus and validation

Primary full text: https://journals.iucr.org/paper?ic5116=,
doi:10.1107/S2059798322001978. Checked formulation, sources of bias, FSC and
independent-half discussion, conclusions; supporting experiments not yet audited.
The paper directly warns that agreement can persist under shared processing
bias, and advocates examining parameter estimates across runs or algorithms.
It also notes that coincident estimates can share bias. Its claims about the
relative importance of parameter error are tied to its experimental examples;
our fixed-pose linear controls still exhibit a conventional bias--variance tradeoff.
We should cite this work when motivating the ambient audit, not present shared-
bias failure of half maps as a new discovery. Parameter agreement is a diagnostic,
not an externally calibrated bound on pose error.

### Released Ullrich code audit

Inspected the author repository at commit
`013149a0927183b5969095cb774d43b387644771`; file hashes are recorded in
`ullrich-code-manifest.json`. The README lists Python 3.6/PyTorch 1.1 and GPL-3.0.
The notebook includes generative projections and reconstruction training. Its
shown training cell minimizes sampled negative image log likelihood, with the
log scale initialized at -100; that cell does not include the prior/entropy
terms required for the paper's full variational objective. Thus merely running
that tutorial cannot be labeled a faithful complete uncertainty reproduction.
The project currently uses explicitly identified mathematical covariance
controls. No external-code reproduction is claimed from this inspection.

## Batlle et al.: simultaneous constrained inverse inference

Source: https://arxiv.org/abs/2510.11708v1, primary full text checked for setup,
single-versus-multiple functional distinction, strict bounds, test inversion,
and the Bonferroni comparison. The paper develops constrained finite-sample
regions for multiple linear functionals and distinguishes their joint geometry
from a product of marginal intervals. A generic data-consistency region can
cover all functionals simultaneously; applying a second Bonferroni adjustment
to projections of that same region is unnecessary. Its sharper test-inversion
regions go beyond the elementary correction appropriate to separate marginal
certificates. We cannot present simultaneous feature inference or test inversion
as an unexplored contribution. Our present nonlinear pose treatment uses an
explicit bounded perturbation relaxation; their principal forward-model setup
is linear with structural parameter constraints.

## Experiment-guided AlphaFold3: the atomic-ensemble target

Source: https://www.nature.com/articles/s41587-026-03166-5, published 29 June
2026. Checked cryo-EM examples, inputs, and ESP forward/guidance models. Cryo-EM
inputs are reconstructed ESP maps plus sequence/reference-model information;
the guidance uses density fit and optimal transport, with separate alignment.
This is atomic ensemble inference conditioned on maps, not raw-particle density
coverage. The multimodal examples show that a better map fit need not imply a
better fit to independent NMR constraints. This supports keeping measurement
fit, structural ensembles and calibrated reconstruction uncertainty distinct.

## CryoDiff: evidence-access boundary

The primary bioRxiv abstract at
https://www.biorxiv.org/content/10.64898/2026.06.04.730282v1 was retrieved through
search on 29 September 2026. It describes diffusion map enhancement and voxel
confidence from Monte Carlo sampling. Direct full-text access returned HTTP403.
The full calibration target and procedures remain unverified; do not equate
its confidence score with frequentist raw-particle density coverage or claim
to have critically reviewed/reproduced the full method from this abstract.

## EMReady2: local quality conditioning is not density interval calibration

Primary publisher text checked on 29 September 2026:
https://doi.org/10.1038/s41467-026-71794-1 (16 May 2026 publication; publisher also
lists a 17 July version-of-record date). Read the input/output description,
local-resolution target construction, evaluation and discussion limitations.
The model enhances a supplied map, training against atom-derived maps with
locally varying blur informed by Q-scores. Its 136-map evaluation uses map/model
and interpretation metrics. This is not a raw-particle density confidence
interval procedure. The authors explicitly distinguish real-space map defects
from conformational/compositional variability and discuss learned-prior
hallucination risk. Therefore we should not label it an uncertainty baseline
merely because training is local-quality-aware. A matched use would test whether
particle evidence supports selected enhanced-map features, with target selection
separated from evaluation. That experiment is not yet executed here.

## Stochastic inverse reconstruction over conformational distributions

Primary full text: https://arxiv.org/html/2509.05541v1, checked 29 September 2026.
The page labels arXiv v1 (2025), while its displayed manuscript date is August
2026; preserve the exact URL/hash rather than inventing a v2 publication.
Checked the random forward model, variational formulation, Proposition 5,
particle-discretization lemma and conclusions. The estimand is a distribution
of structures whose random projections reproduce an image distribution. The
paper relates distributional optimization to MAP discretizations and gives
conditional consistency statements using continuity, compactness and empirical
convergence assumptions. This is distinct from finite-sample coverage of a
fixed homogeneous density functional. Its conclusion identifies experimental
data benchmarks and fuller particle-method convergence analysis as future work.
Continuous/discretized distinctions are already prior art; our proposed
continuous adjoint norm is a specific uncertainty audit, not the first
continuous-space formulation of cryo-EM.
