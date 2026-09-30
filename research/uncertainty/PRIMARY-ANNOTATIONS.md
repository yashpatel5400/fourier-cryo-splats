# Primary-source annotations (living; targeted reading, not a full-screening count)

These notes record which claims and sections were checked, rather than treating
PDF retrieval as completed review. The source/hash manifest identifies the local
versions. General survey synthesis is in SURVEY.md.

## Cai and Low: adaptive confidence intervals (revision reading)

Primary electronic reprint: arXiv:math/0503662v1, published Annals of Statistics
32(5), 1805–1840 (2004), DOI 10.1214/009053604000000049. Checked introduction,
ordered/between-class modulus definitions, Proposition 1 and Theorem 1 with their
proofs, and Remark 1 (PDF pages 1–7). This is targeted reading, not a claim to
have checked all 37 pages. Retrieval bytes are hashed in
`revision-foundations-manifest.json`.

Coverage over a union of classes constrains expected length even at parameters
in an easier subclass. Their two-point Gaussian argument and between-class
modulus quantify that obstruction. Their framework includes variable-length
adaptive intervals; our fixed-length design audit is much narrower. The current
two-point proposition and generic efficiency observation cannot be presented as
new general confidence-interval theory. A pilot that appears smooth does not by
itself license replacing the advertised density class by a smaller one.

Low (1997), DOI 10.1214/aos/1030741084, was also requested. The primary publisher
download returns a security page instead of a PDF. We have not critically read
that full text and do not count it as such. The directly read Cai--Low source
supports the adaptation/expected-length discussion used in this revision.

## Kuczynski and Wozniakowski: random-start spectral estimation

Verified published record: SIAM J. Matrix Anal. Appl. 13(4), 1094–1122 (1992),
DOI 10.1137/0613066. Downloaded the authors' March 1989 Columbia technical report
CUCS-465-89 from the university host, 37 scanned pages. Rendered and OCR-read
pages 1–5 (abstract, introduction and start of problem definition); this is not
a full proof audit of the report or an assertion that it is identical to the
published version. Its hash is in `revision-foundations-manifest.json`.

The report analyzes average and probabilistic power/Lanczos errors for random
starts and emphasizes the difficulty of reliable stopping. A deterministic
start can miss an unseen eigendirection. We use this as prior art for randomized
norm estimation, with our short Gaussian-projection upper-bound derivation
stated separately in full. We do not call that construction a new general
spectral-estimation theorem or use a Ritz value as an upper bound.

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

## DynaMight: validation of estimated deformation fields

Primary published PDF, doi:10.1038/s41592-024-02377-5; checked training,
model-bias experiments, half-set reconstruction and deformation-error discussion.
Separate variational autoencoders and consensus Gaussian models are fitted to
halves. A reserved particle subset is embedded through the separate models;
disagreement of the resulting displacement estimates is an error diagnostic.
The paper explicitly demonstrates that atomic-model regularization can transmit
incorrect structural features, and warns that a shared deformation model can
inflate half-map agreement. These concerns are prior art, not failures newly
identified by our audit. The diagnostic targets motion estimates and depends on
its split/regularization design; it is not a uniform confidence interval for a
fixed homogeneous density functional. Its comparison also reports cases where
multi-body refinement has better local resolution. Our validation should retain
similarly unfavorable comparisons rather than interpreting conservatism as
universal superiority.

## 3DFlex: motion learning and frequency-separated validation

Primary publisher text: https://doi.org/10.1038/s41592-023-01853-8. Checked forward
model, training/refinement split and experimental validation description.
The method deforms a canonical density through a mesh-based neural flow,
assuming supplied poses and CTFs in the described experiments. Low-resolution
images train the motion model; fixed learned motions then support full-resolution
half-map reconstruction. Latent-coordinate noise is a training regularizer,
not a calibrated interval for the physical density. The authors' validation
argument relies on evaluating recovered information beyond motion-training
frequencies. DynaMight raises an additional shared-bias concern; the appropriate
comparison must preserve both methods' actual protocols rather than asserting
that any half split is automatically independent. This is a heterogeneity and
resolution baseline, not a matched homogeneous confidence-coverage method.

## RECOVAR: uncertainty of latent coordinates and population deconvolution

Primary final article: https://pmc.ncbi.nlm.nih.gov/articles/PMC11892586/,
doi:10.1073/pnas.2419140122. Publisher/PMC-indexed method and discussion passages
were retrieved; direct PMC access returned a browser challenge. The separately
available PMC10634927 version is a preprint and must not be mislabeled the final
publication. The method estimates conformational covariance, embeds images in a
linear principal-component space, models noisy latent coordinates, and uses
adaptive kernel regression to reconstruct states. It then corrects the latent
population density for the observation-noise effect. This separates covariance
of physical structures, uncertainty of per-image coordinates, and populations;
none should be equated automatically with a confidence interval for a single
homogeneous density. Its linear embedding permits statistical corrections that
a general neural latent coordinate system need not preserve. Our homogeneous
scope does not reproduce its population task, so a simplistic one-number
coverage ranking would be misleading.

## Evans et al.: counting and population inference

Primary published Comment: doi:10.1038/s42003-026-09859-6. Checked two-state and
continuous examples, the real-data mixture construction and scope discussion.
Particle assignments require the ensemble's base rates; simply counting hard
or soft image assignments can distort populations as noise increases. The
comparison uses whole-dataset reweighting/deconvolution to address that inverse
problem. Its real spike example uses a constructed 80/20 mixture of heavily
classified subsets and explicitly acknowledges that this baseline is not known
biological truth. In the continuous example, supplying true candidate structures
to a comparator is also disclosed. Thus population uncertainty deserves a
separate benchmark with candidate-structure error and pose/contrast uncertainty;
our homogeneous two-state assumption-stress test is not a reproduction of this
population-estimation result.

## Mattingly et al.: choosing a learnable conformational resolution

Primary arXiv:2606.14449v1. Checked model, mutual-information objective, Gaussian
information approximation, RNA example and its conditional imaging parameters.
The method chooses a representative conformational ensemble using expected
information about its mixture weights under a Dirichlet prior and probabilistic
image model. Its RNA example uses a Gaussian/Fisher-information approximation,
conditional on specified imaging parameters, to compare ensemble sizes and
selection schemes. The quantity is prior-averaged learnability of populations,
not a uniform finite-sample confidence interval for arbitrary density. The
paper directly studies indistinguishable nearby states and measurement-induced
coarse graining. Therefore our fine-feature uncertainty cannot be claimed as
the first theory of measurement-limited structural resolution. A useful future
comparison would connect sensitivity of density functionals to population
resolution while retaining their distinct targets and assumptions.

## Cryo-forum: orientation dispersion as a quality score

Read the locally archived arXiv:2307.09847v1, sections 3.6 and 5, rather than
silently treating it as the final 2024 journal text. The quaternion quadratic
representation defines Bingham dispersion statistics; eigenvalue-derived
scores rank images and permit quantile-based filtering. Its uncertainty claim
is a proxy for orientation error and particle quality. A score that improves
reconstruction after filtering does not establish a prescribed simultaneous
coverage probability for per-particle rotation balls. The proposed pipeline
also uses conventional refinement to obtain training orientations. This makes
it relevant to practical pose-quality assessment, while leaving the absolute
pose bounds in our density audit as assumptions to be justified separately.

## ARCHER: transferable reference-conditioned pose scoring

Read arXiv:2608.22029v1, including its discrete-pose classifier, continuous
refinement, independent-reference controls and supplemental quality-score
experiments. The score conditions on a supplied volume; transfer to a new
specimen is not reference-free joint reconstruction. Learned temperature and
soft angular training labels define the grid probabilities. The reported
behavioral calibration aligns synthetic operating conditions with experimental
pose-recovery performance, and negative entropy separates constructed junk
from curated particles. Neither experiment alone calibrates density confidence
intervals or a simultaneous set of small rotation balls. Reported experimental
pose differences use benchmark assignments as reference labels. This is a
useful practical comparator for pose estimation, but importing its posterior
scores as guaranteed nuisance radii would require an additional argument.

## QUTCC: image-level calibration distribution versus fixed-density coverage

Read arXiv:2507.14760v2 (24 May 2026), sections 3.1--3.3 and 5. Simultaneous
quantile regression produces spatially varying intervals; calibration adjusts
quantile inputs using paired reference images and an aggregate pixel-error
criterion. The paper explicitly limits its guarantees to pixelwise marginal
coverage and discusses rare-event and distribution-shift risks. Its imaging
experiments are denoising, MRI and quantitative phase microscopy, not cryo-EM
particle reconstruction. A conformal control here would need a stated
exchangeable population of molecules/geometries plus calibration truth; treating
voxels from one molecule as independent calibration specimens would not supply
that design. Learned conditional distributions and marginal interval coverage
should also be kept distinct. This is relevant general inverse-problem prior
art, not an unexecuted cryo-EM software baseline under another name.

## Lai et al.: empirical Bayes, dynamic reconstruction and hybrid resampling

The earlier Python TLS download failure was resolved with the Mac's verified
curl transport. Read the 18-page published main article, doi:10.5705/ss.202021.0419,
especially section 5; the supplement remains separately unreviewed. This is
direct uncertainty prior art: it connects empirical Bayes and sequential MCMC
filtering, discusses discretized random fields, and describes exact, bootstrap
and hybrid-resampling confidence sets for scalar functionals. Hybrid sets use
a fitted resampling family and approximate coverage, distinct from an exact
uniform statement over our fixed bounded density/nuisance class. The main
article does not provide a directly matched three-stack Fourier-feature coverage
benchmark. We must cite its actual confidence-set discussion rather than
dismissing it as merely a generic dynamic-imaging title. Monte Carlo approximation
error, posterior uncertainty and repeated-data coverage require separate checks.

## CryoBench: matching a metric to the reconstruction task

Read the archived primary manuscript's sections 4.1, 6 and appendix C.1. Its
synthetic heterogeneous datasets support several distinct FSC evaluations:
per-conformation reconstruction, representative samples matched to reference
states, and per-image reconstruction at an inferred latent coordinate. The
last combines state assignment and reconstruction, while maximizing over
reference states answers a different question. The paper explicitly identifies
shared-bias limitations of half-map FSC. Its discussion also identifies more
realistic noise, combined compositional/conformational variation and
nonstructural heterogeneity as extensions. Our homogeneous acquisition-geometry
simulations cannot be relabeled a CryoBench heterogeneity reproduction. Likewise,
an independent map used as a known simulation generator does not become the
unknown experimental density merely because the accession matches.

## CryoRes: learned resolution labels are not density coverage

Read the archived accepted-manuscript version of doi:10.1016/j.jmb.2023.168059,
including training-label construction and the single-map/half-map comparisons.
The network learns deposited global resolution, ResMap-derived relative local
resolution, and masks derived from atomic models; these are specified proxy
labels. Evaluation compares resolution estimates with established methods,
including half-map FSC where available. This is a map-resolution predictor,
not an estimator of a conditional density-error interval. Its use of experimental
maps does not remove uncertainty in the training targets. A fair comparison
would ask whether resolution scores predict feature recoverability under a
matched task, rather than counting them as failed confidence intervals.

## Narnhofer et al.: variance-bin conformal calibration

Read the archived manuscript dated 1 August 2024, arXiv:2212.12499, particularly
section 3.2, Proposition 3.5, Corollary 3.7 and the image-sampling caveat. The
method calibrates squared reconstruction errors within bins of estimated
posterior variance using independent signal--observation pairs. It explicitly
permits approximate posterior means and variances; sampler accuracy affects
efficiency without being the source of its rank-based coverage result. Bin-
conditional coverage is not arbitrary pointwise conditional coverage. The
experimental use of all pixels is discussed separately because same-image
pixels need not be independent. Thus a faithful cryo-EM adaptation would require
a defensible calibration sampling unit. We must not imply that posterior
variance can never be calibrated, or that its calibration automatically transfers
to every fixed molecule and imaging design.

## CESPED, cryoPARES and CryoFastAR: three pose-transfer settings

CESPED arXiv:2311.06194v2, methods and evaluation, was read in primary HTML and
archived as a PDF. Pose labels come from consistently reprocessed RELION
refinements. The authors acknowledge uncertainty in these labels and use
confidence-weighted angle errors plus reconstructed-volume comparisons. These
are useful standardized benchmarks rather than experimental ground truth.

For cryoPARES, the primary bioRxiv v6 abstract (10 August 2026) and author
repository documentation were checked; PDF retrieval returned 429. It transfers
supervised pose assignments between related samples and includes pruning and
local-refinement machinery. Full-paper validation has not been audited here.
The source limitation must accompany detailed claims about its uncertainty.

CryoFastAR arXiv:2506.05864v1, method, experimental comparisons and limitations,
was read and archived. It predicts relative Fourier-plane geometry from
multiple views, using synthetic supervision and real-data fine-tuning.
Comparisons distinguish direct inference from subsequent cryoSPARC refinement.
Its discussion reports domain-gap and image-batch limitations and excludes
heterogeneous reconstruction. These different transfer/reference assumptions
prevent interpreting all three approaches as the same ab initio or confidence-
calibration task. No original-code benchmark for these methods is claimed.


## EMMIVox: map-conditioned atomic inference

Read the 2024 PLOS Computational Biology primary article, especially the
pre-filtering, noise-prior and benchmark sections (doi:10.1371/journal.pcbi.1012180).
EMMIVox combines force fields with map likelihoods and half-map-informed
uncertainty parameters. It subsamples correlated map voxels before assuming
independence; this is not proof of exact whitening. Half-map differences guide a
lower noise-scale prior, while the authors acknowledge additional errors.
Its evaluations concern stereochemistry, map agreement and atomic ensembles.
Those outputs cannot be scored as raw-particle density confidence intervals.
The released PLUMED/GROMACS workflow would require a separately matched atomic
inference task. The molecular-dynamics examples are not executed in this project.

## Gold-standard local validation: half-map phase and amplitude

Read the Zenodo record 20730619 (v3, 17 June 2026), whose attached 47-page
manuscript filename is v5: overview, methods, threshold rationale, discussion and
supplementary SNR derivation. LocSpiral2, LocFOM, LocSNR, LocQ, LocAnisotropy and
LocResMap compare local phase/amplitude information from half maps. Their SNR
derivation explicitly assumes a shared signal, independent equal-variance noise
and unbiased reconstructions. Operational quality thresholds are not advertised
as fixed-parameter density interval coverage. Shared reconstruction bias is thus
a distinct audit question. The source is a versioned repository manuscript;
no peer-reviewed publication status was verified. Its Scipion software was
identified but has not been executed here.

## Radiative Gaussian uncertainty: adjacent CT evidence

Primary version: https://arxiv.org/abs/2607.13682v4, revised 12 September 2026.
Targeted reading: sections 1--5.4, including density moments, evaluation scope,
shared-error decomposition and ranking controls. This is a preprint, not a
verified conference publication. The current version substantially revises the
initial title and claims. Fixed-geometry density uncertainty propagates
analytically, but shared reconstruction errors can remain invisible to spread.
Its scale-only regularizer is explicitly distinguished from an ELBO with a fixed
prior. Foreground evaluation and transferred calibration are separated from
whole-volume and oracle diagnostics. This adjacent CT work limits novelty
claims about analytic Gaussian variance or the distinction between agreement
and accuracy; it does not establish raw-particle cryo-EM confidence coverage.

## cryoSENSE (CVPR 2026): compressed images versus density confidence

Primary accepted PDF, sections 2.3, 3.1 and 4, and supplementary implementation
passages were read; this is targeted reading, not a full critical review of all
experiments. [Primary source](https://openaccess.thecvf.com/content/CVPR2026/html/Shabeeb_cryoSENSE_Compressive_Sensing_Enables_High-throughput_Microscopy_with_Sparse_and_Generative_CVPR_2026_paper.html).
The method reconstructs compressed two-dimensional images using sparse or
approximate diffusion-posterior guidance. Downstream evaluation uses supplied
poses and structural/image metrics. These measurements do not establish
coverage of density-feature intervals. The stated particle-level masking versus
micrograph acquisition limitation matters when assessing experimental realism.
It motivates checking uncertainty after learned preprocessing, not treating
compression fidelity as a competing density-confidence estimator.

## CalPro: coordinate-error prediction and exchangeability unit

Primary: https://proceedings.mlr.press/v306/shihab26b.html (ICML 2026); accepted
PDF and hash in calpro-source-manifest.json. Targeted Sections 3--4 and the
coverage-semantics paragraph were read; this is not a full critical review.

The estimand is per-residue coordinate error after structural alignment. The
method applies final split-conformal calibration after fitting an evidential
head; a training surrogate does not supply its formal guarantee. Proteins are
the exchangeability unit, with averaged scores. Our reading is that the
resulting average-score event must not be silently interpreted as simultaneous
residue coverage or fixed-density coverage from particles. This is adjacent
structure-prediction work, not a direct density-interval baseline.
