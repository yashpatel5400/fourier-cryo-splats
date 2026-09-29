# Uncertainty and validation in cryo-EM: living research survey

Research cutoff: 2026-09-29. This is an ongoing critical survey, not a claim
that every search hit has been read. Bibliographic candidates are in
`reading-list.tsv`; downloaded materials have hashes in `background-manifest.json`.
The broad Europe PMC search and its small retrieval discrepancies are preserved
in `europepmc-search-ledger.json`. Primary full texts, abstracts and current
preprints require different evidentiary weight. Conclusions below distinguish
the quantity a method estimates from the guarantee it actually establishes.

## 1. The uncertainty target must come first

Five quantities are often discussed together but answer different questions:

1. **Measurement noise:** how repeated exposures or independent particles would
   change an estimate under the imaging model.
2. **Reconstruction uncertainty:** which density features are consistent with
   finite data, acquisition geometry, regularization, and uncertain nuisance
   parameters.
3. **Structural heterogeneity:** the distribution of actual molecular states,
   including their populations. Physical variation is not an error bar around a
   homogeneous mean structure.
4. **Map-processing uncertainty:** how a sharpened or learned enhanced map may
   differ from a chosen reference representation.
5. **Atomic-model uncertainty:** ambiguity in coordinates, occupancy, chemistry,
   conformational ensembles, or sequence assignment given maps and other priors.

A metric can be useful for one of these without being calibrated for another.
Neither a high FSC nor a held-out image likelihood, by itself, establishes
frequentist coverage of an unknown local density. Conversely, a confidence
interval for a fixed density functional does not validate a conformational
population or an atomic interpretation.

## 2. Statistical reconstruction and explicit density uncertainty

[Scheres' Bayesian formulation](https://doi.org/10.1016/j.jmb.2011.11.010)
and [RELION](https://doi.org/10.1016/j.jsb.2012.09.006) are foundational. Their
likelihood, regularization and pose marginalization mean that Bayesian reasoning
is established cryo-EM methodology. RELION also discusses orientation accuracy
and independent-half resolution assessment. A new neural or Gaussian
representation alone is not a new statistical treatment of uncertainty.

[Ullrich et al.](https://proceedings.mlr.press/v115/ullrich20a.html) are a particularly
close comparison: a differentiable Fourier-slice model with variational density
uncertainty, including high uncertainty in unobserved directions. The paper's
discussion explicitly distinguishes variance assessment from model bias and
proposes held-out likelihood evaluation. Its diagonal covariance approximation
and its inference target must be represented fairly in a baseline. Comparing
only against a deterministic reconstruction would miss this prior art.

[Rangan et al.](https://arxiv.org/abs/2411.13263v2) study small-curvature directions
of the joint reconstruction likelihood, including coupling between density and
alignment parameters. Their examples make pose–volume ambiguity a central
issue. Projecting out nuisance Jacobians or reporting low Hessian curvature
cannot be claimed as a first recognition of that issue. The potential distinction
of the present candidate is a bounded-error confidence statement with an
explicit nonlinear remainder, rather than a local curvature diagnostic alone.

[Lai et al.](https://doi.org/10.5705/ss.202021.0419) discuss empirical Bayes and
MCMC filtering for dynamic reconstruction. Its broader dynamic-imaging framework
requires careful separation from fixed-pose single-particle density coverage;
the title alone is insufficient to establish a matched experimental baseline.

The [bootstrap variance work](https://doi.org/10.1016/j.jsb.2006.01.003) predates
modern neural methods. Resampling units matter: particles from the same source
image may share preprocessing and imaging errors. Conditional fixed-pose
resampling and rerunning the complete alignment pipeline measure different
uncertainties. A particle bootstrap does not automatically capture shared
regularization bias or a misestimated acquisition model.

## 3. Poses, nuisance parameters and identifiability

[Cryo-forum](https://doi.org/10.1016/j.jsb.2023.108058) attaches an uncertainty
measure to orientation recovery and uses it for data cleanup. The
[2026 Bayesian orientation paper](https://doi.org/10.1107/S2059798326001415)
compares posterior-averaged rotation estimation with maximum-likelihood choices.
[ARCHER](https://arxiv.org/abs/2608.22029) and
[equivariant amortized inference](https://arxiv.org/abs/2406.01630) belong in the
current pose-method comparison. Their uncertainty outputs must not be assumed
to be calibrated per-particle confidence regions without checking the actual
claim and validation.

For a local density estimator, uncertainty in rotations, translations, CTF,
image scale, and optics groups can enter differently. Independent random pose
errors can average differently from a coherent bias or a global gauge rotation.
A single combined variance cannot demonstrate robustness to both. A density
claim also needs a coordinate convention: global rotation, translation and hand
ambiguities cannot be ignored by an absolute pointwise coverage metric.

The [CAHRA challenge](https://heterogeneity.notion.site/The-Community-Wide-Assessment-of-Cryo-EM-Heterogeneous-Reconstruction-Algorithms-CAHRA-23831153834580b98df3fcba63a2f7b9)
explicitly separates compositional variation, conformational variation with
atomic modeling, and pose entanglement. Its September 2026 preprint describes
the construction of three benchmark datasets; the organizers say a later update
will analyze full competition results. This is direct community evidence that
distinguishing structural variation from imaging ambiguity remains important.

## 4. Heterogeneity: structures, coordinates and populations

[CryoDRGN](https://doi.org/10.1038/s41592-020-01049-4) learns a continuous neural
density family. Its latent encoder distribution should not be casually equated
with uncertainty in a reconstructed density or with a calibrated population
distribution. [CryoDRGN-AI](https://doi.org/10.1038/s41592-025-02720-4) addresses
ab initio reconstruction; classical `cryodrgn backproject_voxel` is not a
substitute for that neural baseline.

[RECOVAR](https://doi.org/10.1073/pnas.2419140122) uses regularized covariance
estimation and kernel regression, accounting for noisy latent coordinates when
estimating populations. [3DFlex](https://doi.org/10.1038/s41592-023-01853-8) and
[DynaMight](https://doi.org/10.1038/s41592-024-02377-5) model structural motion.
DynaMight's independent-half validation and held-out deformation comparisons
are especially relevant: shared learned components can contaminate validation,
and disagreement between deformation estimates targets a different quantity
from density confidence coverage.

Gaussian and structural representations are already substantial prior art:
[mixed-dimensional GMMs](https://doi.org/10.1038/s41592-021-01220-5),
[cryoSPHERE](https://arxiv.org/abs/2407.01574),
[CryoSTAR](https://doi.org/10.1038/s41592-024-02486-1),
[CryoSPIRE](https://arxiv.org/abs/2506.09063),
[CryoSplat](https://arxiv.org/abs/2508.04929v5), and
[GEM](https://arxiv.org/abs/2509.25075v2) preclude claiming Gaussian cryo-EM
reconstruction as the new contribution. The last two use real-space Gaussian
representations; the current implementation uses conjugate-paired Fourier
kernels. That representation difference needs demonstrated statistical or
computational value.

[Evans et al.](https://doi.org/10.1038/s42003-026-09859-6) show that hard or soft
particle assignments can give distorted populations under noise. Whole-stack
likelihood/reweighting and deconvolution address the population inverse problem.
Their real mixed-particle experiment uses a constructed baseline population,
which the authors distinguish from unknowable biological truth.
[Mattingly et al.](https://arxiv.org/abs/2606.14449) analyze measurement-limited
heterogeneity and information-based coarse graining. Identifiability and the
appropriate granularity of a structural ensemble are already active theory
topics, not unexplored consequences of a new latent model.

## 5. FSC, local significance and validation beyond agreement

Independent halves and [overfitting prevention](https://doi.org/10.1038/nmeth.2115),
[high-resolution noise substitution](https://doi.org/10.1016/j.ultramic.2013.06.004),
and [modified FSC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7642792/) address
different failure modes. An FSC confidence interval concerns the correlation
statistic and its effective sampling assumptions; it is not a pointwise density
confidence interval. Common masks, fitted poses, reference maps and regularizers
can change the interpretation of half-map agreement.

[ResMap](https://doi.org/10.1038/nmeth.2727) tests local sinusoidal structure;
[MonoDir](https://pmc.ncbi.nlm.nih.gov/articles/PMC6940361/) quantifies directional
resolution. [FDR confidence maps](https://doi.org/10.1107/S2052252518014434) and
[subsequent confidence-map inference](https://pmc.ncbi.nlm.nih.gov/articles/PMC7137106/)
provide significance-based map interpretation. False-discovery control, familywise
coverage, and a posterior probability are different statements. Resolved versus
unresolved signal is also not identical to a bound on reconstruction error.

[Independent-particle validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC7385033/)
is a useful complement to map comparison. The
[bias/variance analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8972802/),
[Lander's validation perspective](https://doi.org/10.1016/j.sbi.2024.102918), and
[Fourier-space deposition proposal](https://pmc.ncbi.nlm.nih.gov/articles/PMC12676377/)
motivate preserving uncertainty-relevant intermediate data, rather than reporting
only a sharpened map and a scalar resolution.

## 6. Learned priors and uncertainty in enhanced maps

[Blush](https://doi.org/10.1038/s41592-024-02304-8) brings learned regularization
into reconstruction; [diffusion priors](https://arxiv.org/abs/2412.14897) expand
the Bayesian imaging toolkit. Any learned prior can improve average error while
making a particular unsupported feature look stable. The relevant controls are
held-out measurements, prior mismatch and deliberate unsupported-feature tests,
not just visual sharpness or agreement between identically regularized halves.

[LocScale-2.0](https://doi.org/10.1038/s41467-026-75327-8) uses Monte Carlo dropout
and calibrates uncertainties against a phase-conservative baseline map. That
baseline is a proxy, not the unknown true density. Its scores address differences
between the enhanced and baseline maps and support interpretation of map
processing. A matched evaluation must preserve this target instead of scoring
its output as though it promised conditional coverage of particle-derived
linear density functionals.

[CryoDiff](https://doi.org/10.64898/2026.06.04.730282) is a recent uncertainty-aware
map-enhancement preprint; full-text access and exact calibration details remain
to be resolved before a detailed comparison. [EMReady2](https://doi.org/10.1038/s41467-026-71794-1)
uses local-quality-aware learning and evaluates map interpretability. Its quality
scores and reconstruction improvements do not automatically define calibrated
uncertainty intervals. Model-derived evaluation maps can share structural priors
with an enhancement network, a dependency that a validation design must inspect.

## 7. Atomic ensembles and tomography

[Habeck's Bayesian modeling](https://doi.org/10.3389/fmolb.2017.00015),
[EMMIVox](https://doi.org/10.1371/journal.pcbi.1012180), and
[cryoENsemble](https://doi.org/10.1038/s41598-024-68468-7) concern molecular models
or ensembles constrained by maps. These are scientifically important but are
not interchangeable with confidence for a density reconstructed from raw
particle images. Map-to-model tools such as FSC-Q and DAQ also assess distinct
parts of the inference chain.

[tomoDRGN](https://doi.org/10.1038/s41592-024-02210-z) and
[CryoDRGN-ET](https://doi.org/10.1038/s41592-024-02340-4) extend learned structural
variation to tomography. Missing wedges, shared tilt-series noise and alignment
dependencies make particle independence particularly consequential there.
Tomographic examples would be an extension of the current single-particle model,
not evidence that a single-particle proof already covers in situ reconstruction.

## 8. Benchmarks and the contribution threshold

[CryoBench](https://cryobench.cs.princeton.edu/) supplies diverse known-truth
heterogeneity tasks. [Roodmus](https://doi.org/10.1107/S2052252524009321) provides
a simulation/benchmarking toolkit. The
[Flatiron challenge](https://pmc.ncbi.nlm.nih.gov/articles/PMC12637659/) evaluates
heterogeneous reconstructions using multiple complementary views of performance.
These and CAHRA argue for reporting structure, population, measurement fit,
calibration and compute separately. One scalar FSC cannot answer all of them.

The general inference construction here is grounded in
[Donoho](https://doi.org/10.1214/aos/1176325367) and
[Armstrong–Kolesar](https://doi.org/10.3982/ECTA14434).
[Recent simultaneous inverse-problem calibration](https://arxiv.org/abs/2510.11708)
is also directly relevant. Conformal imaging methods, including
[posterior-variance calibration](https://arxiv.org/abs/2212.12499), must be assessed
with their exchangeability assumptions and calibration targets intact. Calibration
against simulated or paired ground-truth maps cannot be silently transferred
to an arbitrary new molecule with uncertain poses.

The strongest current candidate is **structural-feature uncertainty that exposes
pose ambiguity and reconstruction bias under an explicit sensitivity class**.
To justify a paper, it must offer more than the classical support-function formula:
useful nonlinear error bounds, scalable computation, defensible choice/reporting
of sensitivity bounds, and compelling comparisons on known-truth and real-data
tasks. The pilot's initially vacuous bounds demonstrate why this remains work
to be done.

Other open questions worth tracking are calibrated population inference under
pose misspecification; validation of learned structural detail absent from the
measurements; cross-specimen uncertainty transfer; data acquisition chosen for
specific structural questions; and practical uncertainty propagation through
movie correction, CTF estimation, alignment, classification and model building.
These are research directions inferred from the sources above, not claims of
unclaimed priority.

## 9. Additional targeted search, 29 September 2026

Four additional queries combined cryo-EM with conformal inference, posterior
calibration, Gaussian/ bootstrap uncertainty, and feature validation. They did
not identify a directly matched conformal density-reconstruction paper in the
returned results; that is a search outcome, not a proof of absence. Three further
candidates brought the then-current ledger to 89 entries (historical checkpoint;
the current ledger has 93).

[SIMPLE's 2025 probabilistic ab initio method](https://doi.org/10.1107/S2059798325005686)
uses coupled orientation assignments and adaptive spatial regularization. Its
sampling neighborhoods guide optimization, and should not be treated as
calibrated confidence sets for individual poses. The paper also ties acquisition
and reconstruction speed to useful online processing, an important practical
comparison beyond a final FSC.

[CryoETGS](https://doi.org/10.1016/j.jsb.2025.108281) adds published cryo-ET Gaussian
representation work to the existing single-particle Gaussian literature. The
primary abstract, rather than a verified full-text reproduction, supports its
inclusion here. [Experiment-guided AlphaFold3](https://doi.org/10.1038/s41587-026-03166-5)
belongs to the atomic-ensemble branch: agreement of ensemble-averaged observables
with experiments addresses a different target from particle-derived density CIs.
Its full text has been checked for inputs, cryo-EM guidance and independent-modality limitations; see the critical annotations.


## 10. Target-specific reading updates

Additional primary-source annotations now cover DynaMight, 3DFlex, RECOVAR,
particle-counting bias and measurement-limited ensemble selection. These refine
the distinction between per-image latent-coordinate error, variability of
physical structures, population uncertainty, shared-reference bias and the
homogeneous density-functional target used here. RECOVAR's final PMC record and
its preprint record are distinguished; direct full-page access limitations are
recorded rather than silently treating all versions as interchangeable.

The present implementation has progressed from finite dictionary/voxel audits
to continuous L2 density and bounded nonlinear-pose post-audits. Completed
assumption controls show failures from support exclusion, heterogeneous images
and misspecified noise, while some tested defocus/heavy-tail controls remain
conservative. A frozen additional-exposure prediction study now favors the
neural baseline slightly on all three stacks. Those results do not change the
survey's novelty threshold or supply experimental density-coverage labels.
Both frozen continuous-uncertainty studies are now complete. Their analytic
coverage calculations verify conditional implementations; poor sign power and
conservative widths limit practical inference.

The pose-transfer follow-up adds [CESPED](https://arxiv.org/abs/2311.06194v2),
[cryoPARES](https://www.biorxiv.org/content/10.1101/2025.03.04.641536v6), and
[CryoFastAR](https://arxiv.org/abs/2506.05864v1), bringing the curated ledger to
92 entries at that historical checkpoint (93 currently). CESPED standardizes refinement-derived pose labels and explicitly
acknowledges their uncertainty. CryoPARES reuses alignments of related specimens;
CryoFastAR learns multiview pose prediction from synthetic training with real
fine-tuning. These are distinct transfer assumptions. Their quality scores,
reconstruction comparisons, and inference speed do not by themselves calibrate
the bounded nuisance sets used in this project.

The Lai et al. main article is now available and critically annotated. It
explicitly discusses scalar confidence sets and hybrid resampling in a multilevel
empirical-Bayes treatment, so its relationship to our construction is more
substantive than a shared title. The broader survey still distinguishes targeted
primary reading from exhaustive screening or original-code reproduction.
