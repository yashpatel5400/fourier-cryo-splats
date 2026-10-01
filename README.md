# Fourier Cryo Splats

**Active research; not an acceptance-ready paper.** The
[current ICML-format manuscript](output/pdf/fourier-cryo-splats.pdf) studies
uncertainty in cryo-EM density features through continuous Fourier-slice bias
and pose auditing. The original Gaussian reconstruction work remains available
in the [earlier paper](paper/reconstruction-v0.1.0.pdf) and v0.1.0 release.

The fourth authentic **Claude Fable 5.1 full review recommends rejection**, with
confidence 4/5, and does not consider the work a strong ICML contender. The
[unaltered round-4 review](research/uncertainty/reviews/round-04/review.md),
[current response plan](research/uncertainty/reviews/response-to-round-04-development.md)
and [focused mathematical audits](research/uncertainty/reviews/README.md) are
public. Focused audits do not replace a full acceptance assessment. Experimental
pose/noise/class calibration, useful fine-scale inference and substantive novelty
remain open; neither passing numerical checks nor approximate-map inclusion
establishes experimental density coverage.

The compact [evidence and open-decisions table](research/uncertainty/CURRENT-EVIDENCE.md)
distinguishes completed outcomes from unresolved scientific limitations.

Round 4 calls for a substantive change, not more variants of the current moment
test. Those variants are frozen while a matched information/nuisance ledger and
the three stacks' recorded imaging conditions are investigated. The existing
33-page PDF is the reviewed v0.7.6-dev snapshot, not a revision addressing this
review. All four full reviews reject. The living survey now has 139 candidates,
with [new population-UQ primary reading](research/uncertainty/POPULATION-UQ-PRIMARY-READING.md)
and [CAHRA version-2 reading](research/uncertainty/CAHRA-V2-READING.md); candidate
counts are not full-reading counts.

Post-review-2 work includes the [complete matched continuous Gaussian comparison](research/uncertainty/CONTINUOUS-GAUSSIAN-V2-RESULTS.md): all 48 fits converge and all 672 conditional cases are retained. The [600-dataset local-refinement study](research/uncertainty/END-TO-END-LOCAL-POSE-RESULTS.md) also finishes, following a separate 384-dataset calibration batch. All 12,000 weight solves converge. Every one of 216 procedure/control cells covers its coarse truth in 200/200 trials, including the simpler baselines; this is **not a demonstrated coverage advantage**. The three pose-audit procedures always select the no-data interval. Their raw median widths span 642–1,242 times its width, with nonlinear remainder dominating numerically; even removing it leaves all first-order bounds uninformative. Near-truth local pose initialization deteriorates substantially during the declared coarse fit; these are not production reconstruction results.

The [rewritten manuscript](paper/focused-main.tex) reports these failures in its abstract and main results. The [registered replay](research/uncertainty/REGISTERED-DICTIONARY-RESULTS.md), [breakdown analysis](research/uncertainty/BREAKDOWN-RADIUS-RESULTS.md), and [phase independence control](research/uncertainty/PHASE-SPLIT-CONTROL-RESULTS.md) are complete. A [raw-movie pilot](research/uncertainty/RAW-MOVIE-PILOT-RESULTS.md) retains all acquisition diagnostics and the interrupted/recovered transfer, but does not establish independent pure noise. None resolves the experimental inputs. [Current reproduction instructions](research/uncertainty/REPRODUCE-REVISION3.md) distinguish immutable archived studies from new runs.

The [post-review-3 saved-array diagnosis](research/uncertainty/REFITTING-BIAS-REANALYSIS-RESULTS.md) reproduces every original estimate and adds 18 converged true-pose fits. Noise-only coverage falls to .720 for the pilot-aligned 10049 contrast, versus .955 with true poses; the original broad bias-aware intervals still cover every trial. Exact realized-pose class envelopes are much smaller than the old bounds but remain large. These post hoc results are retained in the current paper and do not resolve the rejected novelty or usefulness assessment.

The new [translation-invariant moment diagnostic](research/uncertainty/paired-power-v1/BISPECTRUM-RESULTS.md) retains all 24 cells and 160 continuous-orientation searches. Two ranged-amplitude contrasts fail on 10028; the four 10049 contrasts retain sampled positive margins without a global validity certificate; no positive direction was found on 10076. All twelve removal fits reach their iteration limit. An independent complex-moment calculation checks 48 noise-variance values within 3.6e-15. This is oracle method development, not new experimental coverage or a new reconstruction. The updated [literature reading](research/uncertainty/ALIGNMENT-LITERATURE-REVIEW3.md) distinguishes collective recovery, individual alignment and prior moment-based posterior methods; the ledger now has 138 candidates, with reading depth stated separately.

A subsequent [global-cover cost calculation](research/uncertainty/paired-power-v1/BISPECTRUM-GLOBAL-BOUND-RESULTS.md) rules out that particular derivative-envelope grid implementation, requiring 1e13–2e16 rotations. The separate [Monte Carlo study](research/uncertainty/paired-power-v1/MONTE-CARLO-VIEW-LAW-RESULTS.md) changes the model to a bounded viewing density and known noise. At 10,000 Haar-view particles, ranged-amplitude 10049 rejection projections are .700 for power and .826 for combined moments, with overlapping intervals. Both ranged 10028 contrasts fail, 10076 remains uninformative, and the guarantee becomes ineffective for 10049 at a viewing-density ratio bound of 1.1. These are binomial projections from independent simulations, not observed experimental power or a non-Haar experiment. All 192 projections, controls, scores and numerical checks are retained; [reproduction instructions](research/uncertainty/REPRODUCE-REVISION6.md) state the dependencies and assumptions.

The [conditional-noise follow-up](research/uncertainty/paired-power-v1/VIEW-VARIANCE-RESULTS.md) and [classical risk comparisons](research/uncertainty/paired-power-v1/VIEW-RISK-BASELINES-RESULTS.md) separate conditional noise from viewing variability. The added CVaR comparators are looser at this budget, but the paired variance bounds are dominated by Monte Carlo concentration slack. Grouped methods require amplitude conditionally independent of noise given view. A [fresh preferred-view study](research/uncertainty/paired-power-v1/PREFERRED-VIEW-RESULTS.md) uses 1,179,648 rotations, retains 1,944 projections and 648 actual repeated-group cells, and finds strong amplitude dependence. This remains a selected two-stack study; 10076 is not newly resolved. The [authentic focused Fable audit and response](research/uncertainty/reviews/view-variance-audit-01/response.md) identify remaining calibration, novelty and usefulness gaps. The v0.7.4-dev paper is 30 pages; [reproduction instructions](research/uncertainty/REPRODUCE-REVISION7.md) cover all new numerical arrays and checks. All three full reviews remain rejections.

The [candidate-derived score study](research/uncertainty/paired-power-v1/CANDIDATE-SCORE-RESULTS.md) now tests all three fitted Gaussian maps without an external reference for region or score design. All 18,900 scalar projections and 6,300 repeated-group cells are retained. Removing scale means helps, but small deletions remain difficult, viewing uncertainty strongly reduces power, and classical CVaR outperforms the paired bound on some candidates. The largest correct-null repeated-group count is 7/128, with a wide pointwise interval; calibration is fixed. The v0.7.5-dev paper is 32 pages. [Reproduction instructions](research/uncertainty/REPRODUCE-REVISION8.md) include the actual Gaussian candidate maps, every numerical array, independent checks and the simulator limitations. A candidate rejection does not uniquely localize structural error or supply regional-density coverage.

A [fresh classical Fisher comparison](research/uncertainty/paired-power-v1/FISHER-SCORE-RESULTS.md) adds 18,900 projections, 6,300 repeated-group cells and 9,450 conservative difference intervals on all three stacks. Covariance-aware scores improve the largest 10028 case but leave smaller changes difficult. An [equal-noise-budget allocation study](research/uncertainty/paired-power-v1/REPLICA-ALLOCATION-RESULTS.md) adds 37,800 projections and 12,600 group cells on another independent test sample: more replicas reduce estimated nuisance envelopes, but fewer views often loosen confidence bounds. Neither allocation uniformly improves power. All outcomes and full training inputs are retained; 25 targeted tests pass and separate implementations check covariances, constrained directions, bounds and outcome records. The current paper is 33 pages; [reproduction instructions](research/uncertainty/REPRODUCE-REVISION9.md) describe the four numerical bundles. These are conditional simulator diagnostics, and all three full Fable reviews remain rejections.

Earlier evidence includes:

- [Direct folded-width optimization](research/uncertainty/FOLDED-RIDGE-REVIEW3-RESULTS.md)
  completes six targets and 119 converged linear solves, reducing widths only
  .062--.098%. The updated paper retains restricted and global gaps separately.
- [Paired-exposure validation gates](research/uncertainty/paired-power-v1/METHOD-GATES-RESULTS.md)
  retain diagonal-cone failures, all full-frequency/common-shift cases, and
  42 verified continuous-pose violations. The exact-model focused audit and
  response are public. These candidates are not calibrated experimental tests.

- Original-code CryoLike scoring completes all 24 declared three-stack cases.
  Both metrics and both viewing grids are retained; their rankings differ.
  These reused-particle scores do not establish density calibration.

- Two completed frozen conditional-uncertainty studies and an additional-exposure
  prediction comparison on three EMPIAR stacks. Stock neural prediction wins
  all six prescribed contrasts; this is prediction, not density uncertainty.
- A 10,000-particle, radius-12 matrix-free audit using 1.95 GB peak memory on
  the Mac. The completed eighteen-case pose-aware design grid retains failed
  convergence and unfavorable power results; its two-degree bounds still miss
  the full review's usefulness and tightness criteria.
- All twelve locked pilot-selected fixed-pose feature fits are complete. Their
  known-noise reference sign power depends on correct map registration. Their nonlinear
  audits and 48 experimental intervals are complete. All twelve
  estimators and 116 files were frozen and published before downloading the
  new calibration cohort: 128 particles from unused exposures per stack.
  All three downloads and the fresh recalibration are complete. Six fixed-pose
  intervals and four shift-only intervals exclude zero; none with rotational
  uncertainty do. A later, separately declared centered-noise procedure raises
  those counts to ten and six, with rotational exclusions still zero.
  [All centered outcomes](research/uncertainty/CENTERED-NOISE-CALIBRATION-RESULTS.md)
  retain the original procedure and five new disagreements with the Class A
  reference on heterogeneous 10076. Approximate-reference inclusion is not density coverage.
- A completed cubic audit reduces one selected one-degree width from 0.473
  to 0.252 of no-data width, without useful sign power. Optimized cubic weights
  reduce this to 0.180, with a large remaining optimization gap. A separately
  declared reference-frame check changes its simulated minimum sign power from
  0.00654 to 0.95725 without changing the estimator. Its subsequent experimental
  interval still contains zero: the experimental SD upper bound is much larger
  than the supplied simulation SD and is not a pure-noise measurement.
  [All three fixed-estimator outcomes](research/uncertainty/REGISTERED-TARGET-SENSITIVITY-RESULTS.md)
  are retained. The coordinate follow-up completed with a
  worse width of 0.223. The reduced convex design selected the original weights
  and left the full-space optimization gap unresolved. The completed enriched
  and joint designs give relative widths .178365 and .175192, modest .79% and
  2.55% improvements. All ten raw/centered experimental intervals for the five
  designs contain zero; the joint centered interval is 2.46896 ± 3.01088.
  [All new outcomes](research/uncertainty/JOINT-ENRICHED-EXPERIMENTAL-RESULTS.md)
  remain conditional on unverified pose, density and noise assumptions.
- Matched bootstrap and [Fourier Gaussian baselines](research/uncertainty/FOURIER-VARIATIONAL-BASELINE.md),
  including favorable broader-prior results. A completed
  [local Gaussian pose comparison](research/uncertainty/FOURIER-POSE-BASELINE-RESULTS.md)
  widens its intervals only slightly; this pilot linearization omits nonlinear
  and density/pose interaction terms. Positivity/support controls,
  two-pose ambiguity witnesses, higher-band and nuisance-sensitivity failures
  remain reported; they do not rescue the current practical limitations.
- A [known-map information diagnostic](research/uncertainty/ORACLE-POSE-INFORMATION-RESULTS.md)
  distinguishes individual pose information from collective reconstruction.
  The native RELION CPU comparison has completed a converged unknown-pose
  reconstruction on all three stacks. Registered-reference mean FSC is
  .702/.202/.122 on 10028/10049/10076, respectively, despite convergence flags
  on all three; every original and registered curve is retained.
  Pilot-only reference registration improves cross-method reference FSC on
  10028 and 10049; the Class A reference on heterogeneous 10076 remains a
  weak match to the consensus. Both original and registered curves are kept.
- A separate [likelihood-validation candidate](research/uncertainty/MIXTURE-VALIDATION-CANDIDATE.md)
  targets structural compatibility. Its discrete-view, oracle-predictor screen
  and shared-scale follow-up detect full local removal on two of three simulated
  geometries, with all weaker-change failures retained. The first full-orientation
  denominator calculation and two diagnostics completed with unresolved large
  gaps; the completed refinement still leaves gaps of 12,898--15,892. A reviewer-
  motivated disk diagnostic does not improve sampled coarse cells. These are development feasibility tests;
  a practical independent predictor and experimental calibration remain missing.

The [research plan](research/uncertainty/PLAN.md),
[critical survey](research/uncertainty/SURVEY.md),
[theory](research/uncertainty/THEORY.md),
[development log](research/uncertainty/DEVELOPMENT-LOG.md) and
[reproduction commands](research/uncertainty/REPRODUCE-DEVELOPMENT.md) give the
full methods, sources, protocols and outcomes. The survey distinguishes curated
candidates, targeted reading and retrieval; a download is not a full reading.

Exact early frozen models are in the
[v0.2 uncertainty checkpoint](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.2.0-dev).
The immutable [v0.6 development checkpoint](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.6.0-dev)
contains 905 verified arrays and its historical 50-page preprint. The
[v0.6.1 increment](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.6.1-dev)
adds seven verified arrays and the corrected 52-page manuscript with registered
and centered-calibration outcomes. The historical v0.6 main-text statement that
no converged ab initio comparison had run is corrected: 10028's completed
RELION result was already described in that release's appendix. The
[v0.6.2 completed increment](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.6.2-dev)
adds 128 verified numerical/map/diagnostic files and the 55-page manuscript,
including all three completed RELION comparisons and all new cubic outcomes.
It requires the prior bundles described in its manifest.
The [v0.7.0 focused revision](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.7.0-dev) adds five verified numerical bundles (1.19 GB total) and the 26-page rewritten paper. It retains all 600 refitting trials, 384 calibration datasets, favorable matched-prior comparisons and the raw-movie diagnostics; the proposed pose bounds fail the informativeness criterion.

On macOS, installing the optional uncertainty tools requires the FINUFFT build
configuration in [COMPUTE.md](research/uncertainty/COMPUTE.md) to avoid conflicting
OpenMP runtimes. Numerical verification: `OPENBLAS_NUM_THREADS=4 pytest -q`.

A research implementation of **conjugate-paired Gaussian kernels centered in
Fourier space** for single-particle cryo-EM reconstruction. The pair structure
ensures a real inverse transform. An anisotropic Gaussian has an exact analytic
restriction to each particle's central Fourier plane; CTF and translation are
applied in that plane.

This repository includes an ICML-format manuscript, actual experimental-data
reconstruction code, matched voxel and cryoDRGN classical backprojection
comparisons, numerical tests, provenance, and complete FSC curves.

[Earlier reconstruction results](RESULTS.md) ·
[Reconstruction artifacts](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.1.0)

**Original reconstruction scope:** a three-accession feasibility study using 8,192 extracted experimental
particles per accession, Fourier-cropped to 64×64, with supplied consensus poses.
It is not ab initio, a detector-movie processing pipeline, a full-stack benchmark,
or evidence of near-atomic heterogeneous reconstruction. Half-map FSC is
conditional on the supplied poses. Cross-method FSC measures agreement, not
independent accuracy. Gaussian cryo-EM methods already exist; see
[the related-work review](background/README.md).

## Method

For real Gaussian kernels `g`, represent a Fourier volume by

```
F(k) = sum_j [ c_j g(k; mu_j, Lambda_j)
            + conj(c_j) g(k; -mu_j, Lambda_j) ].
```

For a central plane `k = E q`, the restricted precision is `Eᵀ Lambda E`.
Its center and attenuation follow by completing the square. This is a
restriction of a spectral Gaussian, not the covariance projection used to
integrate a spatial Gaussian.

The main benchmark fits complex coefficients in a stationary Gaussian lattice
by preconditioned conjugate gradients, with a coefficient ridge penalty. A
separate differentiable implementation supports learned centers and full
anisotropic precision; its small real-data demonstration is reported separately.
No claim of a validated adaptive high-resolution splatting engine is made.

## Reproduce

The completed run used Python 3.13 on an Apple M4 Pro with 24 GiB unified memory.
The main sparse solves use CPU; the nonlinear demonstration uses PyTorch MPS.
Linux/CUDA is not required for the main reconstruction. Allow approximately
12 GB of network transfer and several GB of local storage; the exact downloaded
byte counts are in `provenance/`. Original selected particle bytes are downloaded,
hashed, Fourier-cropped, and discarded; the downsampled images remain locally.

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
pip install -e .
git clone https://github.com/zhonge/cryodrgn_empiar background/cryodrgn_empiar
git -C background/cryodrgn_empiar checkout 886c5d106f98a7d9647c1cc60e48a4cefa20ffaa
python scripts/gather_background.py
pytest -q
python scripts/validate_synthetic.py
```

For the complete measured experimental sequence after installation and metadata setup,
run `bash scripts/reproduce.sh`. The individual commands below explain its stages.

Download the three selections. The downloader verifies HTTP Content-Range and
hashes every response; exact source indices and byte ranges are saved. A
documented truncated EMPIAR-10028 source file is excluded automatically; the
same seeded selection fills its 32 affected positions with usable particles.
Interrupted downloads can be resumed with `--resume`.

```sh
python scripts/download_data.py 10028
python scripts/download_data.py 10049
python scripts/download_data.py 10076
python scripts/package_provenance.py
```

The final settings use all 1,410 nonredundant Fourier pixels inside radius 30,
an image window tapering from 0.85 to 0.99 of the half-box, Gaussian width 0.5
Fourier bins, truncation radius 2 bins, and regularization factor 1. The regularizer
was selected using only the RAG validation set from candidates 0.1, 1, and 10,
and then frozen across accessions and both representations.

```sh
python scripts/run_experiment.py 10028 --samples 2048 --window --reg 1 --tag full-reg1
python scripts/run_experiment.py 10049 --samples 2048 --window --reg 1 --tag full-reg1
python scripts/run_experiment.py 10076 --samples 2048 --window --reg 1 --tag full-reg1
```

`--samples 2048` is capped to all 1,410 available independent frequencies. Both
halves are fitted separately; 409 particles are reserved for validation and 409
for test, leaving 3,687 particles in each training half. Global amplitude scaling
uses training images only. Test images never enter the normal equations.

Run the external classical backprojection for each accession with the same
particle partitions; this invokes installed cryoDRGN rather than reimplementing
or relabeling it. It is **not** cryoDRGN's neural or ab initio command.

```sh
python scripts/run_cryodrgn_baseline.py 10028 --tag full-reg1
python scripts/run_cryodrgn_baseline.py 10049 --tag full-reg1
python scripts/run_cryodrgn_baseline.py 10076 --tag full-reg1
```

The original experimental sequence retained an initial 256-Fourier-sample,
unwindowed run under `results/main/`. The external runs used its identical
particle partitions, and final analyses reused those external maps. See saved
commands and timing records for exact provenance. This order affects filenames,
not which particles were reconstructed.

Controls and nonlinear demonstration:

```sh
python scripts/run_experiment.py 10049 --box 32 --samples 256 --window --reg 1 --phase-randomize --tag noise
python scripts/fit_adaptive_demo.py --dataset 10049 --steps 1500
python scripts/plot_results.py --tag final
```

The phase-randomized control is also run on the other two accessions. The
synthetic phantom is generated from displaced spatial Gaussians, a different
basis from the reconstruction dictionary. Unit tests cover analytic restriction,
Hermitian symmetry, adjoints, CTF symmetry, FSC, inverse-transform reality, and
nonlinear geometry gradients.

## Outputs and interpretation

Binary maps and coefficient checkpoints are supplied in the GitHub release
archive. Extract it into the repository root to populate the paths below.
The smaller JSON metrics, FSC CSVs, figures, and provenance records are in Git.

- `paper/`: official ICML 2026 template, LaTeX source, bibliography, and initial draft.
- `output/pdf/fourier-cryo-splats.pdf`: current uncertainty development manuscript.
- `results/final/`: reported metrics, unmasked FSC curves, predictions and maps.
- `results/main/`: initial sparse-observation ablation, retained transparently.
- `results/full-reg*/10049/`: regularization-selection runs.
- `results/noise/`: phase-randomized controls.
- `results/adaptive-demo/`: a separate low-frequency nonlinear demonstration.
- `provenance/`: download manifests, SHA-256 hashes, and original particle indices.
- `background/`: source catalog and research-positioning notes.

FSC crossings use a sustained two-shell crossing. If no crossing occurs, the
reported value is censored at the sampled limit: it must not be presented as a
measured finer resolution. The image window is a conventional preprocessing
approximation. No real-space mask, sharpening, or test-based alignment is applied
to FSC maps. Plotting-only smoothing of slice figures is labeled separately.

The current Gaussian solver is slower than the matched trilinear solver. Sparse
kernel geometry learning and a compiled renderer are possible future extensions;
this implementation does not establish a speed advantage. Common data, poses,
preprocessing, and smoothness can increase cross-method agreement, which must
not be mistaken for an independent resolution estimate.

## License and publication

Project code is MIT licensed. Third-party templates retain their notices.
Original publications, cryoDRGN inputs, and EMPIAR particles remain governed by
their respective sources and are not redistributed as project-owned material.
This is a research preprint in ICML format, not an ICML submission or acceptance.
