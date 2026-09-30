# Fourier Cryo Splats

**Active research; not an acceptance-ready paper.** The
[current ICML-format manuscript](output/pdf/fourier-cryo-splats.pdf) studies
uncertainty in cryo-EM density features through continuous Fourier-slice bias
and pose auditing. The original Gaussian reconstruction work remains available
in the [earlier paper](paper/reconstruction-v0.1.0.pdf) and v0.1.0 release.

The first authentic **Claude Fable 5.1 full review recommends rejection**, with
confidence 4/5, and does not consider the work a strong ICML contender. The
[unaltered review](research/uncertainty/reviews/round-01/review.md),
[current response](research/uncertainty/reviews/response-to-round-01-development.md)
and [focused mathematical audits](research/uncertainty/reviews/README.md) are
public. Focused audits do not replace a full acceptance assessment. Experimental
pose/noise/class calibration, useful fine-scale inference and substantive novelty
remain open; neither passing numerical checks nor approximate-map inclusion
establishes experimental density coverage.

The compact [evidence and open-decisions table](research/uncertainty/CURRENT-EVIDENCE.md)
distinguishes completed outcomes from the remaining running fits.

Current evidence includes:

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
  and left the full-space optimization gap unresolved.
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
  reconstruction on 10028 and a weaker converged result on 10049; the 10076
  continuation remains in progress. The second result has registered-reference
  mean FSC 0.202 despite its convergence flag.
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
RELION result was already described in that release's appendix.

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
