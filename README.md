# Fourier Cryo Splats

**Active research; not an acceptance-ready paper.** The completely rewritten
[16-page ICML-format manuscript](output/pdf/fourier-cryo-splats.pdf),
*Auditing Population Uncertainty in Cryo-EM: Numerical Integration, Nuisance Laws,
and Identifiability*, presents the completed three-stack diagnostics and their
negative scientific decision. It does not claim a validated population-UQ method.
[Manuscript source](paper/diagnostic-main.tex) and
[proofs, survey and reproduction appendix](paper/diagnostic-appendix.tex) are public.

The [evidence and open-decisions table](research/uncertainty/CURRENT-EVIDENCE.md)
distinguishes completed outcomes from unmet research requirements:

- The final [orientation-integration repair](research/uncertainty/CATALOG-POSE-INTEGRATION-RESULTS.md)
  passes the fixed numerical gate on EMPIAR-10028, 10049 and 10076: 128/128,
  124/128 and 119/128 cases meet the between-bank tolerance. All 384 cases and
  independent checks are retained. The numerical branch is closed.
- The [population materiality screen](research/uncertainty/POPULATION-MATERIALITY-RESULTS.md)
  fails on all three stacks. Maximum primary modeled biases are .000214,
  .000122 and .003493, below .02. All thirty prelisted comparisons are retained.
  Source halves and filter groups supply viewing-law proxies, not biological
  state labels. The proposed compact-region method is not built.
- The [matched information comparison](research/uncertainty/MATCHED-INFORMATION-LEDGER-RESULTS.md)
  and [score-error report](research/uncertainty/FIXED-BANK-SEPARATION-RESULTS.md)
  distinguish descriptive discrimination from accurate numerical likelihoods.
  The [planted-truth audit](research/uncertainty/PLANTED-TRUTH-INTEGRATION-AUDIT.md)
  explains why between-bank agreement cannot certify absolute integrals.
- The [recorded-nuisance inventory](research/uncertainty/STACK-NUISANCE-INVENTORY-RESULTS.md),
  [corner spectra](research/uncertainty/BACKGROUND-SPECTRUM-INVENTORY-RESULTS.md)
  and [acquisition groups](research/uncertainty/ACQUISITION-BACKGROUND-RESULTS.md)
  document unsupported simulator assumptions. They do not calibrate pure noise.
- The [author-code population replay](research/uncertainty/POPULATION-BASELINE-REPLAY-RESULTS.md)
  completes all 23 published two-state likelihood conditions and retains
  numerical discrepancies. It does not establish experimental coverage.

All four authentic **Claude Fable 5.1 full reviews recommend rejection**, with
confidence 4/5. The [unchanged fourth review](research/uncertainty/reviews/round-04/review.md)
and [response](research/uncertainty/reviews/response-to-round-04-development.md)
remain public. The later [focused methods consultation](research/uncertainty/reviews/post-round04-method-consultation/critique.md)
and [audited response](research/uncertainty/reviews/post-round04-method-consultation/response.md)
are complete; that consultation is not a fifth full acceptance review. Useful
uncertainty, calibrated experimental nuisances and substantive novelty remain open.

The [living survey](research/uncertainty/SURVEY.md) has 141 candidates, with
source/version and reading-depth records. Retrieval counts are not full-reading
counts. [Population-UQ primary reading](research/uncertainty/POPULATION-UQ-PRIMARY-READING.md)
and [CAHRA v2 notes](research/uncertainty/CAHRA-V2-READING.md) cover direct prior
art and external controls. The downloaded CAHRA supplement checks coordinates;
the separate noisy population benchmark has not been evaluated.

The [v0.7.8 numerical evidence release](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.7.8-dev)
contains four independently verified archives and their dependency manifests.
[Reproduction instructions](research/uncertainty/REPRODUCE-POST-REVIEW4-GATES.md)
cover frozen runners, saved-array checks and the stopped branches. That release
predates this manuscript rewrite; its statement that the PDF was unchanged
was correct at its own checkpoint.

Earlier reconstruction and uncertainty studies remain in immutable releases.
The [v0.7.6 reviewed manuscript](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.7.6-dev)
is 33 pages. The [600-dataset local-refinement study](research/uncertainty/END-TO-END-LOCAL-POSE-RESULTS.md)
finds uninformative intervals, while the [matched continuous-prior comparison](research/uncertainty/CONTINUOUS-GAUSSIAN-V2-RESULTS.md)
finds no coverage advantage. [Earlier reconstruction results](RESULTS.md),
[original paper](paper/reconstruction-v0.1.0.pdf), the
[development log](research/uncertainty/DEVELOPMENT-LOG.md) and the historical
sections of [current evidence](research/uncertainty/CURRENT-EVIDENCE.md) preserve
the full record. No external GPU rental is needed for the completed diagnostics.

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

A [source-group metadata correction](research/uncertainty/BACKGROUND-SOURCE-GROUP-CORRECTION.md)
supersedes the initial missing-micrograph narrative. Existing joins recover
229 / 137 / 351 groups; all selected joins and saved split hashes replay. The
original numerical archives remain immutable and require this erratum.
