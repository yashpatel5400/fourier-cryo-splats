# Completed experimental study

The repository contains the proposed Fourier Gaussian method, an ICML-format
research preprint, background materials, executed three-accession reconstructions,
matched voxel and actual cryoDRGN classical backprojection comparisons, full FSC
curves, controls, provenance, and rerun scripts.

## What was actually run

- 8,192 experimental particles each from EMPIAR-10028, 10049 and 10076.
- 24,576 particles total, Fourier-cropped to 64×64.
- Public cryoDRGN consensus rotations, translations, and CTF parameters.
- Separate coefficient fits on 3,687-particle halves; 409 validation and 409 test
  particles per accession.
- Gaussian lattice coefficients fitted with converged preconditioned CG;
  identical-data trilinear Fourier voxel baseline.
- cryoDRGN 4.3.1 `backproject_voxel` invoked for the same halves.
- Full-frequency, windowed final runs; initial sparse-frequency runs retained.
- Regularization chosen by RAG validation prediction and frozen for other datasets.
- Phase-randomized controls for all three accessions.
- A separate 709-pair, low-frequency nonlinear MPS fit on RAG, learning centers,
  complex coefficients, and full anisotropic precisions.
- Numerical tests, source identity audit, kernel-cutoff check, and PDF visual QA.

Mean cross-method FSC against cryoDRGN backprojection is 0.9941, 0.9835, and
0.9949. Half-map FSC stays above 0.143 through sampled limits of 16.08, 7.87,
and 13.97 Å respectively. These are sampling-limit bounds, not measured finer
resolutions. The Gaussian CPU reference implementation is slower than the
matched voxel solver.

## Scientific scope

The experiment is a feasibility study on selected extracted particle images,
not raw detector-movie processing, a full-stack study, ab initio recovery,
validated adaptive high-resolution reconstruction, or near-atomic heterogeneity
analysis. Supplied consensus poses make the half-map FSC conditional, not fully
gold standard. Cross-method FSC shares data and poses and measures agreement.

Gaussian representations in cryo-EM already have substantial prior art. The
specific proposal concerns conjugate-paired kernels centered in Fourier space.
No exhaustive novelty or state-of-the-art claim is made.

## Data recovery

One EMPIAR-10028 source is truncated to 4 KB despite a header declaring 113
360×360 float32 particles. Its 32 selected particles were replaced by continuing
the same seeded selection. The public exclusion record and source-index manifests
make this reproducible. An SHA-256 audit verified all 8,007 previously recorded
particle images retained their correct source identity after reindexing.

Hardware: Apple M4 Pro, 24 GiB unified memory. The main solver used CPU; the small
nonlinear fit used MPS. macOS caffeinate kept the machine awake during the work.

Publication target: https://github.com/yashpatel5400/fourier-cryo-splats
