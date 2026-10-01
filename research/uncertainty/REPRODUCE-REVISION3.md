# Reproducing the focused revision

Use a separate checkout for reproduction. Do not overwrite the archived results or alter frozen scientific settings. Repository records contain source snapshots, not a claim of bitwise equivalence across BLAS, FFT, compiler or platform versions.

## Environment and artifacts

Use the repository's locked dependencies and the macOS FINUFFT instructions in `COMPUTE.md`. The completed studies used one BLAS/NUFFT thread per process on the M4 Pro Mac. They require no rented GPU.

The v0.7.0-dev artifact specifications describe three refitting bundles, one comparison bundle and one raw-movie diagnostic bundle. The refitting bundles contain every calibration/test record and array, including the continuous generators, poses, observations and weights. Their archived calibration generators are sufficient to rerun the test fits; constructing the calibration generators afresh also needs the original particle metadata and earlier pilot/registration dependencies. The comparison manifest identifies its prior v0.6.* release dependencies. A completed collection can contain an interrupted historical fit; retain each record's individual completion and convergence flags.

Validate every archive member against its accompanying manifest before use. `scripts/verify_uq_release_artifacts.py` performs a streaming read without extracting or modifying the archive. All five archives and their manifests are published in [v0.7.0-dev](https://github.com/yashpatel5400/fourier-cryo-splats/releases/tag/v0.7.0-dev). All eleven remote assets, including the PDF, have byte counts and GitHub-provided SHA-256 digests matching the local files.

## Re-estimated local poses and intervals

After restoring the three archived calibration directories, run all 200 prescribed test datasets per geometry into a distinct output directory:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/run_uq_end_to_end_local.py \
  --datasets 10028,10049,10076 --output independent-refitting-reproduction
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/summarize_uq_end_to_end_local.py \
  --source independent-refitting-reproduction \
  --output independent-refitting-reproduction-summary
```

The fixed seed roots, both alignment templates, all nine procedures and both image controls are unchanged. A repeated seed reproduces a trial; it is not a new independent Monte Carlo sample. The summarizer requires all attempts and retains exceptions in the planned denominator. Compare scientific quantities and numerical residuals rather than runtime, timestamp, source-HEAD or whole-file hashes across independent executions.

To independently replay calibration itself, use `calibrate_uq_local_alignment.py --output independent-calibration-reproduction` with its recorded dependencies. The test runner deliberately reads the archived `local-alignment-calibration-v1` inputs; changing calibration is a separate study, not an interchangeable rerun.

## Continuous Gaussian comparison

With the locked targets and earlier continuous-audit arrays restored:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python scripts/benchmark_uq_continuous_gaussian_v2.py \
  --output independent-gaussian-reproduction
```

This executes all 48 fits and 672 frame/scenario comparisons. It uses the declared rank-8,192 numerical amendment. The interrupted rank-1,024 version is retained as a failed numerical attempt, not silently rerun. The report generator targets the archived study by name; use each new per-geometry JSON and arrays when comparing an independent output directory.

## Raw movie and paper

The movie specification provides its original EMPIAR URL, byte count and SHA-256. The release contains diagnostics, not a duplicate of the 1 GiB movie. The initial interrupted download, byte-range recovery and completed analysis are all retained. For a new analysis, preserve the existing `raw-movie-pilot-v1` output directory and execute the archived analysis source only after verifying the downloaded file. No picking, motion correction, poses or density inference are part of this acquisition pilot.

Build the current manuscript according to `paper/BUILD.md`. The older `write_uq_results.py` recipe reproduces a superseded table and must not be run over the current paper. A fresh build still needs rendered-page inspection; successful LaTeX compilation is not visual verification.
