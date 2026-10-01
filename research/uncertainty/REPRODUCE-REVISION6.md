# Reproducing invariant-moment and bounded-view development

The v0.7.3-dev increment preserves every completed outcome of the moment-hull,
adversarial-pose, global-cover cost and bounded-view Monte Carlo gates. It is
not a favorable acceptance review or experimental calibration. The 28-page
manuscript includes the new theory, outcomes and limitations. All three full
Fable reviews remain rejections.

Use the tagged source and pinned Python environment. Download the
`v0.7.3-dev-invariant-moments.tar.gz` asset and adjacent manifest, plus the
v0.7.2-dev method-gate archive and its manifest. Verify before extracting at
the repository root:

```bash
.venv/bin/python scripts/verify_uq_release_artifacts.py output/artifacts/v0.7.2-dev-method-gates.tar.gz output/artifacts/v0.7.3-dev-invariant-moments.tar.gz
```

The earlier bundle supplies the exact saved 4,160-view/10,000-view arrays,
three EMDB maps and preflight geometry. The new increment supplies all eight
new NPZ files, all summaries/CSVs, source, protocols, verification records
and logs. Third-party literature PDFs are not redistributed: their URLs,
reading scope, failed retrievals and downloaded-file hashes are retained.
Older raw-particle reconstruction reproduction uses the preceding releases.

Saved-array verification requires no new EMPIAR download:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_bispectrum.py --output /tmp/bispectrum-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_moment_mc.py --output /tmp/moment-mc-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_bispectrum.py tests/test_moment_rotation_bound.py tests/test_moment_mc.py
bash scripts/build_paper.sh
```

Choose unused verification output paths. The first verifier checks all twelve
saved hull mixtures and independently enumerates 48 complex Gaussian variance
calculations. The second checks 64 complete score-array event counts, 192
binomial critical values, 16 calibration tails and 32 independently reconstructed
first-view polynomials, using direct physical-cell means. One event count
changes under a 1e-8 threshold perturbation; the minimum distance is 8.60e-9.
These are numerical checks, not validated interval arithmetic. Fifteen new
tests pass across the three separately logged test runs; the earlier full-suite
count is historical and has not been relabeled as a new full-suite run.

To rerun original computations, use a clean checkout at the source commit
recorded below, restore its dependencies, and leave its output directory
absent. The runners guard source equality and refuse to overwrite outcomes.

| Computation | Script | Source commit |
|---|---|---|
| All 24 moment-hull cases | `scripts/probe_uq_bispectrum.py` | `baf6b29` |
| All 160 continuous-pose searches | `scripts/probe_uq_bispectrum_poses.py` | `66c08c2` |
| Analytic global-cover cost gate | `scripts/probe_uq_moment_global_bound.py` | `5c0f406` |
| Independent bounded-view Monte Carlo | `scripts/probe_uq_moment_mc.py` | `ada3c36` |

The Monte Carlo study uses two independent 131,072-draw stages per nonzero
stack, fresh continuous Haar orientations and physical-cell projections. It
does not sample the old finite view catalog. The four zero 10076 contrasts
are retained without new simulations. Scores, seeds and RNG states are saved.
Correct-null controls and all 192 projected rejection probabilities remain
in the CSV. Projections and their pointwise intervals are not observed
repeated-dataset power or experimental density coverage. The viewing-density
ratio and noise/amplitude assumptions are explicit, uncalibrated inputs.
