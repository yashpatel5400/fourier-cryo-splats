# Reproducing conditional-noise and preferred-view development

The v0.7.4-dev increment contains the 30-page ICML-format paper's new
conditional-noise study, classical risk comparators, preferred-view controls,
authentic focused Fable audit and response, and all numerical arrays. This
is not a favorable full review or experimental density coverage.

Use the tagged source and pinned Python environment. Restore the v0.7.2-dev
method-gate bundle (maps, geometry and original orbit arrays), the v0.7.3-dev
invariant-moment bundle (frozen directions and earlier calibration), and this
increment. Verify every archive against its adjacent manifest before extracting
at the repository root. Third-party PDFs are not redistributed; primary URLs,
reading scope and downloaded-file hashes are retained.

```bash
.venv/bin/python scripts/verify_uq_release_artifacts.py output/artifacts/v0.7.4-dev-viewing-uncertainty.tar.gz
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_view_variance.py --output /tmp/view-variance-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_view_variance_sample.py --output /tmp/view-variance-one-percent-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_preferred_views.py --output /tmp/preferred-view-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_view_risk.py --output /tmp/view-risk-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_view_variance.py tests/test_preferred_views.py tests/test_view_risk.py tests/test_moment_mc.py
bash scripts/build_paper.sh
```

Choose unused output paths; verifiers refuse to overwrite records. The 1%
replay independently evaluates separable physical-cell Fourier sums and
moment polynomials at verifier-selected views, with companion-root amplitude
profiling. It checks 167,936 cubics and 10,496 grouped/individual counts, not
all views. Other replays check all saved event counts, critical values and
repeated-group vectors, first-view physical scores, non-antipodal frequencies
and exact transfer-profile provenance. These are floating-point checks.
Nineteen targeted tests pass; the earlier full-suite count is historical.

For new simulations, use a clean checkout at the recorded source commit,
restore its dependencies, and leave the designated output directory absent:

| Computation | Script | Source |
|---|---|---|
| Conditional noise, both stacks | `probe_uq_view_variance.py --dataset 10028` and `--dataset 10049` | `96934a8` |
| Nine preferred laws per stack | `probe_uq_preferred_views.py --dataset 10028` and `--dataset 10049` | `1dc3a70` |
| Four classical saved-array comparators | `probe_uq_view_risk.py` | `82f8403` |

Scripts are in `scripts/` and run with the same Python environment variables
shown above. Each runner requires committed source equality and preserves all
attempts. Seeds, initial/final RNG states, source hashes and input hashes are
recorded in every summary. The first two computations were run as separate,
overlapping processes per stack; summed run seconds are not elapsed wall time.

Grouped methods require amplitude conditionally independent of noise given
view. All tests assume the stated noise/CTF simulator and independent particles.
Viewing caps are supplied, not experimentally estimated. Actual 64-group
outcomes hold calibration fixed; binomial power projections are distinct from
those observed outcomes. No optimized allocation, calibration-repeat,
extremal-view-law or new three-stack positive result is claimed. The full
Fable verdicts remain rejections; this increment includes only a focused audit.
