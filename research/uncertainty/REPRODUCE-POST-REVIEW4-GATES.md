# Reproducing the numerical-repair and population-screen checkpoint

The v0.7.8-dev archives are an increment to v0.7.7-dev and v0.7.6-dev. Extract
all required numerical archives at a clean repository root. Each new manifest
records exact dependency hashes; upstream acquisition metadata and particle
pixels are obtained through the existing pinned download manifests. Third-party
papers and CAHRA raw supplement pixels are not redistributed here.

The three stack bundles contain every catalogue-proposal point, family and
component label, catalogue weight, physical likelihood, observation, summary,
independent check, and all 65,536-rotation population-screen arrays. Common
evidence contains all reports, protocols, implementations, compact CSVs,
both saved-draw autopsies, acquisition-group descriptions and the focused
consultation. The exact external prompt and visible final text are preserved;
private provider reasoning is not included.

## Existing-array verification

Use the project `.venv` and `PYTHONPATH=src`. The independent checks are:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_catalog_pose_integration.py --dataset 10028
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_population_materiality.py
```

Repeat the first command for 10049 and 10076. The second covers all stacks.
Scripts intentionally refuse to overwrite an original check: for a new audit,
change only its output destination in a separate checkout and retain both
source hashes. Do not remove old failures to create a clean-looking history.
The catalogue density check independently evaluates 32 fixed draws in every
bank/image, not all density points; all integral/ESS/gate summaries replay.
Six direct physical residual checks per stack cover both maps. The population
check uses explicit rotation-matrix Euler formulas, real-coordinate projected
residuals and a distinct 160-node direct mixture-score integral.

## Frozen runners and reports

The original source/protocol commits are recorded in every summary. Their
commands are:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 .venv/bin/python scripts/probe_catalog_pose_integration.py --dataset 10028
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/screen_population_materiality.py --dataset 10028
```

Both cover one dataset per command; all three were run, with unchanged seeds,
cohorts and criteria. They refuse an occupied output directory. They require
committed source/protocol, and the original dependencies must match hashes.
Actual reruns belong in an explicitly new output version, not these immutable
attempts. Report generators likewise refuse to overwrite completed reports.

`autopsy_pose_integration.py --attempt adaptive` and `--attempt catalog` use
only saved draws and the matching generating-state truth. `report_fixed_bank_separation.py`
uses all 96 saved bank/level/deletion statistics, no new data. The score-SE
calculation has 288 independent finite-difference derivative checks. The
population and planted-identity tests are `test_population_screen.py` and
`test_planted_importance_identity.py`; six new targeted tests pass. Previously
recorded catalogue/importance tests are retained as a separate six-test check.

The public CAHRA quarter-turn check requires the officially linked supplement;
its download hash and all 432 comparison rows are retained. It is not a noisy
population validation dataset and supplies no experimental coverage result.

## Decisions and limits

The numerical branch ends after its single allowed repair, which passes.
The proposed compact-region population direction ends at its failed
materiality gate; sensitivity grids, new regions or relaxed thresholds do not
rescue it. The state/view law proxies come from source halves and filtering,
not biological labels. Nonflat corner spectra do not determine in-band pure
noise. All four full manuscript reviews reject. Neither these arithmetic
gates nor a correct conditional theorem establishes an experimentally useful
uncertainty method. The main PDF has not been updated to imply otherwise.
