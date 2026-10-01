# Reproducing the Fisher and replica-allocation comparisons

This increment has four numerical archives: one common source/report bundle
and one for each stack (10028, 10049, 10076). Restore all four at the repository
root after verifying each manifest. The per-stack archives retain full training
Fourier means and noise draws, covariance estimates, calibration envelopes,
held-out scores and all outcomes. No training data are omitted to reduce size.
All individual assets are below 2GB. The original Gaussian MRC inputs are in
v0.7.5-dev, which records its own earlier dependencies.

Use the tagged source and pinned environment. Verification output paths must
be absent; every verifier preserves earlier attempts.

```bash
.venv/bin/python scripts/verify_uq_release_artifacts.py output/artifacts/v0.7.6-dev-score-calibration*.tar.gz
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_fisher_scores.py --output /tmp/fisher-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_fisher_comparisons.py --output /tmp/fisher-comparison-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_replica_allocation.py --output /tmp/allocation-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_fisher_score.py tests/test_candidate_score.py tests/test_view_variance.py tests/test_view_risk.py tests/test_moment_mc.py tests/test_preferred_views.py
bash scripts/build_paper.sh
```

Twenty-five targeted tests pass. The first verifier recomputes all three
feature covariances with separate feature arithmetic and a full-array sample
covariance, and uses a null-space formulation to check score optimality. It
replays all training noises and terminal RNG states, 420 probability bounds,
1,680 critical values, 18,900 projections and 6,300 group vectors. The second
checks 90 joint event tables and 9,450 conservative pointwise difference
intervals. The third checks unchanged inherited scores, all new bounds and
outcomes under both allocations, and first-view physical scores/cubics. These
are ordinary numerical checks, not certified arithmetic or full operator replay.

For new simulations use source `669df5a` with
`scripts/probe_uq_fisher_scores.py --dataset 10028` (and 10049/10076), then
source `ec73892` with `scripts/probe_uq_replica_allocation.py --dataset 10028`
(and 10049/10076). Use the environment variables above and absent output
directories. The Fisher comparison-interval addendum is at `3824f60`, declared
during calibration before held-out outcomes. Scripts require committed source
equality; source/input hashes, seeds, batches and full RNG states are saved.
Per-stack runs overlap; summing their durations is not elapsed wall time.

The Fisher construction is a classical equality-constrained Rayleigh design,
with fixed 0.1 covariance shrinkage and the same training means as its matched
comparator. It is not a new optimization theorem. The allocation comparison
keeps the noise-draw count fixed (2,097,152 per stack): 32,768 views with two
32-noise groups versus 8,192 with two 128-noise groups. Fourier costs and wall
times differ. Both allocations use the same new test draws, distinct from
those previously inspected in the Fisher comparison. No cross-allocation or
cross-method minimum is presented as a combined test.

All simulator assumptions remain supplied. Tests concern the whole candidate,
not regional attribution; groups hold calibration fixed. Unconditional
calibration repetition, realistic nuisance calibration, useful small-error
sensitivity and substantive methodological novelty remain open. No new full
Fable review or acceptance verdict is part of this release.
