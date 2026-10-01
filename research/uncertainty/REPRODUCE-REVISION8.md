# Reproducing candidate-derived score development

The v0.7.5-dev increment supplies the 32-page paper, the three actual fitted
Gaussian candidate maps, all candidate-derived simulation arrays and results,
and their checks. It is a conditional known-simulator study, not a new
experimental reconstruction or validated regional-density confidence method.

Restore the tagged source and pinned environment. The v0.7.2-dev and v0.7.3-dev
archives supply the inherited frequency/transfer fixture and physical operator
inputs; v0.7.4-dev supplies the previous comparison studies and paper context.
Verify each archive with its adjacent manifest before extraction at repo root.
This increment also includes all three original MRC inputs, with hashes in each
run's summary. Third-party PDFs are not redistributed.

```bash
.venv/bin/python scripts/verify_uq_release_artifacts.py output/artifacts/v0.7.5-dev-candidate-scores.tar.gz
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_candidate_scores.py --output /tmp/candidate-score-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_candidate_score.py
bash scripts/build_paper.sh
```

Choose an unused verifier output path. The independent verifier checks region
selection via FFT convolution, score directions via separate feature arithmetic
and least-squares projection, all saved events and minimal critical counts,
every projection and actual repeated-group vector. It also replays first-view
physical sums, 768 calibration polynomials and 180 held-out direct scores.
It does not independently replay all physical views or establish validated
arithmetic. Three score-design tests pass; the broader suite is historical.

For new simulations, use source `932ea23` with
`scripts/probe_uq_candidate_scores.py --dataset 10028` (and 10049/10076),
restored input dependencies and absent designated output directories. The
classical-comparator addendum was frozen at `f10f4a7` while calibration was
running and before held-out outcomes were computed; run
`scripts/analyze_uq_candidate_scores.py` at that source after all stacks finish.
The scripts check committed-source equality and preserve every attempt.
Seeds, batch sizes, RNG states, input/source hashes and timings are saved.
The three run times overlap and must not be summed as wall time.

All seven calibration methods are separate tests; there is no unadjusted
minimum across scores or methods. Projections assume fresh independent Haar
views and the supplied amplitude/noise/transfer model. Actual 128-group tests
hold calibration fixed and share draws across alternatives and procedures.
They cannot be pooled as independent calibration repetitions. The largest
correct-null count 7/128 has pointwise exact interval [.0223,.1094].

The regions and directions use candidates only, but the deletion family is
hand-specified. Whole-candidate rejection cannot localize the discrepancy;
non-rejection is not a regional occupancy confidence statement. Viewing-density
caps and experimental noise calibration remain unvalidated inputs. All three
full Fable reviews remain rejections; this increment has no new full review.
