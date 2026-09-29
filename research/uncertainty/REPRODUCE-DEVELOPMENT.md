# Reproduce the uncertainty development comparisons

These commands reproduce development studies, not an untouched confirmatory
benchmark. Start from the project environment, downloaded particle selections,
reference maps and metadata described in the repository README and provenance.
Install uncertainty dependencies and use the macOS FINUFFT build in COMPUTE.md.
Run from the repository root with its virtual environment active. The existing
scripts for splits, reference downloads, representation fits and initial audits
remain prerequisites; their recorded commands/configurations are in the result
files and development log.

```sh
export OPENBLAS_NUM_THREADS=4
pytest -q
for target in contrast center; do
  python scripts/benchmark_uq_ambient_pose.py --target "$target" \
    --quadratic-pose --gaussian-baselines --maxiter 500 \
    --output "ambient-quadratic-$target.json"
done
for reg in 0.001 0.01 0.1 1.0; do
  python scripts/benchmark_uq_group_bootstrap.py --replicates 1000 \
    --bootstrap 999 --regularization "$reg" --output "group-bootstrap-reg-$reg"
done
python scripts/audit_uq_grid_refinement.py
for target in center contrast; do
  for angle in 0.5 2.0; do
    python scripts/audit_uq_pose_grid.py --backend direct --target "$target" \
      --angle "$angle" --output "pose-grid-$target-$angle"
    python scripts/stress_uq_nonlinear_bias.py --steps 150 --restarts 4 \
      --source "pose-grid-$target-$angle" --output "nonlinear-stress-$target-$angle"
    python scripts/audit_uq_shared_density.py --source "pose-grid-$target-$angle" \
      --stress "nonlinear-stress-$target-$angle" --output "shared-density-$target-$angle"
  done
done
python scripts/summarize_uq_comparisons.py
python scripts/write_uq_results.py
python scripts/write_uq_citations.py
bash scripts/build_paper.sh
```

The source-group reconstruction runner, neural training/extension commands and
full checkpoint/evaluation records are under `group-reconstruction`,
`neural-reconstruction`, and `reconstruction-comparison`. Recreating MRCs and
weights requires executing the runs; large local arrays/checkpoints and
third-party source PDFs are intentionally not committed to Git. JSON records,
configuration files, FSC curves, bibliography metadata, figures and the PDF
are committed. The original acquisition and map hashes remain in provenance.

Computational availability is not scientific completion. Higher-bandwidth,
class-violation and fresh frozen studies are recorded separately as they finish.
