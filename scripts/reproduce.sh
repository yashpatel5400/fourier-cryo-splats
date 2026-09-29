#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=4
python -m pytest -q
python scripts/validate_synthetic.py
for dataset in 10028 10049 10076; do
  if [ ! -f "data/$dataset/manifest.json" ]; then
    python scripts/download_data.py "$dataset" --resume
  fi
  python scripts/run_experiment.py "$dataset" --samples 256 --tag main
  python scripts/run_cryodrgn_baseline.py "$dataset" --tag main
  python scripts/run_experiment.py "$dataset" --samples 2048 --window --reg 1 --tag full-reg1
  python scripts/run_experiment.py "$dataset" --box 32 --samples 256 --window --reg 1 --phase-randomize --tag noise
done
for factor in 0.1 10; do
  python scripts/run_experiment.py 10049 --samples 2048 --window --reg "$factor" --tag "full-reg$factor"
done
python scripts/fit_adaptive_demo.py --dataset 10049 --steps 1500
python scripts/select_final.py
python scripts/package_provenance.py
python scripts/check_truncation.py
python scripts/plot_results.py --tag final
python scripts/build_report.py
