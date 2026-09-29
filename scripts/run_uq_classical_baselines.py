#!/usr/bin/env python3
"""Run classical voxel/Gaussian fits on the neural baseline's exposure halves."""
import os
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for dataset in ['10028','10049','10076']:
    command=[str(ROOT/'.venv/bin/python'),str(ROOT/'scripts/run_experiment.py'),dataset,
             '--samples','2048','--window','--reg','1','--iterations','120','--group-splits',
             '--tag','uncertainty/development/group-reconstruction']
    with (ROOT/'logs/uncertainty'/f'classical-{dataset}.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,
                       env={**os.environ,'OPENBLAS_NUM_THREADS':'4','OMP_NUM_THREADS':'4'})
    print(dataset,'classical fits completed',flush=True)
