#!/usr/bin/env python3
"""Bounded stock cryoDRGN neural training profile, not a convergence result."""
import csv
import hashlib
import json
import os
import pickle
import re
import subprocess
import time
from pathlib import Path
import numpy as np
import cryodrgn

ROOT=Path(__file__).resolve().parents[1]
dataset='10028';seed=609320
rows=list(csv.DictReader((ROOT/'research/uncertainty/splits'/f'{dataset}.csv').open()))
indices=np.array([int(r['output_index']) for r in rows if r['split']=='pilot'])
indices=np.random.default_rng(seed).choice(indices,512,replace=False)
out=ROOT/'results/uncertainty/development/neural-runtime-probe';out.mkdir(parents=True,exist_ok=True)
ind=ROOT/'data/uncertainty/neural-runtime-probe-indices.pkl';ind.parent.mkdir(parents=True,exist_ok=True)
with ind.open('wb') as f:pickle.dump(indices,f)
source=ROOT/'data'/dataset/'cryodrgn'
command=[str(ROOT/'.venv/bin/cryodrgn'),'train_nn',str(source/'particles.mrcs'),
         '--poses',str(source/'poses.pkl'),'--ctf',str(source/'ctf.pkl'),'--ind',str(ind),'-o',str(out),
         '--layers','3','--dim','128','--num-epochs','1','--batch-size','8','--no-amp','--seed',str(seed),'--shuffle-seed',str(seed)]
manifest=json.loads((ROOT/'data'/dataset/'manifest.json').read_text())
if manifest['data_sign']==1:command.append('--uninvert-data')
log=ROOT/'logs/uncertainty/cryodrgn-neural-profile.log'
start=time.perf_counter()
with log.open('w') as f:
    completed=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4'})
result={'stage':'one-epoch timing profile only; not converged reconstruction or uncertainty baseline',
        'dataset':dataset,'particles':512,'box':64,'indices_sha256':hashlib.sha256(ind.read_bytes()).hexdigest(),
        'cryodrgn_version':cryodrgn.__version__,'command':command,'seconds':time.perf_counter()-start,
        'returncode':completed.returncode,'stock_device_policy':'CUDA if available, else CPU; this Mac uses CPU',
        'input_preprocessing':'existing v0.1.0 per-image background normalization; stock real-space window retained'}
contents=log.read_text()
def seconds(value):
    hours,minutes,seconds=map(float,value.split(':'));return hours*3600+minutes*60+seconds
epoch=re.search(r'Epoch: 1 Average loss = .*?Finished in ([0-9:.]+)',contents)
total=re.search(r'Finished in ([0-9:.]+) \([^\n]+per epoch\)',contents)
parameters=re.search(r'(\d+) parameters in model',contents)
if epoch:result['training_seconds']=seconds(epoch.group(1))
if total:result['training_and_output_seconds']=seconds(total.group(1))
if parameters:result['parameters']=int(parameters.group(1))
(out/'profile.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
raise SystemExit(completed.returncode)
