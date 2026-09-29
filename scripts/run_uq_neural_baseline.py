#!/usr/bin/env python3
"""Train actual stock neural cryoDRGN on audited development exposure halves.

This is fixed-pose homogeneous training, not ab initio or posterior uncertainty.
It establishes a reconstruction baseline whose sampling uncertainty/validation
must be evaluated separately. All epochs and architecture choices are recorded.
"""
import argparse
import csv
import hashlib
import json
import os
import pickle
import subprocess
import time
from pathlib import Path
import mrcfile
import numpy as np
import cryodrgn
from fourier_splats.fsc import fsc,resolution
from fourier_splats.physics import fft_volume_center

ROOT=Path(__file__).resolve().parents[1]


def run(args,dataset):
    data=ROOT/'data'/dataset;source=data/'cryodrgn';out=ROOT/'results/uncertainty/development/neural-reconstruction'/dataset
    out.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'research/uncertainty/splits'/f'{dataset}.csv').open()))
    manifest=json.loads((data/'manifest.json').read_text());halves=[];records=[]
    for half in [0,1]:
        ids=np.array([int(r['output_index']) for r in rows if r['split']==f'inference_half{half}'])
        ind=data/f'uq-neural-half{half}.pkl'
        with ind.open('wb') as f:pickle.dump(ids,f)
        dest=out/f'half{half}'
        command=[str(ROOT/'.venv/bin/cryodrgn'),'train_nn',str(source/'particles.mrcs'),'--poses',str(source/'poses.pkl'),
                 '--ctf',str(source/'ctf.pkl'),'--ind',str(ind),'-o',str(dest),'--layers','3','--dim',str(args.dim),
                 '--num-epochs',str(args.epochs),'--checkpoint','5','--batch-size','8','--no-amp',
                 '--seed',str(args.seed+half),'--shuffle-seed',str(args.seed+half)]
        if manifest['data_sign']==1:command.append('--uninvert-data')
        start=time.perf_counter();log=ROOT/'logs/uncertainty'/f'neural-{dataset}-half{half}.log'
        with log.open('w') as f:
            result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,
                                  env={**os.environ,'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4'})
        record={'half':half,'particles':len(ids),'source_groups':len({r['source_group'] for r in rows if r['split']==f'inference_half{half}'}),
                'indices_sha256':hashlib.sha256(ind.read_bytes()).hexdigest(),'command':command,'returncode':result.returncode,
                'seconds':time.perf_counter()-start}
        records.append(record)
        (out/'runs.json').write_text(json.dumps({'stage':'development fixed-pose actual neural training; convergence not yet assessed',
            'dataset':dataset,'cryodrgn_version':cryodrgn.__version__,'config':vars(args),'runs':records},indent=2)+'\n')
        if result.returncode:raise RuntimeError(f'Neural training failed: {log}')
        with mrcfile.open(dest/'reconstruct.mrc') as m:halves.append(fft_volume_center(m.data.copy()))
        print(dataset,half,record['seconds'],flush=True)
    curve=fsc(*halves,manifest['pixel_size_A'])
    np.savetxt(out/'neural-half-fsc.csv',curve,delimiter=',',header='shell,frequency_inverse_A,fsc,count',comments='')
    metrics={'stage':'development reconstruction; supplied consensus poses; half FSC is conditional',
             'dataset':dataset,'architecture':{'layers':3,'dim':args.dim,'domain':'fourier','positional_encoding':'gaussian'},
             'epochs':args.epochs,'half_fsc_resolution':resolution(curve),
             'preprocessing':'v0.1.0 per-image background normalization, stock cryoDRGN window, per-half normalization',
             'training_convergence_assessed':False}
    (out/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--epochs',type=int,default=20)
    p.add_argument('--dim',type=int,default=256);p.add_argument('--seed',type=int,default=609322);args=p.parse_args()
    for dataset in args.datasets.split(','):run(args,dataset)
