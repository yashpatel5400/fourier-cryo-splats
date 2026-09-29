#!/usr/bin/env python3
"""Continue saved stock neural fits without replacing initial run provenance."""
import argparse,json,os,subprocess,time,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--from-epoch',type=int,default=20)
p.add_argument('--to-epoch',type=int,default=60);args=p.parse_args()
for dataset in args.datasets.split(','):
    base=ROOT/'results/uncertainty/development/neural-reconstruction'/dataset
    records=[]
    for half in [0,1]:
        folder=base/f'half{half}';source=folder/f'weights.{args.from_epoch}.pkl'
        if not source.exists():raise FileNotFoundError(source)
        if (folder/f'weights.{args.to_epoch}.pkl').exists():raise RuntimeError('Refusing to replace an existing endpoint checkpoint')
        shutil.copyfile(folder/'config.yaml',folder/f'config.before-resume-{args.to_epoch}.yaml')
        prior=json.loads((base/'runs.json').read_text())['runs'][half]
        command=prior['command'].copy();command[command.index('--num-epochs')+1]=str(args.to_epoch)
        command[command.index('--checkpoint')+1]='10';command+=['--load',str(source)]
        start=time.perf_counter();log=ROOT/'logs/uncertainty'/f'neural-extend{args.to_epoch}-{dataset}-half{half}.log'
        with log.open('w') as f:
            result=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,
                env={**os.environ,'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4'})
        records.append({'half':half,'command':command,'seconds':time.perf_counter()-start,'returncode':result.returncode,
            'from_epoch':args.from_epoch,'to_epoch':args.to_epoch,'note':'Optimizer checkpoint resumed; stock training RNG reinitialized from recorded seed.'})
        (base/f'extension-{args.to_epoch}.json').write_text(json.dumps(records,indent=2)+'\n')
        if result.returncode:raise RuntimeError(f'Training failed: {log}')
        print(dataset,half,records[-1]['seconds'],flush=True)
    metrics=json.loads((base/'metrics.json').read_text());metrics['epochs']=args.to_epoch
    metrics['training_convergence_assessed']=False
    (base/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
