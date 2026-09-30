#!/usr/bin/env python3
"""Phase B of declared reference registration; never refit a transform."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import mrcfile
import numpy as np
from fourier_splats.fsc import fsc,resolution
from fourier_splats.physics import fft_volume_center,volume_from_fourier
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--dataset',choices=['10028','10049','10076'],required=True)
    parser.add_argument('--method',choices=['gaussian','voxel','neural','relion'],required=True)
    args=parser.parse_args();ds=args.dataset;method=args.method
    protocol=ROOT/'research/uncertainty/REFERENCE-REGISTRATION-PROTOCOL.md'
    # All three transforms, not just a favorable case, must precede phase B.
    for dataset in ['10028','10049','10076']:
        p=BASE/f'reference-registration-v1/{dataset}.json'
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Publish every transform before phase-B comparison')
    for p in [Path(__file__),protocol]:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit comparison code and protocol')
    out=BASE/'reference-registered-comparison-v1';out.mkdir(exist_ok=True)
    path=out/f'{ds}-{method}.json'
    if path.exists():raise RuntimeError('Preserve previous comparisons and failures')
    start=time.perf_counter();record=dict(complete=False,dataset=ds,method=method,
        source_snapshot=source_snapshot(ROOT,Path(__file__),[str(protocol.relative_to(ROOT))]),
        scope='Post-outcome registered reference agreement, pilot-selected transform. Not independent accuracy or density coverage.')
    def save():path.write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        source=BASE/f'reference-registration-v1/{ds}.json';registration=json.loads(source.read_text())
        if not registration['scientific_run_complete']:raise ValueError('Reference alignment failed')
        array_path=ROOT/registration['arrays_file'];files={str(source.relative_to(ROOT)):sha(source)}
        if sha(array_path)!=registration['arrays_sha256']:raise ValueError('Reference transform changed')
        files[str(array_path.relative_to(ROOT))]=sha(array_path);references=np.load(array_path)
        if method=='relion':
            p=BASE/f'relion-evaluation-v2/{ds}/metrics.json';r=json.loads(p.read_text())
            if not r['complete'] or r.get('skipped') or r.get('error'):raise ValueError('Completed RELION evaluation required')
            a=p.parent/'metrics-maps.npz'
            if sha(a)!=r['array_sha256']:raise ValueError('Changed RELION maps')
            files[str(p.relative_to(ROOT))]=sha(p);files[str(a.relative_to(ROOT))]=sha(a)
            volume=np.load(a)['relion'];record['relion_converged']=r['refinement_converged']
        else:
            lockpath=ROOT/'research/uncertainty/confirmation/prediction-v1/locked-models.json'
            lock=json.loads(lockpath.read_text());files[str(lockpath.relative_to(ROOT))]=sha(lockpath)
            values=[]
            for half in [0,1]:
                if method=='neural':
                    p=BASE/f'neural-reconstruction/{ds}/half{half}/reconstruct.{lock["epochs"][ds]}.mrc'
                    with mrcfile.open(p) as m:value=m.data.astype(float)
                else:
                    p=BASE/f'group-reconstruction/{ds}/{method}-half{half}.npz'
                    value=volume_from_fourier(np.load(p)['fourier'])
                if sha(p)!=lock['models'][ds]['files'][str(p.relative_to(ROOT))]:raise ValueError('Changed locked reconstruction')
                files[str(p.relative_to(ROOT))]=sha(p);values.append(value)
            volume=sum(values)/2
        curves={};metrics={}
        for name in ['original','registered']:
            ref=references['reference_'+name];curve=fsc(fft_volume_center(volume),fft_volume_center(ref),registration['field_A']/64)
            curves[name]=curve
            metrics[name]=dict(threshold_0143=resolution(curve,.143),threshold_05=resolution(curve,.5),
                mean_fsc=float(np.nanmean(curve[:,2])))
        p=out/f'{ds}-{method}-curves.npz';np.savez_compressed(p,**curves)
        record.update(complete=True,scientific_run_complete=True,seconds=time.perf_counter()-start,
            input_hashes=files,fsc=metrics,arrays_file=str(p.relative_to(ROOT)),arrays_sha256=sha(p))
    except Exception as error:
        record.update(complete=True,scientific_run_complete=False,error=repr(error),seconds=time.perf_counter()-start)
    save();print(ds,method,record.get('fsc'),record.get('error'),flush=True)


if __name__=='__main__':main()
