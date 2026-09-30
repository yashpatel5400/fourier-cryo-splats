#!/usr/bin/env python3
"""Phase A: declared independent-pilot registration of deposited map frames."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_map_alignment import align_map_to_pilot,apply_map_alignment
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/REFERENCE-REGISTRATION-PROTOCOL.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=[Path(__file__),PROTOCOL,ROOT/'src/fourier_splats/uq_map_alignment.py',
           ROOT/'tests/test_uq_map_alignment.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Commit phase-A source before registration')
    out=BASE/'reference-registration-v1'
    if out.exists():raise RuntimeError('Preserve all previous alignments')
    out.mkdir();summaries=[]
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        start=time.perf_counter();path=out/f'{ds}.json'
        record=dict(complete=False,dataset=ds,emd=emd,
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]))
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            pilot_path=BASE/f'representation/{ds}/real_particles-spacing-2.0.npz'
            ref_path=ROOT/f'data/uncertainty/references/emd_{emd}.map'
            ck=np.load(pilot_path);reference=VoxelReference.from_mrc(ref_path,box=64)
            grid=np.arange(-32,32)/64.;z,y,x=np.meshgrid(grid,grid,grid,indexing='ij')
            xyz=np.column_stack([x.ravel(),y.ravel(),z.ravel()]);mask=np.sum(xyz*xyz,axis=1)<=.35**2
            pilot=np.zeros(64**3);active=np.flatnonzero(mask)
            for first in range(0,len(active),2048):
                ids=active[first:first+2048]
                pilot[ids]=density_functionals(xyz[ids],ck['centers'],float(ck['sigma']))@ck['coefficients']
            pilot=pilot.reshape((64,)*3);pixel=reference.field_A/64
            alignment=align_map_to_pilot(reference.volume,pilot,pixel)
            aligned=apply_map_alignment(reference.volume,alignment['matrix'],alignment['offset'])
            array_path=out/f'{ds}-arrays.npz'
            np.savez_compressed(array_path,reference_original=reference.volume,reference_registered=aligned,
                pilot=pilot,matrix=alignment['matrix'],offset=alignment['offset'])
            record.update(complete=True,scientific_run_complete=True,field_A=reference.field_A,
                input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [pilot_path,ref_path]},
                alignment={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in alignment.items()},
                arrays_file=str(array_path.relative_to(ROOT)),arrays_sha256=sha(array_path),
                seconds=time.perf_counter()-start,
                scope='Pilot-only local reference registration. No inference maps, intervals or FSC read or optimized.')
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error),seconds=time.perf_counter()-start)
        save();summaries.append(dict(dataset=ds,record_sha256=sha(path),
            scientific_run_complete=record['scientific_run_complete'],seconds=record['seconds']))
        print(ds,summaries[-1],flush=True)
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summaries,
        scope='All phase-A outcomes retained. Publish transforms before phase-B map comparison.'),indent=2)+'\n')


if __name__=='__main__':main()
