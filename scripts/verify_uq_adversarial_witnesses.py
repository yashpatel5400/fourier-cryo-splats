#!/usr/bin/env python3
"""Verify four worst pose witnesses with direct physical-cell exponential sums."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import VoxelReference

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main(output=None):
    def load(directory):
        p=BASE/directory/'summary.json';record=json.loads(p.read_text())
        row=record['arrays'] if directory.endswith('10049-v1') else record['arrays']['10049']
        ap=ROOT/row['path'];assert hashlib.sha256(ap.read_bytes()).hexdigest()==row['sha256']
        with np.load(ap) as f:arrays={k:f[k].copy() for k in f.files}
        return record,arrays
    outcomes,witnesses=load('paired-covariance-adversarial-10049-v1')
    _,source=load('paired-covariance-shifted-v1')
    geometry=json.loads((BASE/'mixture-validation-preflight-v2/10049.json').read_text())
    rho=VoxelReference.from_mrc(ROOT/'data/uncertainty/references/emd_6487.map',box=64).volume.ravel();rho/=np.linalg.norm(rho)
    centers=(np.indices((64,)*3).reshape(3,-1).T[:,::-1]+.5)/64-.5
    center=np.asarray(geometry['target']['center_fraction_field'])
    density=rho*(1-np.exp(-.5*np.sum((centers-center)**2,axis=1)/(20/geometry['field_A'])**2))
    q=source['q'];transfer=source['transfer'];plane=np.pad(q,((0,0),(0,1)));rows=[]
    for case in outcomes['cases']:
        sigma=case['shift_sd_pixels'];key=f'sigma{sigma:g}'
        index=int(np.argmax([r['verified_cell_score'] for r in case['starts']]))
        rotation=witnesses[key+'_selected_rotations'][index];shift=witnesses[key+'_selected_shifts'][index]
        k=plane@rotation;values=np.empty(len(k),complex)
        for begin in range(0,len(k),8):
            points=k[begin:begin+8]
            values[begin:begin+len(points)]=(np.exp(-2j*np.pi*(points@centers.T))@density)/64**1.5*np.prod(np.sinc(points/64),axis=1)
        observed=values*transfer*np.exp(-2j*np.pi*(q@shift)/64)
        mean=np.r_[observed.real,observed.imag];direction=source[key+'_union_direction']
        score=float(mean@direction@mean/(mean@mean));old=case['starts'][index]['verified_cell_score']
        assert score>1e-4 and abs(score-old)<1e-9
        rows.append(dict(shift_sd_pixels=sigma,start=index,direct_physical_cell_score=score,
            nufft_score=old,absolute_difference=abs(score-old)))
        print(sigma,score,abs(score-old),flush=True)
    record=dict(complete=True,witnesses=rows,
        scope='Independent direct cell-center exponential sums, including exact cell extent; ordinary floating-point verification of four positive violations.')
    path=ROOT/'provenance/uncertainty/adversarial-pose-direct-verification.json' if output is None else Path(output)
    if path.exists():raise RuntimeError('Preserve prior verification')
    path.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args();main(args.output)
