#!/usr/bin/env python3
"""Cost gate for an explicit global rotation certificate; evaluates no grid."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_moment_rotation_bound import (
    cell_absolute_radial_moments,contrast_rotation_curvature,euler_cover_for_remainder)

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
SOURCES=['scripts/probe_uq_moment_global_bound.py','src/fourier_splats/uq_moment_rotation_bound.py',
    'research/uncertainty/paired-power-v1/BISPECTRUM-GLOBAL-BOUND-PROTOCOL.md']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol first')
    out=BASE/'bispectrum-global-bound-v1'
    if out.exists():raise ValueError('Preserve every attempt')
    out.mkdir()
    result=dict(complete=False,grid_scores_evaluated=0,cases=[],input_hashes={},
        sources={n:sha(ROOT/n) for n in SOURCES},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
        return json.loads(p.read_text())
    prior=read(BASE/'bispectrum-hull-v1/summary.json')
    poses=read(BASE/'bispectrum-adversarial-pose-v1/summary.json')
    for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
        mp=ROOT/f'data/uncertainty/references/emd_{emd}.map'
        ap=ROOT/prior['arrays'][ds]['path']
        assert sha(ap)==prior['arrays'][ds]['sha256']
        for path in [mp,ap]:result['input_hashes'][str(path.relative_to(ROOT))]=sha(path)
        geometry=read(BASE/f'mixture-validation-preflight-v2/{ds}.json')
        with np.load(ap) as f:data={k:f[k].copy() for k in f.files}
        rho=VoxelReference.from_mrc(mp,box=64).volume.ravel();rho/=np.linalg.norm(rho)
        x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij')
        xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
        center=np.array(geometry['target']['center_fraction_field'])
        density=rho*(1-np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/geometry['field_A'])**2))
        moments=cell_absolute_radial_moments(density,64)
        for case in poses['cases']:
            if case['dataset']!=ds:continue
            key=case['key'];direction=data[key+'_direction']
            row=dict(dataset=ds,key=key,absolute_radial_moments=moments.tolist())
            if not np.linalg.norm(direction):
                row['skipped']='Zero direction; no separator';result['cases'].append(row);continue
            triads=data['triads'] if key.startswith('power_bispectrum') else np.empty((0,3),int)
            curvature=contrast_rotation_curvature(data['q'],data['transfer'],triads,direction,
                moments,max(abs(a) for a in case['amplitudes']))
            margin=case['remaining_mean_separation']
            row.update(curvature=curvature,sampled_margin=margin)
            if margin>0:
                row['proposed_grid']=euler_cover_for_remainder(curvature['total_curvature'],margin/2)
                row['exceeds_10_million_gate']=row['proposed_grid']['points']>10_000_000
            else:row['skipped_grid']='Already nonpositive sampled separation'
            result['cases'].append(row)
            print(ds,key,'H',curvature['total_curvature'],'grid',row.get('proposed_grid',{}).get('points'),'sampled margin',margin,flush=True)
    result['complete']=True
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('COMPLETE',len(result['cases']),'cases; no rotation grids evaluated',flush=True)


if __name__=='__main__':main()
