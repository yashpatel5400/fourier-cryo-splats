#!/usr/bin/env python3
"""Frozen known-pose screen on existing means and explicitly defined view proxies."""
import argparse
import hashlib
import json
import pickle
import subprocess
import time
import warnings
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_population_screen import (
    amplitude_tangent_information, pseudo_true_fraction, histogram_law_weights,
    stratified_linear_summary)
from probe_uq_bispectrum_poses import direct_cells

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
SOURCES = ['scripts/screen_population_materiality.py',
           'src/fourier_splats/uq_population_screen.py',
           'research/uncertainty/POPULATION-MATERIALITY-PROTOCOL.md',
           'research/uncertainty/POPULATION-SCREEN-ALGEBRA.md']


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def bins_for(rotations, b):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        angles = Rotation.from_matrix(rotations).as_euler('ZYZ')
    u = np.column_stack([(angles[:,0]%(2*np.pi))/(2*np.pi),
                         (np.cos(angles[:,1])+1)/2,
                         (angles[:,2]%(2*np.pi))/(2*np.pi)])
    ix = np.clip((u*b).astype(int),0,b-1)
    return ix[:,0]*b*b + ix[:,1]*b + ix[:,2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', choices=['10028','10049','10076'], required=True)
    ds = parser.parse_args().dataset
    for name in SOURCES:
        if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT) != (ROOT/name).read_bytes():
            raise ValueError('Commit source and protocol before execution')
    out = BASE/'population-materiality-v1'/ds
    if out.exists():
        raise ValueError('Preserve every attempt')
    out.mkdir(parents=True)
    start = time.perf_counter()
    prior = BASE/'candidate-fisher-score-v1'/ds
    source = ROOT/f'background/cryodrgn_empiar/empiar{ds}/inputs'
    cs = next(source.glob('*.cs'))
    files = [prior/'arrays.npz', prior/'summary.json', cs, source/'poses.pkl']
    if (source/'filtered.ind.pkl').exists():
        files.append(source/'filtered.ind.pkl')
    old = json.loads((prior/'summary.json').read_text())
    assert old['complete'] and sha(prior/'arrays.npz') == old['array']['sha256']
    result = dict(complete=False,dataset=ds,source_hashes={s:sha(ROOT/s) for s in SOURCES},
                  git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  inputs={str(p.relative_to(ROOT)):sha(p) for p in files},
                  interpretation='Known-pose white-noise diagnostic; law proxies are not biological class labels',
                  population=.75,nominal_particles=10000,primary_grid=8,
                  physical_checks=[],grids=[])
    def save():
        (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    save()
    with np.load(prior/'arrays.npz') as f:
        m = f['training_candidate_means']; g = f['training_region_means']
        q=f['q']; transfer=f['transfer']; rho=f['candidate_density']; region=f['region_density']
        first=f['stage0_first_rotation']
    assert m.shape == g.shape == (65536,220)
    rng=np.random.default_rng(261017+int(ds)); rotations=[]
    for _ in range(128):
        rotations.append(Rotation.random(512,random_state=rng).as_matrix())
        rng.normal(size=(512,220)); rng.normal(size=(512,220))
    rotations=np.concatenate(rotations)
    np.testing.assert_allclose(rotations[0],first,atol=1e-14,rtol=0)
    plane=np.pad(q,((0,0),(0,1)))
    for index in [0,511,512,32767,32768,65535]:
        for label,density,means in [('full',rho,m),('region',region,g)]:
            direct=direct_cells(density.ravel(),plane@rotations[index])*transfer
            error=float(abs(direct-means[index]).max())
            if error>1e-8:
                raise ArithmeticError('Later-batch physical replay failed')
            result['physical_checks'].append(dict(index=index,map=label,maximum_error=error))
    s,perp,perp_b=amplitude_tangent_information(m,g)
    inverse=np.divide(1,s,out=np.full_like(s,np.inf),where=s>0)
    inverse_perp=np.divide(1,perp,out=np.full_like(s,np.inf),where=perp>0)
    inverse_perp_b=np.divide(1,perp_b,out=np.full_like(s,np.inf),where=perp_b>0)
    values=dict(energy=s,amplitude_residual_full=perp,amplitude_residual_deleted=perp_b,
                inverse_energy=inverse,inverse_amplitude_residual_full=inverse_perp,
                inverse_amplitude_residual_deleted=inverse_perp_b)
    result['haar']={k:dict(mean=float(v.mean()) if np.isfinite(v).all() else None,
        mc_standard_error=float(v.std(ddof=1)/np.sqrt(len(v))) if np.isfinite(v).all() else None,
        quantiles={str(q):float(np.quantile(v,q)) for q in [0,.01,.05,.5,.95,.99,1]} if np.isfinite(v).all() else None)
        for k,v in values.items()}
    data=np.load(cs,allow_pickle=False)
    recorded=np.asarray(pickle.load((source/'poses.pkl').open('rb'))[0],float)
    assert len(recorded)==len(data)
    np.testing.assert_allclose(Rotation.from_rotvec(data['alignments3D/pose']).as_matrix().transpose(0,2,1),recorded,atol=5e-6,rtol=0)
    result['metadata_fields']={name:dict(unique=np.unique(data[name]).tolist()) for name in ['alignments3D/class','alignments3D/class_posterior']}
    masks={'all_metadata':np.ones(len(data),bool), 'source_half_0':data['alignments3D/split']==0,
           'source_half_1':data['alignments3D/split']==1}
    pairs=[('source_half_0','source_half_1')]
    if (source/'filtered.ind.pkl').exists():
        keep=np.zeros(len(data),bool)
        keep[np.asarray(pickle.load((source/'filtered.ind.pkl').open('rb')),int)]=True
        masks.update(published_filter=keep,published_filter_complement=~keep)
        pairs.append(('published_filter','published_filter_complement'))
    archive=dict(rotations=rotations,**values)
    for b in [4,8,12]:
        train_bins=bins_for(rotations,b); recorded_bins=bins_for(recorded,b)
        archive[f'training_bins_{b}']=train_bins
        counts={name:np.bincount(recorded_bins[mask],minlength=b**3) for name,mask in masks.items()}
        row=dict(bins_per_axis=b,cells=b**3,laws={},pairs=[],available=True)
        result['grids'].append(row)
        if any(np.any((c>0)&(np.bincount(train_bins,minlength=b**3)==0)) for c in counts.values()):
            row.update(available=False,reason='Occupied metadata cell with no training rotation'); save(); continue
        weights={k:histogram_law_weights(train_bins,c) for k,c in counts.items()}
        for name,c in counts.items():
            archive[f'counts_{b}_{name}']=c
            row['laws'][name]=dict(particles=int(c.sum()),counts=c.tolist(),
                summaries={k:stratified_linear_summary(v,train_bins,c) for k,v in values.items()})
        for left,right in pairs:
            for a,z in [(left,right),(right,left)]:
                estimate=pseudo_true_fraction(s,weights[a],weights[z])
                harmonic=float((.75*weights[a]+.25*weights[z])@inverse_perp) if np.isfinite(inverse_perp).all() else None
                known_harmonic=float((.75*weights[a]+.25*weights[z])@inverse) if np.isfinite(inverse).all() else None
                passed=bool(estimate['quadrature_check_passed'] and abs(estimate['bias'])>=.02 and harmonic is not None and harmonic<=8.8125)
                row['pairs'].append(dict(state_a_law=a,state_b_law=z,**estimate,
                    amplitude_tangent_harmonic=harmonic,
                    amplitude_tangent_surrogate_se=np.sqrt((.1875+harmonic)/10000) if harmonic is not None else None,
                    known_amplitude_harmonic=known_harmonic,
                    known_amplitude_oracle_se=np.sqrt((.1875+known_harmonic)/10000) if known_harmonic is not None else None,
                    screen_passed=passed))
                save(); print(ds,'grid',b,'pair',a,z,'complete',flush=True)
        save()
    primary=next(r for r in result['grids'] if r['bins_per_axis']==8)
    result['primary_screen_passed']=bool(primary['available'] and any(r['screen_passed'] for r in primary['pairs']))
    path=out/'arrays.npz';np.savez_compressed(path,**archive)
    result.update(complete=True,seconds=time.perf_counter()-start,array=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bytes=path.stat().st_size))
    save();print('COMPLETE',ds,result['primary_screen_passed'],flush=True)


if __name__=='__main__':
    main()
