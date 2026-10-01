#!/usr/bin/env python3
"""Replay law bins and projected residuals independently; direct expected scores."""
import hashlib
import json
import pickle
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.spatial.transform import Rotation

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def explicit_bins(r,b):
    # scipy is used only to orthogonalize stored float32 matrices, not for Euler angles.
    r=Rotation.from_matrix(r).as_matrix()
    beta=np.arccos(np.clip(r[:,2,2],-1,1))
    alpha=np.arctan2(r[:,1,2],r[:,0,2])
    gamma=np.arctan2(r[:,2,1],-r[:,2,0])
    low=np.sin(beta)<1e-12
    zero=low&(r[:,2,2]>=0);pi=low&~zero
    alpha[zero]=np.arctan2(r[zero,1,0],r[zero,0,0])
    alpha[pi]=np.arctan2(-r[pi,1,0],-r[pi,0,0]);gamma[low]=0
    unit=np.column_stack([(alpha%(2*np.pi))/(2*np.pi),(np.cos(beta)+1)/2,(gamma%(2*np.pi))/(2*np.pi)])
    ix=np.minimum(np.maximum(np.floor(b*unit).astype(int),0),b-1)
    return ix[:,0]*b*b+ix[:,1]*b+ix[:,2]


def expected_score(energy,weight_a,weight_b,t):
    # Different quadrature implementation and direct density-ratio score under
    # both Gaussian states, rather than the J_t identity used in the runner.
    z,w=hermgauss(160);z=z*np.sqrt(2);w=w/np.sqrt(np.pi)
    score=0.
    for begin in range(0,len(energy),2048):
        s=energy[begin:begin+2048,None]
        for sign,law,pop in [(1,weight_a,.75),(-1,weight_b,.25)]:
            logratio=np.sqrt(s)*z+sign*s/2
            ratio=np.exp(logratio)
            conditional=((ratio-1)/(1-t+t*ratio))@w
            score+=pop*(law[begin:begin+len(s)]@conditional)
    return float(score)


def main():
    for ds in ['10028','10049','10076']:
        directory=BASE/'population-materiality-v1'/ds
        out=directory/'independent-check.json'
        if out.exists():raise ValueError('Preserve every verification attempt')
        j=json.loads((directory/'summary.json').read_text());assert j['complete']
        assert sha(directory/'arrays.npz')==j['array']['sha256']
        with np.load(directory/'arrays.npz') as f:archive={k:f[k] for k in f.files}
        with np.load(BASE/'candidate-fisher-score-v1'/ds/'arrays.npz') as f:
            m=f['training_candidate_means'];g=f['training_region_means']
        mr=np.concatenate([m.real,m.imag],axis=1);gr=np.concatenate([g.real,g.imag],axis=1)
        energy=np.einsum('ij,ij->i',gr,gr)
        errors=[]
        for point,name in [(mr,'amplitude_residual_full'),(mr-gr,'amplitude_residual_deleted')]:
            coefficient=np.einsum('ij,ij->i',point,gr)/np.einsum('ij,ij->i',point,point)
            residual=gr-coefficient[:,None]*point
            actual=np.einsum('ij,ij->i',residual,residual)
            error=float(abs(actual-archive[name]).max());errors.append(error)
            np.testing.assert_allclose(actual,archive[name],atol=1e-12,rtol=1e-12)
        np.testing.assert_allclose(energy,archive['energy'],atol=1e-12,rtol=1e-12)
        source=ROOT/f'background/cryodrgn_empiar/empiar{ds}/inputs'
        data=np.load(next(source.glob('*.cs')),allow_pickle=False)
        poses=pickle.load((source/'poses.pkl').open('rb'))[0]
        masks={'all_metadata':np.ones(len(data),bool),'source_half_0':data['alignments3D/split']==0,'source_half_1':data['alignments3D/split']==1}
        if (source/'filtered.ind.pkl').exists():
            keep=np.zeros(len(data),bool);keep[np.asarray(pickle.load((source/'filtered.ind.pkl').open('rb')),int)]=True
            masks.update(published_filter=keep,published_filter_complement=~keep)
        record=dict(complete=False,dataset=ds,source_sha256=sha(Path(__file__)),summary_sha256=sha(directory/'summary.json'),
                    projected_residual_maximum_errors=errors,grids=[],expected_scores=[])
        primary_passes=[]
        for grid in j['grids']:
            assert grid['available'];b=grid['bins_per_axis'];size=b**3
            bins=explicit_bins(archive['rotations'],b)
            np.testing.assert_array_equal(bins,archive[f'training_bins_{b}'])
            rbin=explicit_bins(poses,b);counts=np.bincount(bins,minlength=size)
            weights={}
            for name,law in grid['laws'].items():
                c=np.bincount(rbin[masks[name]],minlength=size)
                np.testing.assert_array_equal(c,law['counts'])
                weights[name]=c[bins]/c.sum()/counts[bins]
                for key,summ in law['summaries'].items():
                    actual=weights[name]@archive[key]
                    np.testing.assert_allclose(actual,summ['mean'],rtol=1e-12,atol=1e-12)
            record['grids'].append(dict(bins_per_axis=b,training_rotations_checked=len(bins),metadata_rotations_checked=len(rbin),law_count=len(weights)))
            for pair in grid['pairs']:
                wa=weights[pair['state_a_law']];wb=weights[pair['state_b_law']]
                value=expected_score(energy,wa,wb,pair['fraction'])
                if abs(value)>1e-7:raise ArithmeticError('Independent expected mixture score fails')
                harmonic=(.75*wa+.25*wb)@archive['inverse_amplitude_residual_full']
                np.testing.assert_allclose(harmonic,pair['amplitude_tangent_harmonic'],atol=1e-12,rtol=1e-12)
                passed=bool(abs(pair['fraction']-.75)>=.02 and harmonic<=8.8125 and pair['quadrature_check_passed'])
                assert passed==pair['screen_passed']
                if b==8:primary_passes.append(passed)
                record['expected_scores'].append(dict(grid=b,state_a_law=pair['state_a_law'],state_b_law=pair['state_b_law'],score_160=value))
        assert any(primary_passes)==j['primary_screen_passed']
        record.update(complete=True,primary_screen_passed=any(primary_passes))
        out.write_text(json.dumps(record,indent=2)+'\n');print('VERIFIED',ds,len(record['expected_scores']),flush=True)


if __name__=='__main__':main()
