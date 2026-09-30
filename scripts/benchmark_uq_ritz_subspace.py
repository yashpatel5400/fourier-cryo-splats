#!/usr/bin/env python3
"""Bounded profile of the clustered-mode eigensolver; no interval certificate."""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import finufft
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_provenance import source_snapshot
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    path=BASE/'ritz-subspace-profile.json'
    if path.exists():raise RuntimeError('Preserve old results')
    if 'tmp/nufft-openmp' not in finufft.__file__ or 'torch' in sys.modules:raise RuntimeError('Isolated CPU runtime required')
    initial_path=BASE/'continuous-quadrature-optimized/10028-center-0.07-weights.npz'
    saved=np.load(initial_path);g=particle_geometry(ROOT,'10028','inference_half0',radius=5,count=128,seed=609315)
    np.testing.assert_array_equal(saved['indices'],g['indices'])
    checkpoint=BASE/'pose-aware-optimized-shift05/10028-center-0.5-checkpoint.npz'
    weights=np.load(checkpoint)['weights'];op=PolynomialPoseFieldOperator(g['k'],g['q'],g['ctf'],saved['weights'],float(saved['noise_std']),np.deg2rad(.5),.5/g['field_A'],order=32,nthreads=2)
    op.establish_group_scaling();op.set_weights(weights)
    result={'complete':False,'scope':'Design-only eigensolver timing and residuals; Ritz values are not upper bounds',
            'source_snapshot':source_snapshot(ROOT,Path(__file__)),
            'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),'records':[]}
    def save():path.write_text(json.dumps(result,indent=2)+'\n')
    save();v0=np.random.default_rng(609861).normal(size=op.shape[0])
    for ncv in [33,49,17]:
        calls=[0]
        def matvec(v):calls[0]+=1;return op.spatial_gram(v)
        operator=LinearOperator((op.shape[0],)*2,matvec=matvec,dtype=float);start=time.perf_counter();error=None
        try:
            values,vectors=eigsh(operator,k=8,which='LA',v0=v0,tol=1e-4,ncv=ncv,maxiter=15)
            converged=True
        except ArpackNoConvergence as exc:
            values,vectors=exc.eigenvalues,exc.eigenvectors;converged=False;error=str(exc)
        seconds=time.perf_counter()-start;iterations=calls[0]
        residuals=[float(np.linalg.norm(op.spatial_gram(v)-lam*v)) for lam,v in zip(values,vectors.T)] if len(values) else []
        row={'ncv':ncv,'seconds':seconds,'gram_calls':iterations,'converged':converged,'available_modes':len(values),
             'eigenvalues':values.tolist(),'absolute_residual_norms':residuals,'error':error}
        result['records'].append(row);save();print(row,flush=True)
    result['complete']=True;save()


if __name__=='__main__':main()
