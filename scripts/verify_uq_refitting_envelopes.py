#!/usr/bin/env python3
"""Independent sinc check and higher-order quadrature of post hoc pose diagnostics."""
from pathlib import Path
import csv, json, hashlib
import numpy as np
from fourier_splats.uq_continuous import fourier_field_norm2, gaussian_target_integrals
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_refitting_diagnostics import dephase_adjoint_weights
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
rows=[]
for ds in ['10028','10049','10076']:
 with np.load(BASE/'local-alignment-calibration-v1'/ds/'generators.npz') as f:g={k:f[k].copy() for k in f.files}
 original=BASE/'end-to-end-local-pose-v1'/ds
 saved=list(csv.DictReader((BASE/'refitting-bias-reanalysis-v1'/ds/'realized-envelopes.csv').open()))
 for rep in [0,99,199]:
  with np.load(original/f'replicate-{rep:03d}.npz') as f:a={k:f[k].copy() for k in f.files}
  for template in ['oracle_reference','independent_pilot']:
   k=g['k'];kh=k@a[template+'_rotations'];noise=float(g['noise_std'])
   gt=QuadratureObservationGram(k,g['ctf'],noise,order=48,eps=1e-13,preconditioner_rank=0)
   gf=QuadratureObservationGram(kh,g['ctf'],noise,order=48,eps=1e-13,preconditioner_rank=0)
   for target,centers,signs in [('center',[[0.,0.,0.]],[1.]),('contrast',[[0.,0.,.08],[0.,0.,-.08]],[1.,-1.])]:
    w=a[template+'_'+target+'_audit_weights'];v=dephase_adjoint_weights(w,g['q'],a[template+'_shifts_A'],float(g['field_A']))
    ct=gt.coefficients(v);cf=gf.coefficients(w);ft=gt.field(v);ff=gf.field(w)
    at,ell2=gt.target(centers,signs,.07);res2=ell2-2*v@at+gt.weights@(ft*ft);pose2=gt.weights@((ft-ff)**2)
    row=next(r for r in saved if int(r['replicate'])==rep and r['template']==template and r['target']==target)
    result=dict(dataset=ds,replicate=rep,template=template,target=target,
      residual_order40_48_difference=abs(float(row['residual_norm2_unpadded'])-res2),pose_order40_48_difference=abs(float(row['pose_norm2_unpadded'])-pose2))
    if rep==0 and template=='independent_pilot' and target=='center':
     norm2,_=fourier_field_norm2(k,ct)
     dense_res2=ell2-2*v@at+norm2
     dense_pose2,_=fourier_field_norm2(np.concatenate([k.reshape(-1,3),kh.reshape(-1,3)]),np.concatenate([ct.ravel(),-cf.ravel()]))
     result.update(residual_dense_sinc_difference=abs(dense_res2-res2),pose_dense_sinc_difference=abs(dense_pose2-pose2))
    if max(value for key,value in result.items() if key.endswith('difference'))>1e-7:raise ArithmeticError(result)
    rows.append(result)
 print(ds,'verified',flush=True)
out=ROOT/'provenance/uncertainty/refitting-realized-envelope-verification.json'
assert not out.exists()
out.write_text(json.dumps(dict(complete=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  scope='36 preset higher-order/tighter-tolerance checks, including three independent full sinc integrations. Numerical verification, not interval arithmetic.',records=rows),indent=2)+'\n')
