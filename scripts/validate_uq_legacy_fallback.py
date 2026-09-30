#!/usr/bin/env python3
"""Archive-dependent regression: reconstruct raw centers from legacy metadata.

The copied metadata deliberately claims a prior fallback and omits raw means.
No scientific fit, image, weight, pose, or stored diagnostic is modified.
"""
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_joint_bias import recover_reference_centers

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    source=BASE/'pose-exchange-conic-duals/10049-center-2.json'
    diagnostic=BASE/'pose-optimized-diagnostics/pose-exchange-conic-duals/10049-center-2.json'
    output=BASE/'audit-regressions/legacy-fallback-reconstruction.json';output.parent.mkdir(exist_ok=True)
    if output.exists():raise RuntimeError('Preserve regression evidence')
    prior=json.loads(source.read_text());original=json.loads(diagnostic.read_text());previous=copy.deepcopy(original)
    previous['uses_no_data']=True;previous.pop('pilot_target',None)
    for row in previous['reference_checks']:row.pop('raw_expected_center',None)
    g=particle_geometry(ROOT,'10049','inference_half0',radius=5,count=128,seed=prior['source_geometry_config']['seed'])
    wp=source.with_name(source.stem+'-weights.npz');saved=np.load(wp)
    angle=np.full(128,np.deg2rad(2));shift=np.full(128,.5/g['field_A'])
    result={'stage':'Forced legacy-metadata regression, not an additional scientific experiment',
        'complete':False,'source_snapshot':source_snapshot(ROOT,Path(__file__),['scripts/audit_uq_joint_bias.py','scripts/audit_uq_grid_refinement.py']),
        'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,diagnostic,wp]}}
    def save():output.write_text(json.dumps(result,indent=2)+'\n')
    save()
    try:
        pilot,means=recover_reference_centers(previous,g,saved['weights'],float(saved['noise_std']),'10049','center',angle,shift)
        old=np.array([r['expected_center'] for r in original['reference_checks']])
        new=np.array([means[r['scenario']] for r in original['reference_checks']])
        np.testing.assert_allclose(new,old,rtol=1e-10,atol=1e-10)
        result.update(complete=True,pilot_target=pilot,recovered_raw_means=means,maximum_absolute_difference=float(np.max(abs(new-old))))
        save();print('PASS',result['maximum_absolute_difference'])
    except Exception as exc:
        result.update(complete=False,error=repr(exc));save();raise


if __name__=='__main__':main()
