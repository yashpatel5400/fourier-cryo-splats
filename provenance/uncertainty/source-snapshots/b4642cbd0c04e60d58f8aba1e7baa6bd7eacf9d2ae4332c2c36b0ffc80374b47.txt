#!/usr/bin/env python3
"""Archive-scale direct-sum checks for completed known-pilot pose audits."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_cell_moments import check_cell_pose_pairings
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--folder',default='joint-bias-pilot-sharp-audit')
    parser.add_argument('--output',default='pilot-pairing-direct-checks.json')
    args=parser.parse_args();out=BASE/'audit-regressions'/args.output
    if out.exists():raise RuntimeError('Preserve previous check')
    result={'complete':False,'stage':'Selected-column numerical regression checks; not an all-column floating-point certificate',
            'config':vars(args),'source_snapshot':source_snapshot(ROOT,Path(__file__),['scripts/audit_uq_grid_refinement.py']),
            'checks':[]};start=time.perf_counter()
    def save():out.write_text(json.dumps(result,indent=2)+'\n')
    save()
    try:
        for audit in sorted((BASE/args.folder).glob('*/*.json')):
            evidence=json.loads(audit.read_text())
            if not evidence.get('complete') or evidence.get('error'):raise ValueError('All inventoried audits must be complete')
            source=ROOT/evidence['config']['fit'];fit=json.loads(source.read_text())
            assert hashlib.sha256(source.read_bytes()).hexdigest()==evidence['source_fit_sha256']
            saved=np.load(source.with_name(source.stem+'-weights.npz'));noise=float(saved['noise_std'])
            geometry=particle_geometry(ROOT,fit['dataset'],'inference_half0',radius=5,count=128,seed=fit['source_geometry_config']['seed'])
            np.testing.assert_array_equal(saved['indices'],geometry['indices'])
            op=PolynomialPoseFieldOperator(geometry['k'],geometry['q'],geometry['ctf'],saved['weights'],noise,
                    np.deg2rad(fit['rotation_radius_degrees']),fit['translation_radius_A']/geometry['field_A'],order=2,
                    column_denominators=evidence['pose_scaling']['column_denominators'])
            checkpoint_path=BASE/'representation'/fit['dataset']/'real_particles-spacing-2.0.npz'
            assert hashlib.sha256(checkpoint_path.read_bytes()).hexdigest()==evidence['pilot_refinement']['pilot_checkpoint_sha256']
            cellop,pilot,_,check_noise=model(geometry,np.load(checkpoint_path),24);pilot=cellop.expand(pilot)
            np.testing.assert_allclose(noise,check_noise)
            pair_path=audit.with_suffix('.npz');pair=np.load(pair_path)['pilot_pose_cross']
            check=check_cell_pose_pairings(op,pilot,24,pair)
            result['checks'].append({'audit':str(audit.relative_to(ROOT)),
                'audit_sha256':hashlib.sha256(audit.read_bytes()).hexdigest(),
                'pairing_array_sha256':hashlib.sha256(pair_path.read_bytes()).hexdigest(),**check})
            save();print(audit.parent.name,audit.stem,check['relative_pairing_difference'],flush=True)
        if not result['checks'] or not all(row['passed'] for row in result['checks']):raise AssertionError('Direct-sum regression failed')
        result.update(complete=True,seconds=time.perf_counter()-start);save()
    except Exception as exc:
        result.update(complete=False,error=repr(exc),seconds=time.perf_counter()-start);save();raise


if __name__=='__main__':main()
