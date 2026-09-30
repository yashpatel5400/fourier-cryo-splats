#!/usr/bin/env python3
"""Audit the actual Fourier-baseline weights in the continuous density class.

Records a known pilot/interpolation offset rather than silently changing the
baseline center. Fixed-pose only; nonlinear pose audit is a separate extension.
"""
import json
import hashlib
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_fourier_variational import HermitianTrilinearOperator
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    out = BASE/'fourier-variational-continuous-audit'; out.mkdir(exist_ok=True)
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    for dataset in ['10028', '10049', '10076']:
        path = out/f'{dataset}.json'
        if path.exists(): raise RuntimeError('Preserve previous audits')
        start = time.perf_counter(); source = BASE/'fourier-variational-baseline'/f'{dataset}.json'
        baseline = json.loads(source.read_text())
        if not baseline['complete']: raise RuntimeError('Baseline is incomplete')
        cfg = baseline['config']; noise = baseline['noise_std_supplied']
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=cfg['radius'], count=cfg['particles'], seed=609841)
        gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=80, preconditioner_rank=0)
        op = HermitianTrilinearOperator(g['k'], g['ctf'], noise, box=cfg['box'])
        ck = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
        old, pilot, _, _ = model(g, ck, 24, noise=noise); pilot = old.expand(pilot)
        pilot_x = op.project_cells(pilot, 24)
        signal_offset = cell_forward(g['k'], g['ctf'], pilot, 24, noise)-op.matrix@pilot_x
        record = {'stage': 'Fixed-pose continuous re-audit of Fourier Gaussian estimators; supplied simulation noise and B=2',
            'dataset': dataset, 'source_snapshot': snapshot, 'baseline_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'complete': False, 'records': []}
        def save(): path.write_text(json.dumps(record, indent=2)+'\n')
        save()
        try:
            for row in baseline['records']:
                target = row['target']; prior = np.sqrt(row['prior_expected_squared_deviation_norm'])
                saved_path = source.parent/f'{dataset}-{target}-prior{prior:g}.npz'
                saved = np.load(saved_path); w = saved['weights']; np.testing.assert_array_equal(saved['indices'], g['indices'])
                width = row['sigma_A']/g['field_A']
                centers = [[0, 0, 0]] if target == 'center' else [[0, 0, width], [0, 0, -width]]
                signs = [1] if target == 'center' else [1, -1]
                a, norm2 = gram.target(centers, signs, width); gw = gram.matvec(w)
                errors = gram.quadrature_error(w)
                pad = 50*np.finfo(float).eps*len(w)*(norm2+2*abs(w@a)+abs(w@gw))
                residual = np.sqrt(max(0., norm2-2*w@a+w@gw+errors['squared_field_norm'])+pad)
                pilot_target = float(cell_target_coefficients(24, centers, signs, width)@pilot)
                offset = float(saved['target_coefficients']@pilot_x-pilot_target+w@signal_offset)
                sd = float(np.linalg.norm(w)); bound = 2*residual+abs(offset)
                half = bias_aware_half_width_stable(sd, bound)
                recentered_half = bias_aware_half_width_stable(sd, 2*residual)
                no_data = 2*np.sqrt(norm2)
                original = next(c for c in row['checks'] if c['method'] == 'full_fourier_gaussian' and c['scenario'] == 'continuous_reference_nominal')
                truth = original['true_target']; mean = original['expected_center']; bias = mean-truth
                result = {'target': target, 'prior_expected_squared_deviation_norm': prior**2,
                    'weights_sha256': hashlib.sha256(saved_path.read_bytes()).hexdigest(),
                    'continuous_residual_norm': float(residual), 'known_pilot_center_offset': offset,
                    'same_estimator_bias_upper': float(bound), 'same_estimator_half_width': half,
                    'recentered_half_width': recentered_half, 'no_data_half_width': float(no_data),
                    'relative_half_width_before_fallback': float(half/no_data),
                    'half_width_over_full_posterior': float(half/original['half_width']),
                    'reference_target': truth, 'reference_expected_center': mean,
                    'reference_same_estimator_coverage': float(norm.cdf((half-bias)/sd)-norm.cdf((-half-bias)/sd)),
                    'reference_same_estimator_correct_sign_probability': float(norm.cdf((np.sign(truth)*mean-half)/sd)),
                    'reference_recentered_correct_sign_probability': float(norm.cdf((np.sign(truth)*(mean-offset)-recentered_half)/sd)),
                    'quadrature_errors': errors, 'squared_roundoff_pad': float(pad)}
                record['records'].append(result); save()
                print(dataset, target, prior, 'width/no-data', result['relative_half_width_before_fallback'],
                      'power', result['reference_same_estimator_correct_sign_probability'], flush=True)
            record.update(complete=True, seconds=time.perf_counter()-start); save()
        except Exception as exc:
            record.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__': main()
