#!/usr/bin/env python3
"""Gaussian two-point lower bounds for every completed continuous development fit."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
OUT = BASE/'continuous-fixed-length-lower'


def main():
    OUT.mkdir(exist_ok=True)
    snapshot = source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py'])
    records = []; sources = {}
    for radius in [5, 12]:
        for dataset in ['10028', '10049', '10076']:
            folder = 'continuous-quadrature-optimized' if radius == 5 else (
                'continuous-quadrature-high-band-probe' if dataset == '10028' else 'continuous-quadrature-high-band-additional')
            path = BASE/folder/f'{dataset}.json'; data = json.loads(path.read_text()); cfg = data['config']
            if len(data['targets']) != (4 if radius == 5 else 2):
                raise AssertionError('Incomplete fit group')
            sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
            g = particle_geometry(ROOT, dataset, 'inference_half0', radius=radius, count=cfg['particles'], seed=cfg['seed'])
            ck = np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
            _, _, _, noise = model(g, ck, 24)
            gram = QuadratureObservationGram(g['k'], g['ctf'], noise, order=cfg['quadrature_order'], preconditioner_rank=0)
            for row in data['targets']:
                name = row['target']; width = row['width_fraction_field']
                saved = np.load(BASE/folder/f'{dataset}-{name}-{width}-weights.npz')
                np.testing.assert_array_equal(saved['indices'], g['indices']); w = saved['weights']
                centers = [[0, 0, 0]] if name == 'center' else [[0, 0, .08], [0, 0, -.08]]
                signs = [1] if name == 'center' else [1, -1]
                a, ell2 = gram.target(centers, signs, width); gw = gram.matvec(w)
                error = gram.quadrature_error(w)
                h2 = ell2-2*w@a+w@gw+error['squared_field_norm']
                # Same numerical padding convention as the solver. Quadrature
                # has an analytic bound; roundoff is not validated arithmetic.
                pad = 20*np.finfo(float).eps*len(w)*(ell2+2*abs(w@a)+abs(w@gw))
                hnorm = np.sqrt(max(0., h2)+pad)
                ahnorm = np.linalg.norm(a-gw)+error['gram_action_norm']
                pair = float(ell2-w@a)
                directions = []
                for label, length, observed_norm, pairing in [
                    ('adjoint_residual', hnorm, ahnorm, pair),
                    ('target', np.sqrt(ell2), np.linalg.norm(a), ell2)]:
                    scale = min(2/length, norm.ppf(.95)/observed_norm if observed_norm > 0 else np.inf)
                    directions.append({'direction': label, 'norm_upper_numerical': float(length),
                                       'observation_norm_upper_numerical': float(observed_norm),
                                       'target_pairing_numerical': float(pairing), 'scale': float(scale),
                                       'half_length_lower_numerical': float(scale*abs(pairing))})
                lower = max(d['half_length_lower_numerical'] for d in directions)
                half = row['fit']['half_width']
                if lower > half+1e-6:
                    raise AssertionError('Lower bound exceeds valid fixed-length upper bound')
                records.append({'dataset': dataset, 'frequency_radius': radius, 'target': name, 'width': width,
                                'directions': directions, 'reported_half_width': half,
                                'fixed_length_half_width_lower_numerical': lower,
                                'lower_over_upper': lower/half, 'quadrature_errors': error})
                print(dataset, radius, name, width, 'lower / upper', lower/half, flush=True)
    summary = {'stage': 'post-hoc design-only classical Gaussian testing lower bounds',
               'source_snapshot': snapshot, 'sources': sources, 'records': records,
               'settings': len(records), 'lower_over_upper_range': [min(r['lower_over_upper'] for r in records), max(r['lower_over_upper'] for r in records)],
               'scope': 'All uniformly covering fixed-length intervals under the fixed-pose continuous L2 class; not random-length pointwise bounds.',
               'numerics': 'Analytic quadrature remainders; ordinary floating point, not validated arithmetic.'}
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
