#!/usr/bin/env python3
"""Lock pilot-only target coordinates before new reference-outcome evaluations."""
import hashlib
import json
from pathlib import Path
import numpy as np
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_continuous import cell_target_coefficients
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]; BASE = ROOT/'results/uncertainty/development'


def main():
    output = ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'
    if output.exists():
        raise RuntimeError('Do not overwrite selected targets')
    box = 24; width_A = 20.; separation_A = 40.; selection_radius = .30
    records = []
    for dataset in ['10028', '10049', '10076']:
        source = BASE/'continuous-quadrature-optimized'/f'{dataset}.json'; prior = json.loads(source.read_text())
        g = particle_geometry(ROOT, dataset, 'inference_half0', radius=5, count=128, seed=prior['config']['seed'])
        checkpoint = BASE/'representation'/dataset/'real_particles-spacing-2.0.npz'
        op, pilot, _, _ = model(g, np.load(checkpoint), box); pilot = op.expand(pilot)
        np.testing.assert_allclose(np.linalg.norm(pilot), 1., rtol=1e-12)
        coordinates = (np.indices((box,)*3).reshape(3, -1).T[:, ::-1]+.5)/box-.5
        coordinates = coordinates[np.linalg.norm(coordinates, axis=1) <= selection_radius]
        width = width_A/g['field_A']
        candidates = [(float(cell_target_coefficients(box, [x], [1], width)@pilot), x) for x in coordinates]
        candidates.sort(key=lambda pair: (-pair[0], *pair[1].tolist()))
        selected = []
        for value, x in candidates:
            if value <= 0:
                break
            if all(np.linalg.norm(x-np.array(r['center_fraction_field']))*g['field_A'] >= separation_A for r in selected):
                selected.append({'name': f'pilot_region_{len(selected)+1}', 'center_fraction_field': x.tolist(),
                                 'center_A': (x*g['field_A']).tolist(), 'pilot_expected_feature': value})
                if len(selected) == 3:
                    break
        control = {'name': 'matched_center', 'center_fraction_field': [0., 0., 0.], 'center_A': [0., 0., 0.],
                   'pilot_expected_feature': float(cell_target_coefficients(box, [[0, 0, 0]], [1], width)@pilot)}
        records.append({'dataset': dataset, 'field_A': g['field_A'], 'width_A': width_A, 'width_fraction_field': width,
                        'candidate_count': len(candidates), 'selected_count': len(selected), 'features': [*selected, control],
                        'pilot_checkpoint_sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                        'source_geometry_config': prior['config'], 'source_record_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                        'planned_radius12_band_endpoint_A': g['field_A']/12})
    result = {'status': 'Target coordinates locked before new reference-map outcomes; full experiment not yet launched',
              'scope': 'Exploratory development target-selection lock, not a new frozen confirmation study',
              'selection': {'width_A': width_A, 'separation_A': separation_A, 'candidate_box': box,
                            'maximum_center_radius_fraction_field': selection_radius,
                            'rule': 'Greedy descending positive pilot score, reject centers within 40 A of a previous choice; lexicographic x,y,z ties'},
              'inputs': 'Independent Gaussian pilot coefficients and geometry metadata only; no deposited map or new inference-image outcome is read',
              'source_snapshot': source_snapshot(ROOT, Path(__file__), ['scripts/audit_uq_grid_refinement.py', 'research/uncertainty/PILOT-SELECTED-TARGETS-PROPOSAL.md']),
              'datasets': records}
    output.parent.mkdir(parents=True, exist_ok=True); output.write_text(json.dumps(result, indent=2)+'\n')
    print('LOCKED', hashlib.sha256(output.read_bytes()).hexdigest())
    for row in records:
        print(row['dataset'], [(f['name'], f['center_A'], f['pilot_expected_feature']) for f in row['features']])


if __name__ == '__main__':
    main()
