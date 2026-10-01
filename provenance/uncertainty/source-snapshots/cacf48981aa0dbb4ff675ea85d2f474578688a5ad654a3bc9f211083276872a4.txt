#!/usr/bin/env python3
"""Verify and summarize complete frozen studies; never omit failed trials."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_end_to_end import binomial_interval
from fourier_splats.uq_end_to_end_summary import aggregate_replicates, distribution, TEMPLATES
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(ds, source, output, snapshot):
    directory = BASE/source/ds; sp = directory/'summary.json'
    study = json.loads(sp.read_text())
    if not study.get('complete') or len(study['records']) != 200:
        raise ValueError(f'{ds}: all 200 prescribed trials must be attempted before final summary')
    hashes = {str(sp.relative_to(ROOT)): sha(sp)}
    records, arrays = [], {}
    for entry in study['records']:
        rep = entry['replicate']; rp = directory/f'replicate-{rep:03d}.json'
        if sha(rp) != entry['record_sha256']:
            raise ValueError(f'{rp}: record changed')
        record = json.loads(rp.read_text()); records.append(record)
        hashes[str(rp.relative_to(ROOT))] = sha(rp)
        suffix = '' if record['complete'] else '-partial'
        ap = directory/f'replicate-{rep:03d}{suffix}.npz'
        expected = record['arrays_sha256' if record['complete'] else 'partial_arrays_sha256']
        if sha(ap) != expected:
            raise ValueError(f'{ap}: arrays changed')
        hashes[str(ap.relative_to(ROOT))] = expected
        if record['complete']:
            with np.load(ap) as data:
                arrays[rep] = {t: np.concatenate([
                    Rotation.from_matrix(data[t+'_rotations']).as_rotvec()*180/np.pi,
                    data[t+'_shifts_A']], axis=1) for t in TEMPLATES}
    # Recheck independent calibration inputs against the launch snapshot.
    cal = BASE/'local-alignment-calibration-v1'/ds
    for name, expected in study['calibration_hashes'].items():
        p = cal/name
        if sha(p) != expected:
            raise ValueError(f'{p}: calibration changed')
        hashes[str(p.relative_to(ROOT))] = expected
    result = aggregate_replicates(records, 200)
    result.update(dataset=ds, source_snapshot=snapshot, input_hashes=hashes,
        protocol_sha256=study['protocol_sha256'], config=study['config'],
        calibration_parameters=study['calibration_parameters'], seconds=study['seconds'],
        peak_resident_bytes=study['peak_resident_bytes'], alignment=[], numerical=[],
        errors=[dict(replicate=r['replicate'], error=r.get('error')) for r in records if not r['complete']])
    for template in TEMPLATES:
        poses = [next(t for t in r['templates'] if t['template'] == template) for r in records if r['complete']]
        error = np.stack([arrays[r['replicate']][template] for r in records if r['complete']]) if poses else None
        count = sum(p['pose_inside_calibrated_ball'] for p in poses)
        diagnostics = dict(template=template, available=len(poses), planned=200,
            pose_bound_covered=count, pose_bound_fraction=count/200,
            pose_bound_exact_binomial_95=binomial_interval(count, 200),
            rotation_rms_degrees=distribution([p['rotation_rms_degrees'] for p in poses]),
            shift_rms_A=distribution([p['shift_rms_A'] for p in poses]),
            joint_score=distribution([p['observed_joint_score'] for p in poses]),
            rotation_boundary_particles=distribution([p['alignment']['rotation_boundary_particles'] for p in poses]),
            shift_boundary_particles=distribution([p['alignment']['shift_boundary_particles'] for p in poses]))
        if error is not None:
            means = error.mean(axis=0); centered = error-means[None]
            # These are unconditional across-replicate diagnostics, not proof of
            # centering conditional on estimated poses, fitted weights or design.
            correlations = []
            for coordinate in range(5):
                e = centered[:, :, coordinate]
                norms = np.linalg.norm(e, axis=0)
                valid = norms > 0
                normalized = e[:, valid]/norms[valid]
                corr = normalized.T@normalized
                correlations.extend(corr[np.triu_indices(len(corr), 1)].tolist())
            diagnostics.update(coordinate_mean=error.mean(axis=(0, 1)).tolist(),
                coordinate_rms=np.sqrt(np.mean(error**2, axis=(0, 1))).tolist(),
                per_particle_coordinate_means=means.tolist(),
                empirical_cross_particle_coordinate_correlations=distribution(correlations),
                centering_scope='Unconditional simulation diagnostics only; not conditional-on-design centering or an independence test.')
        result['alignment'].append(diagnostics)
        for target in ['center', 'contrast']:
            fits = [next(f for f in p['fits'] if f['target'] == target) for p in poses]
            fixed = [f['fixed_audit'] for f in fits]
            result['numerical'].append(dict(template=template, target=target, method='fixed_audit',
                solves=len(fixed), converged=sum(f['converged'] for f in fixed),
                relative_gap=distribution([(f['objective']-f['dual_lower_bound'])/f['objective'] for f in fixed])))
            for name in ['gaussian_tau2_fixed', 'gaussian_tau2_pose', 'gaussian_tau1_fixed', 'gaussian_tau1_pose']:
                values = [next(g for g in f['gaussian_fits'] if g['method'] == name) for f in fits]
                result['numerical'].append(dict(template=template, target=target, method=name,
                    solves=len(values), converged=sum(v['converged'] for v in values),
                    cg_relative_residual=distribution([v['cg_relative_residual'] for v in values]),
                    cg_iterations=distribution([v['cg_iterations'] for v in values]),
                    variance_identity_relative_error=distribution([v['variance_identity_relative_error'] for v in values])))
    output.mkdir(parents=True, exist_ok=True)
    path = output/f'{ds}.json'
    if path.exists():
        raise RuntimeError('Preserve existing final summary')
    path.write_text(json.dumps(result, indent=2)+'\n')
    csv_rows = []
    for row in result['groups']:
        flat = {k: row[k] for k in ['template', 'target', 'method', 'image_mode', 'available_intervals', 'missing_or_failed_intervals']}
        for name in ['covered', 'raw_covered', 'correct_sign_exclusion', 'raw_correct_sign_exclusion', 'fallback']:
            flat[name+'_count'] = row[name]['successes']
            flat[name+'_rate'] = row[name]['fraction']
            flat[name+'_lower'], flat[name+'_upper'] = row[name]['exact_binomial_95']
        for name in ['half_width', 'raw_half_width', 'relative_half_width', 'raw_relative_half_width', 'half_width_over_abs_center']:
            flat[name+'_median'] = row[name].get('median')
        flat['relative_width_below_half_count'] = row['relative_width_below_half_count']
        csv_rows.append(flat)
    with (output/f'{ds}.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0])); writer.writeheader(); writer.writerows(csv_rows)
    print(ds, result['completed_replicates'], 'completed;', result['failed_replicates'], 'failed;', len(hashes), 'verified inputs', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', default='10028,10049,10076')
    parser.add_argument('--source', default='end-to-end-local-pose-v1')
    parser.add_argument('--output', default='end-to-end-local-pose-summary-v1')
    args = parser.parse_args()
    snapshot = source_snapshot(ROOT, Path(__file__), ['src/fourier_splats/uq_end_to_end_summary.py',
        'tests/test_end_to_end_summary.py', 'research/uncertainty/END-TO-END-LOCAL-POSE-PROTOCOL.md'])
    for ds in args.datasets.split(','):
        summarize(ds, args.source, BASE/args.output, snapshot)


if __name__ == '__main__':
    main()
