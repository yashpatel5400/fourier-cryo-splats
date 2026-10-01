#!/usr/bin/env python3
"""Review-2 frame correction for every original finite-dictionary example."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import VoxelObservationOperator, VoxelReference, particle_geometry, local_weights, resample_volume
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_subspace import orthonormalize_columns
from fourier_splats.uncertainty import optimize_certificate, bias_aware_half_width
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_representation import dictionary
from audit_uq_ambient import SupportedOperator

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/REGISTERED-DICTIONARY-PROTOCOL.md'
MAPS = {'10028': '2660', '10049': '6487', '10076': '8434'}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    for p in [Path(__file__), PROTOCOL]:
        if subprocess.check_output(['git', 'show', f'HEAD:{p.relative_to(ROOT)}'], cwd=ROOT) != p.read_bytes():
            raise ValueError('Commit protocol/source before new replay')
    out = BASE/'registered-dictionary-replay-v1'
    if out.exists():
        raise RuntimeError('Preserve previous replay')
    out.mkdir()
    old_path = BASE/'ambient-full-solve.json'; old = json.loads(old_path.read_text())
    snapshot = source_snapshot(ROOT, Path(__file__), [str(PROTOCOL.relative_to(ROOT)),
        'scripts/audit_uq_representation.py', 'scripts/audit_uq_ambient.py'])
    for case in old['cases']:
        ds = case['dataset']; cfg = case['config']; box = cfg['box']; start = time.perf_counter()
        path = out/f'{ds}.json'
        result = dict(complete=False, dataset=ds, source_snapshot=snapshot, historical_config=cfg,
            old_record_sha256=sha(old_path), input_hashes={}, targets=[],
            scope='Reconstructed historical fixed-pose finite-grid dictionary estimator; registered generator sensitivity, not experimental coverage.')
        def save():
            path.write_text(json.dumps(result, indent=2)+'\n')
        def record(p):
            result['input_hashes'][str(p.relative_to(ROOT))] = sha(p)
        save()
        try:
            g = particle_geometry(ROOT, ds, 'inference_half0', radius=cfg['frequency_radius'], count=cfg['particles'], seed=cfg['seed'])
            for p in [ROOT/f'data/{ds}/metadata.npz', ROOT/f'data/{ds}/manifest.json', ROOT/f'research/uncertainty/splits/{ds}.csv']:
                record(p)
            cp = BASE/f'representation/{ds}/real_particles-spacing-2.0.npz'; record(cp)
            if sha(cp) != case['pilot_sha256']:
                raise ValueError('Historical pilot changed')
            ck = np.load(cp); coord = np.arange(-box//2, box//2)/box
            z, y, x = np.meshgrid(coord, coord, coord, indexing='ij')
            xyz = np.column_stack([x.ravel(), y.ravel(), z.ravel()]); mask = (np.sum(xyz**2, axis=1) <= cfg['support']**2).astype(float)
            pilot = density_functionals(xyz, ck['centers'], float(ck['sigma']))@ck['coefficients']*mask
            pilot /= np.linalg.norm(pilot)
            signal = VoxelObservationOperator(g['k'], g['ctf'], box).forward(pilot)
            noise = float(np.sqrt(np.mean(signal**2)/cfg['snr']))
            op = SupportedOperator(VoxelObservationOperator(g['k'], g['ctf'], box, noise), mask)
            centers = dictionary(cfg['dictionary_radius'], cfg['spacing'])
            basis = orthonormalize_columns(density_functionals(xyz, centers, cfg['sigma'])*mask[:, None])
            a = op.forward_columns(basis)
            rp = ROOT/f'data/uncertainty/references/emd_{MAPS[ds]}.map'; record(rp)
            original_volume = VoxelReference.from_mrc(rp, box=box).volume
            reg_path = BASE/f'reference-registration-v1/{ds}.json'; record(reg_path)
            reg = json.loads(reg_path.read_text()); ap = ROOT/reg['arrays_file']; record(ap)
            if sha(ap) != reg['arrays_sha256']:
                raise ValueError('Frozen registration changed')
            saved = np.load(ap)
            registered_volume = resample_volume(saved['reference_registered'], box)
            via64 = resample_volume(saved['reference_original'], box)
            references = {}
            for frame, volume in [('original', original_volume), ('registered', registered_volume)]:
                rho = volume.ravel()*mask; references[frame] = rho/np.linalg.norm(rho)
            arrays = dict(pilot=pilot, reference_original=references['original'], reference_registered=references['registered'], indices=g['indices'])
            result.update(noise_std=noise, historical_noise_relative_error=abs(noise-case['noise_std'])/case['noise_std'],
                original_two_stage_resample_relative_error=float(np.linalg.norm(via64-original_volume)/np.linalg.norm(original_volume)))
            for target in case['targets']:
                name, width = target['target'], target['width_fraction_field']
                locations, signs = ([[0., 0., 0.]], [1.]) if name == 'center' else ([[0., 0., .08], [0., 0., -.08]], [1., -1.])
                ell = sum(s*local_weights(box, np.array(c)*box, width*box).ravel() for c, s in zip(locations, signs))*mask
                fit = optimize_certificate(a, np.zeros((1, len(a), 0)), basis.T@ell, cfg['radius'], maxiter=300, rtol=.0005)
                w = fit.weights; full_bias = cfg['radius']*np.linalg.norm(ell-op.adjoint(w))
                ambient_half = bias_aware_half_width(fit.noise_sd, full_bias); no_data = cfg['radius']*np.linalg.norm(ell)
                array_key = f'{name}-{width:g}-weights'; arrays[array_key] = w
                row = dict(target=name, width_fraction_field=width, array_key=array_key,
                    restricted_half_width=float(fit.half_width), audited_half_width=float(ambient_half),
                    no_data_half_width=float(no_data), noise_sd=float(fit.noise_sd), checks=[])
                errors = [abs(fit.half_width/(target['restricted_width_fraction']*no_data)-1),
                          abs(ambient_half/(target['audited_width_fraction']*no_data)-1)]
                for frame, rho in references.items():
                    delta = rho-pilot; bias = float(w@op.forward(delta)-ell@delta); truth = float(ell@rho)
                    for method, half in [('restricted_dictionary', fit.half_width), ('audited_dictionary', ambient_half)]:
                        coverage = float(norm.cdf((half-bias)/fit.noise_sd)-norm.cdf((-half-bias)/fit.noise_sd))
                        check = dict(frame=frame, method=method, true_target=truth, expected_center=truth+bias,
                            bias=bias, analytic_coverage=coverage, reference_pilot_distance=float(np.linalg.norm(delta)))
                        if frame == 'original':
                            previous = next(r for r in target['coverage'] if r['truth'] == 'independent_EMDB_voxel_reference' and r['method'] == method)
                            check.update(historical_bias_difference=bias-previous['bias'], historical_coverage_difference=coverage-previous['analytic_coverage'])
                            errors.extend([abs(check['historical_bias_difference']), abs(check['historical_coverage_difference'])])
                        row['checks'].append(check)
                row.update(historical_replay_max_error=float(max(errors)), historical_replay_passed=bool(max(errors) <= 1e-6))
                result['targets'].append(row); save()
            arrays_path = out/f'{ds}-arrays.npz'; np.savez_compressed(arrays_path, **arrays)
            result.update(complete=True, arrays_sha256=sha(arrays_path), seconds=time.perf_counter()-start,
                all_historical_replays_passed=all(t['historical_replay_passed'] for t in result['targets']))
            save(); print(ds, result['all_historical_replays_passed'], result['seconds'], flush=True)
        except Exception as exc:
            result.update(error=repr(exc), seconds=time.perf_counter()-start); save(); raise


if __name__ == '__main__':
    main()
