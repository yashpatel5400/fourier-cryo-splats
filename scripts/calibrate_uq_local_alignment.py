#!/usr/bin/env python3
"""Pre-interval local-alignment calibration; no coverage is inspected here."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_local_alignment import CachedCellTemplate, refine_local_poses
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/LOCAL-ALIGNMENT-CALIBRATION-PROTOCOL.md'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(ds, snapshot, args):
    out = BASE/args.output/ds
    if out.exists():
        raise RuntimeError(f'Preserve existing calibration: {out}')
    out.mkdir(parents=True)
    begin = time.perf_counter(); path = out/'summary.json'
    result = dict(complete=False, dataset=ds, protocol_sha256=sha(PROTOCOL),
        source_snapshot=snapshot, config=dict(particles=128, radius=5, geometry_seed=609311,
            replicates=128, seed_root=261001, templates=['oracle_reference', 'independent_pilot'],
            initialization_rotation_rms_degrees=2., initialization_shift_rms_A=.5,
            iterations=10, damping=.01, step_rotation_degrees=2., step_shift_A=1.,
            search_rotation_degrees=20., search_shift_A=10.), input_hashes={}, records=[],
        scope='Local refinement calibration on simulation generators; no intervals/coverage computed. Empirical errors are not experimental or deterministic pose radii.')
    def save():
        path.write_text(json.dumps(result, indent=2)+'\n')
    def read(p):
        result['input_hashes'][str(p.relative_to(ROOT))] = sha(p)
        return json.loads(p.read_text())
    def array(p, expected=None):
        digest = sha(p)
        if expected and digest != expected:
            raise ValueError(f'Changed input: {p}')
        result['input_hashes'][str(p.relative_to(ROOT))] = digest
        return np.load(p)
    save()
    try:
        for source in [ROOT/f'data/{ds}/metadata.npz', ROOT/f'data/{ds}/manifest.json', ROOT/f'research/uncertainty/splits/{ds}.csv']:
            result['input_hashes'][str(source.relative_to(ROOT))] = sha(source)
        g = particle_geometry(ROOT, ds, 'inference_half0', radius=5, count=128, seed=609311)
        reg = read(BASE/f'reference-registration-v1/{ds}.json')
        maps = array(ROOT/reg['arrays_file'], reg['arrays_sha256'])
        rho = maps['reference_registered'].ravel().copy(); rho /= np.linalg.norm(rho)
        noise = float(array(BASE/f'continuous-quadrature-optimized/{ds}-center-0.07-weights.npz')['noise_std'])
        checkpoint = array(BASE/f'representation/{ds}/real_particles-spacing-2.0.npz')
        op, pilot, _, _ = model(g, checkpoint, 24, noise=noise); pilot = op.expand(pilot)
        templates = [('oracle_reference', CachedCellTemplate(rho, 64)),
                     ('independent_pilot', CachedCellTemplate(pilot, 24))]
        mean = templates[0][1].signal(g['k'], g['q'], g['ctf'], np.zeros((128, 2)), noise)
        gp = out/'generators.npz'
        np.savez_compressed(gp, truth=rho, pilot=pilot, k=g['k'], q=g['q'], ctf=g['ctf'],
            indices=g['indices'], noise_std=noise, noiseless_signal=mean, field_A=g['field_A'])
        result.update(generators_sha256=sha(gp), field_A=g['field_A'], noise_std_supplied=noise,
            mean_squared_whitened_signal=float(np.mean(mean**2)), setup_seconds=time.perf_counter()-begin)
        save()
        for replicate in range(128):
            rng = np.random.default_rng(np.random.SeedSequence([261001, int(ds), replicate]))
            observed = mean+rng.normal(size=mean.shape)
            rotation0 = Rotation.from_rotvec(rng.normal(size=(128, 3))*np.deg2rad(2.)/np.sqrt(3)).as_matrix()
            shift0 = rng.normal(size=(128, 2))*.5/np.sqrt(2)
            arrays = dict(initial_rotation=rotation0, initial_shift_A=shift0,
                observed=observed)
            row = dict(dataset=ds, replicate=replicate, complete=False, templates=[])
            rp = out/f'replicate-{replicate:03d}.json'
            def save_row():
                rp.write_text(json.dumps(row, indent=2)+'\n')
            save_row()
            try:
                for name, template in templates:
                    tick = time.perf_counter()
                    fit = refine_local_poses(template, g['k'], g['q'], g['ctf'], observed, noise, g['field_A'],
                        initial_rotation=rotation0, initial_shift_A=shift0)
                    rotations = fit.pop('rotations'); shifts = fit.pop('shifts_A')
                    rotvec = Rotation.from_matrix(rotations).as_rotvec()*180/np.pi
                    angular_error = np.linalg.norm(rotvec, axis=1); shift_error = np.linalg.norm(shifts, axis=1)
                    for key, val in [('rotations', rotations), ('shifts_A', shifts), ('rotation_error_vector_degrees', rotvec),
                            ('initial_objective', fit.pop('initial_objective')), ('final_objective', fit.pop('final_objective'))]:
                        arrays[name+'_'+key] = val
                    metrics = dict(template=name, seconds=time.perf_counter()-tick, **fit,
                        rotation_rms_degrees=float(np.sqrt(np.mean(angular_error**2))),
                        shift_rms_A=float(np.sqrt(np.mean(shift_error**2))),
                        rotation_max_degrees=float(angular_error.max()), shift_max_A=float(shift_error.max()),
                        rotation_error_mean_vector_degrees=rotvec.mean(axis=0).tolist(),
                        shift_error_mean_vector_A=shifts.mean(axis=0).tolist(),
                        rotation_quantiles_degrees=np.quantile(angular_error, [0, .5, .9, .99, 1]).tolist(),
                        shift_quantiles_A=np.quantile(shift_error, [0, .5, .9, .99, 1]).tolist())
                    row['templates'].append(metrics); save_row()
                ap = out/f'replicate-{replicate:03d}.npz'; np.savez_compressed(ap, **arrays)
                row.update(complete=True, arrays_sha256=sha(ap)); save_row()
                result['records'].append(dict(replicate=replicate, record_sha256=sha(rp),
                    metrics=[{k: v for k, v in r.items() if k != 'history'} for r in row['templates']]))
                save()
                print(ds, replicate, [(r['template'], r['rotation_rms_degrees'], r['shift_rms_A']) for r in row['templates']], flush=True)
            except Exception as exc:
                row.update(error=repr(exc)); save_row(); raise
        result.update(complete=True, seconds=time.perf_counter()-begin,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
        save()
    except Exception as exc:
        result.update(error=repr(exc), seconds=time.perf_counter()-begin); save(); raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--datasets', default='10028,10049,10076')
    parser.add_argument('--output', default='local-alignment-calibration-v1')
    args = parser.parse_args()
    files = [Path(__file__), PROTOCOL, ROOT/'src/fourier_splats/uq_local_alignment.py',
        ROOT/'tests/test_local_alignment.py']
    for path in files:
        if subprocess.check_output(['git', 'show', f'HEAD:{path.relative_to(ROOT)}'], cwd=ROOT) != path.read_bytes():
            raise ValueError('Commit calibration protocol and sources before outcomes')
    snapshot = source_snapshot(ROOT, Path(__file__), [str(p.relative_to(ROOT)) for p in files[1:]]+
        ['scripts/audit_uq_grid_refinement.py'])
    for ds in args.datasets.split(','):
        if ds not in ['10028', '10049', '10076']:
            raise ValueError('Undeclared dataset')
        run(ds, snapshot, args)


if __name__ == '__main__':
    main()
