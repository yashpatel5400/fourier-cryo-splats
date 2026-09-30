#!/usr/bin/env python3
"""Evaluate a finished declared RELION attempt, retaining its convergence status."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import time
import mrcfile
import numpy as np
import starfile
from fourier_splats.fsc import fsc, resolution
from fourier_splats.physics import fft_volume_center, volume_from_fourier
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_map_alignment import align_map_to_pilot, apply_map_alignment
from fourier_splats.uq_physics import density_functionals
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def model_diagnostics(source):
    """Retain every saved model's local diagnostics, including failed starts."""
    rows = []
    for file in sorted(source.glob('*_model.star')):
        blocks = starfile.read(file, always_dict=True)
        general = blocks.get('model_general', {})
        classes = blocks.get('model_classes')
        fields = ['rlnClassDistribution','rlnAccuracyRotations',
                  'rlnAccuracyTranslationsAngst','rlnEstimatedResolution',
                  'rlnOverallFourierCompleteness']
        rows.append({'file':file.name, 'sha256':sha(file),
            'general':{k:general.get(k) for k in ['rlnCurrentResolution',
                'rlnAveragePmax','rlnSigmaOffsetsAngst','rlnOriginalImageSize','rlnPixelSize']},
            'classes': [] if classes is None else classes.reindex(columns=fields).replace({np.nan:None}).to_dict('records')})
    return {'records':rows,
        'units':{'rlnAccuracyRotations':'degrees','rlnAccuracyTranslationsAngst':'angstrom',
                 'rlnEstimatedResolution':'angstrom','rlnCurrentResolution':'angstrom'},
        'interpretation':'RELION local model diagnostics; not calibrated uniform pose radii or independent reconstruction accuracy. Missing values remain null.'}


def main():
    p = argparse.ArgumentParser(); p.add_argument('dataset', choices=['10028','10049','10076'])
    p.add_argument('--continuation', action='store_true', help='Evaluate the separate v3 CPU continuation')
    args = p.parse_args()
    ds = args.dataset; source = BASE/('relion-reconstruction-v3' if args.continuation else 'relion-reconstruction-v2')/ds
    fitpath = source/'record.json'; fit = json.loads(fitpath.read_text())
    if not fit.get('complete'): raise RuntimeError('Fit still running; do not evaluate checkpoints adaptively')
    out = BASE/('relion-evaluation-v2' if args.continuation else 'relion-evaluation-v1')/ds
    if out.exists(): raise RuntimeError('Preserve previous evaluation, including failures')
    out.mkdir(parents=True); path = out/'metrics.json'
    protocols = ['research/uncertainty/RELION-BASELINE-PROTOCOL.md',
                 'research/uncertainty/RELION-ALIGNMENT-NOTES.md']
    if args.continuation: protocols.append('research/uncertainty/RELION-CONTINUATION-PROTOCOL.md')
    result = {'complete':False,'dataset':ds,'fit_record':str(fitpath.relative_to(ROOT)),
        'fit_record_sha256':sha(fitpath), 'continuation':args.continuation,
        'refinement_converged':fit.get('refinement_converged',False),
        'half_labels_verified':fit.get('half_labels_verified',False),
        'scope':'Development reconstruction comparisons; FSC is not density confidence coverage.',
        'source_snapshot':source_snapshot(ROOT,Path(__file__),protocols)}
    def save(): path.write_text(json.dumps(result,indent=2)+'\n')
    save(); start = time.perf_counter()
    try:
        result['model_diagnostics'] = model_diagnostics(source)
        # Prefer the converged unfiltered maps; otherwise the last paired iterate.
        halves = [source/f'refine_half{h}_class001_unfil.mrc' for h in [1,2]]
        if not all(p.exists() for p in halves):
            candidates = []
            for file in source.glob('refine_it???_half1_class001.mrc'):
                match = re.fullmatch(r'refine_it(\d{3})_half1_class001.mrc',file.name)
                other = file.with_name(file.name.replace('half1','half2'))
                if match and other.exists(): candidates.append((int(match.group(1)),file,other))
            if not candidates:
                result.update(complete=True, skipped=True, seconds=time.perf_counter()-start,
                    reason='No paired refinement maps; initializer is not a converged baseline')
                save(); return
            iteration, *halves = max(candidates)
            result['evaluated_iteration'] = iteration
        else:
            result['evaluated_iteration'] = 'final_unfiltered'
        if not fit.get('half_labels_verified'):
            raise RuntimeError('Exposure half labels were not verified; do not report independent-half FSC')
        files = {}
        def readmap(file):
            files[str(file.relative_to(ROOT))] = sha(file)
            with mrcfile.open(file) as m:
                assert m.data.shape == (64,)*3 and np.isfinite(m.data).all()
                np.testing.assert_allclose(float(m.voxel_size.x),fit['pixel_size_A'],rtol=1e-5)
                return m.data.astype(float)
        native = [readmap(file) for file in halves]; pixel = fit['pixel_size_A']
        pilotpath = BASE/'representation'/ds/'real_particles-spacing-2.0.npz'
        files[str(pilotpath.relative_to(ROOT))] = sha(pilotpath); ck = np.load(pilotpath)
        grid = np.arange(-32,32)/64.; z,y,x = np.meshgrid(grid,grid,grid,indexing='ij')
        xyz = np.column_stack([x.ravel(),y.ravel(),z.ravel()]); mask = np.sum(xyz*xyz,axis=1) <= .35**2
        pilot = np.zeros(64**3)
        # Small blocks avoid a large voxel-by-basis allocation.
        active = np.flatnonzero(mask)
        for first in range(0,len(active),2048):
            ids = active[first:first+2048]
            pilot[ids] = density_functionals(xyz[ids],ck['centers'],float(ck['sigma']))@ck['coefficients']
        pilot = pilot.reshape((64,)*3)
        alignment = align_map_to_pilot((native[0]+native[1])/2,pilot,pixel)
        aligned = [apply_map_alignment(v,alignment['matrix'],alignment['offset']) for v in native]
        result['alignment'] = {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in alignment.items()}
        # Reference enters only after the pilot-selected transform is fixed.
        emd = {'10028':'2660','10049':'6487','10076':'8434'}[ds]
        refpath = ROOT/f'data/uncertainty/references/emd_{emd}.map'
        files[str(refpath.relative_to(ROOT))] = sha(refpath)
        reference = VoxelReference.from_mrc(refpath,box=64)
        np.testing.assert_allclose(reference.field_A,fit['field_A'],rtol=1e-6)
        lockpath = ROOT/'research/uncertainty/confirmation/prediction-v1/locked-models.json'
        lock = json.loads(lockpath.read_text()); files[str(lockpath.relative_to(ROOT))] = sha(lockpath)
        epoch = lock['epochs'][ds]; result['neural_epoch_from_previous_lock'] = epoch
        maps = {'relion':sum(aligned)/2,'reference':reference.volume,'pilot':pilot}
        arrays = {'relion_native_half0':native[0],'relion_native_half1':native[1],
                  'relion_aligned_half0':aligned[0],'relion_aligned_half1':aligned[1]}
        for kind in ['gaussian','voxel','neural']:
            values = []
            for half in [0,1]:
                if kind == 'neural':
                    file = BASE/'neural-reconstruction'/ds/f'half{half}'/f'reconstruct.{epoch}.mrc'
                    value = readmap(file)
                else:
                    file = BASE/'group-reconstruction'/ds/f'{kind}-half{half}.npz'
                    files[str(file.relative_to(ROOT))] = sha(file)
                    value = volume_from_fourier(np.load(file)['fourier'])
                if sha(file) != lock['models'][ds]['files'][str(file.relative_to(ROOT))]:
                    raise RuntimeError('Earlier comparator lock changed')
                values.append(value)
            maps[kind] = sum(values)/2
        result['fsc'] = {}
        pairs = [('relion_half_native',native[0],native[1]),
                 ('relion_half_aligned',aligned[0],aligned[1])]
        pairs += [(f'{kind}_vs_reference',maps[kind],maps['reference']) for kind in ['relion','gaussian','voxel','neural']]
        pairs += [(f'relion_vs_{kind}',maps['relion'],maps[kind]) for kind in ['gaussian','voxel','neural']]
        for name,a,b in pairs:
            curve = fsc(fft_volume_center(a),fft_volume_center(b),pixel)
            np.savetxt(out/f'{name}.csv',curve,delimiter=',',header='shell,frequency_inverse_A,fsc,voxel_count',comments='')
            result['fsc'][name] = {'threshold_0143':resolution(curve,.143),'threshold_05':resolution(curve,.5),
                'mean_fsc':float(np.nanmean(curve[:,2])),'shared_low_frequency_limit_inverse_A':1/40.,
                'interpretation':'Half FSC below 1/40 A is coupled by fitting; aligned/reference comparisons have additional interpolation/frame-selection dependence.'}
        arrays.update(maps); np.savez_compressed(out/'metrics-maps.npz',**arrays)
        result.update(complete=True,seconds=time.perf_counter()-start,input_hashes=files,
            array_sha256=sha(out/'metrics-maps.npz'), alignment_reference_access=False,
            note='Convergence status is inherited and never inferred from a favorable FSC. Deposited maps are approximate references, not ground truth.')
    except Exception as error:
        result.update(complete=True,error=repr(error),seconds=time.perf_counter()-start)
        save(); raise
    save(); print(ds,result['refinement_converged'],result['fsc']['relion_vs_reference'],flush=True)


if __name__ == '__main__': main()
