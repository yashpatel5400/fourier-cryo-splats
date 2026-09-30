#!/usr/bin/env python3
"""Original CryoLike comparison under the declared exploratory CPU protocol."""
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import resource
import subprocess
import time
import mrcfile
import numpy as np
import torch
from cryolike.metadata import ViewingAngles
from fourier_splats.uq_cryolike_adapter import make_grid,prepare_images,prepare_templates,prepare_ctf,score
from fourier_splats.physics import volume_from_fourier
from fourier_splats.uq_data import VoxelReference
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
AUTHOR=ROOT/'tmp/cryolike-source'
COMMIT='a413ffd265e2815c9f81c492a3d8bfedaf35d737'
PROTOCOL=ROOT/'research/uncertainty/CRYOLIKE-BASELINE-PROTOCOL.md'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def comparison(cases,seed):
    n=len(cases['neural']['cross_correlation_M'])
    rng=np.random.default_rng(seed);draw=rng.integers(n,size=(20000,n));answer={}
    for metric in ['cross_correlation_M','integrated_log_score']:
        rows={}
        for method in ['gaussian','voxel','reference']:
            d=cases[method][metric]-cases['neural'][metric]
            boot=d[draw].mean(axis=1)
            rows[method+'_minus_neural']=dict(estimate=float(d.mean()),
                percentile_95=np.quantile(boot,[.025,.975]).tolist(),
                percentile_bonferroni_18=np.quantile(boot,[.05/36,1-.05/36]).tolist())
        answer[metric]=rows
    return dict(seed=seed,exposure_resamples=20000,contrasts_per_metric_family=18,
        comparisons=answer,scope='Exploratory paired bootstrap approximation on reused exposures, not density calibration.')


def main():
    torch.set_num_threads(2)
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=AUTHOR,text=True).strip()!=COMMIT:
        raise RuntimeError('Wrong CryoLike author commit')
    if subprocess.check_output(['git','diff','HEAD','--','src'],cwd=AUTHOR):
        raise RuntimeError('Author source modified')
    files=[Path(__file__),PROTOCOL,ROOT/'src/fourier_splats/uq_cryolike_adapter.py',ROOT/'tests/test_uq_cryolike_adapter.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise RuntimeError('Commit implementation and protocol before scores')
    out=BASE/'cryolike-comparison-v1'
    if out.exists():raise RuntimeError('Preserve previous scores and failures')
    out.mkdir();start=time.perf_counter();summaries=[]
    author_hashes={str(p.relative_to(AUTHOR)):sha(p) for p in sorted((AUTHOR/'src').rglob('*.py'))}
    (out/'author-source.json').write_text(json.dumps(dict(commit=COMMIT,url='https://github.com/flatironinstitute/CryoLike',
        source_hashes=author_hashes,pydantic=importlib.metadata.version('pydantic'),torch=torch.__version__),indent=2)+'\n')
    lockpath=ROOT/'research/uncertainty/confirmation/prediction-v1/locked-models.json'
    lock=json.loads(lockpath.read_text())
    for ds in ['10028','10049','10076']:
        data=ROOT/f'data/uncertainty/confirmation/prediction-v1/{ds}'
        selection=ROOT/f'research/uncertainty/confirmation/prediction-v1/{ds}-selection.csv'
        rows=list(csv.DictReader(selection.open()));first={}
        for pos,row in enumerate(rows):first.setdefault(row['source_group'],pos)
        groups=sorted(first)[:128];positions=np.array([first[group] for group in groups])
        assert len(groups)==min(128,len(first)) and len(groups)>1
        metadata=np.load(data/'metadata.npz');indices=np.load(data/'indices.npy')
        np.testing.assert_array_equal(indices,[int(row['metadata_source_index']) for row in rows])
        parameters=metadata['ctf'][positions]
        field=float(parameters[0,0]*parameters[0,1]);pixel=field/64
        np.testing.assert_allclose(parameters[:,0]*parameters[:,1],field)
        images=np.load(data/'images.npy',mmap_mode='r')[positions].copy()
        maps={};map_inputs={}
        for method in ['gaussian','voxel','neural']:
            halves=[];hashes={}
            for half in [0,1]:
                if method=='neural':
                    p=BASE/f'neural-reconstruction/{ds}/half{half}/reconstruct.{lock["epochs"][ds]}.mrc'
                    with mrcfile.open(p) as m:volume=m.data.copy()
                else:
                    p=BASE/f'group-reconstruction/{ds}/{method}-half{half}.npz'
                    volume=volume_from_fourier(np.load(p)['fourier'])
                digest=sha(p);assert digest==lock['models'][ds]['files'][str(p.relative_to(ROOT))]
                hashes[str(p.relative_to(ROOT))]=digest;halves.append(volume)
            maps[method]=np.mean(halves,axis=0);map_inputs[method]=hashes
        emd={'10028':'2660','10049':'6487','10076':'8434'}[ds]
        p=ROOT/f'data/uncertainty/references/emd_{emd}.map'
        reference=VoxelReference.from_mrc(p,box=64)
        np.testing.assert_allclose(reference.field_A,field,rtol=1e-6)
        maps['reference']=reference.volume;map_inputs['reference']={str(p.relative_to(ROOT)):sha(p)}
        common_hashes={str(p.relative_to(ROOT)):sha(p) for p in [lockpath,selection,
            data/'images.npy',data/'metadata.npz',data/'indices.npy',data/'manifest.json']}
        grid=make_grid();im=prepare_images(images,pixel,grid);transfer=prepare_ctf(parameters,field,grid)
        for grid_index,distance in enumerate([.4,.2]):
            views=ViewingAngles.from_viewing_distance(distance);successful={};case_records=[]
            for method in ['gaussian','voxel','neural','reference']:
                case_start=time.perf_counter();name=f'{ds}-grid{grid_index}-{method}';path=out/f'{name}.json'
                record=dict(complete=False,dataset=ds,method=method,viewing_distance=distance,
                    source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]),
                    author_source_sha256=sha(out/'author-source.json'),author_commit=COMMIT,
                    input_hashes={**common_hashes,**map_inputs[method]},groups=groups,
                    positions=positions.tolist(),indices=indices[positions].tolist(),particles=len(groups),
                    field_A=field,pixel_size_A=pixel,frequency_radius=12,
                    templates=views.n_angles,inplanes=128,radial_shells=grid.n_shells,
                    max_shift_pixels_per_axis=2.,shift_points_per_axis=5,progress=[],
                    scope='Original CryoLike normalized structural scores on reused experimental particles; no e-value or density coverage.')
                def save():path.write_text(json.dumps(record,indent=2)+'\n')
                save()
                try:
                    templates=prepare_templates(maps[method],pixel,grid,views)
                    record['template_setup_seconds']=time.perf_counter()-case_start
                    def progress(row):
                        if row['image_end']==len(groups):
                            record['progress'].append(row);save()
                    scored=score(im,templates,transfer,wall_seconds=900.,callback=progress)
                    if not all(np.isfinite(v).all() for v in scored.values()):
                        raise FloatingPointError('Nonfinite original CryoLike output')
                    dest=out/f'{name}-arrays.npz'
                    np.savez_compressed(dest,**scored,positions=positions,indices=indices[positions],
                        azimus=views.azimus.numpy(),polars=views.polars.numpy(),gammas=views.gammas.numpy(),
                        viewing_weights=views.weights_viewing.numpy(),polar_radius=grid.radius_shells,
                        polar_angle=grid.theta_shell,polar_weights=grid.weight_points)
                    record.update(complete=True,scientific_run_complete=True,arrays_file=str(dest.relative_to(ROOT)),
                        arrays_sha256=sha(dest),seconds=time.perf_counter()-case_start,
                        means={k:float(scored[k].mean()) for k in ['cross_correlation_M','integrated_log_score']},
                        peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                    successful[method]=scored
                except Exception as error:
                    record.update(complete=True,scientific_run_complete=False,error=repr(error),seconds=time.perf_counter()-case_start)
                save();print(name,{k:record.get(k) for k in ['scientific_run_complete','means','error','seconds']},flush=True)
                case_records.append(dict(method=method,path=str(path.relative_to(ROOT)),sha256=sha(path),
                    scientific_run_complete=record['scientific_run_complete']))
            summary=dict(dataset=ds,viewing_distance=distance,cases=case_records)
            if len(successful)==4:summary.update(comparison(successful,940001+int(ds)+1000*grid_index))
            else:summary['comparisons_skipped']='At least one declared method failed; no partial winner selected.'
            summaries.append(summary)
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summaries,
        seconds=time.perf_counter()-start,scope='Exploratory external scoring comparison; not a reconstruction or confidence guarantee.'),indent=2)+'\n')


if __name__=='__main__':main()
