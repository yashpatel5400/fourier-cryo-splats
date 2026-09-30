#!/usr/bin/env python3
"""Declared local/full-cover disk diagnostic, not new uncertainty calibration."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.special import logsumexp
from fourier_splats.uq_mixture_disks import DiskGaussianOrbit,validate_independent_plane
from fourier_splats.uq_mixture_anchor import scaled_anchor_upper
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
PROTOCOL=ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-DISK-PROTOCOL.md'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    files=[Path(__file__),PROTOCOL,
        ROOT/'research/uncertainty/CONTINUOUS-MIXTURE-DISK-THEORY.md',
        ROOT/'src/fourier_splats/uq_mixture_disks.py',ROOT/'tests/test_uq_mixture_disks.py']
    for p in files:
        if subprocess.check_output(['git','show',f'HEAD:{p.relative_to(ROOT)}'],cwd=ROOT)!=p.read_bytes():
            raise ValueError('Publish source and protocol before calculations')
    out=BASE/'continuous-mixture-disks-v1'
    if out.exists():raise RuntimeError('Preserve previous outcomes')
    out.mkdir();start=time.perf_counter();summary=[]
    def budget():
        if time.perf_counter()-start>1800:raise TimeoutError('Declared total wall budget')
    for ds in ['10028','10049','10076']:
        path=out/f'{ds}.json'
        record=dict(complete=False,dataset=ds,local_cases=[],
            source_snapshot=source_snapshot(ROOT,Path(__file__),[str(p.relative_to(ROOT)) for p in files[1:]]))
        def save():path.write_text(json.dumps(record,indent=2)+'\n')
        save()
        try:
            source=BASE/f'continuous-mixture-v1/{ds}.json'
            local_source=BASE/f'continuous-mixture-curvature-v1/{ds}.json'
            parent=json.loads(source.read_text());local=json.loads(local_source.read_text())
            if not parent['scientific_run_complete'] or not local['scientific_run_complete']:
                raise ValueError('Completed parent diagnostics required')
            array_path=ROOT/parent['arrays_file']
            if sha(array_path)!=parent['arrays_sha256']:raise ValueError('Changed parent arrays')
            arrays=np.load(array_path)
            record['input_hashes']={str(p.relative_to(ROOT)):sha(p) for p in [source,local_source,array_path]}
            validate_independent_plane(arrays['plane'])
            record['frequency_guard_passed']=True
            def orbit(positions):
                return DiskGaussianOrbit(arrays['plane'],arrays['centers'],arrays['sigma'],
                    arrays['coefficients'],arrays['transfer'][positions],arrays['observed'][positions])
            for case in local['cases']:
                budget();begin=time.perf_counter();position=case['position']
                center=np.array(case['center']);half=np.array(case['half_widths_radians'])
                bound=orbit(slice(position,position+1)).cell(center-half,center+half)
                upper=float(bound['envelope'][0]);feasible=case['feasible_log_kernel']
                if upper<feasible-1e-7:raise ArithmeticError('Disk bound excludes existing feasible value')
                record['local_cases'].append(dict(position=position,budget_degrees=case['budget_degrees'],
                    curvature_upper=case['combined_upper'],disk_upper=float(bound['disk_envelope'][0]),
                    retained_upper=min(case['combined_upper'],upper),feasible_log_kernel=feasible,
                    original_gap=case['combined_gap'],new_gap=min(case['combined_upper'],upper)-feasible,
                    disk_primal_gap=float(bound['disk_primal_gap'][0]),seconds=time.perf_counter()-begin))
            save();full=orbit(slice(None))
            envelopes=arrays['leaf_log_envelopes'];anchor=arrays['best_anchor_log_z']
            scores=logsumexp(envelopes-anchor[None,:],axis=1)
            ordered=np.argsort(scores,kind='stable')
            selected=np.unique(np.r_[ordered[np.rint(np.linspace(0,len(ordered)-1,80)).astype(int)],ordered[-16:]])
            new=[];disk=[];gaps=[]
            for index in selected:
                budget();cell=full.cell(arrays['leaf_lower'][index],arrays['leaf_upper'][index])
                new.append(cell['envelope']);disk.append(cell['disk_envelope']);gaps.append(cell['disk_primal_gap'])
            new=np.array(new);updated=envelopes.copy();updated[selected]=np.minimum(updated[selected],new)
            original_upper=scaled_anchor_upper(envelopes.T,anchor)+full.normalization
            updated_upper=scaled_anchor_upper(updated.T,anchor)+full.normalization
            rng=np.random.default_rng(950001+int(ds));angles=rng.uniform(size=(256,3))
            angles[:,0]*=2*np.pi;angles[:,2]*=2*np.pi;angles[:,1]=np.arccos(1-2*angles[:,1])
            random=[]
            for angle in angles:
                budget();random.append(full.log_kernel(angle))
            random=np.array(random);quantiles=np.quantile(random,[.5,.9,.99],axis=0)
            best_center=arrays['support_log_kernels'].max(axis=0)
            local_summary=[]
            for b in [.1,.5,1.,2.,5.,10.]:
                cases=[c for c in record['local_cases'] if c['budget_degrees']==b]
                local_summary.append(dict(budget_degrees=b,
                    old_gap_median=float(np.median([c['original_gap'] for c in cases])),
                    new_gap_median=float(np.median([c['new_gap'] for c in cases]))))
            record.update(selected_cells=len(selected),local_summary=local_summary,
                original_cover_upper=original_upper,updated_cover_upper=updated_upper,
                max_per_image_upper_improvement=float(np.max(envelopes[selected]-updated[selected])),
                positive_improvement_coordinates=int(np.sum(updated[selected]<envelopes[selected]-1e-9)),
                tested_coordinates=int(new.size),
                center_minus_haar_quantile_medians=np.median(best_center[None,:]-quantiles,axis=1).tolist(),
                envelope_minus_best_center_median=float(np.median(envelopes.max(axis=0)-best_center)))
            output_arrays=out/f'{ds}-arrays.npz'
            np.savez_compressed(output_arrays,selected_cells=selected,old_envelopes=envelopes[selected],
                new_envelopes=new,disk_envelopes=np.array(disk),disk_primal_gaps=np.array(gaps),
                haar_angles=angles,haar_log_kernels=random,best_center_log_kernels=best_center,
                original_cell_log_scores=scores)
            record.update(complete=True,scientific_run_complete=True,arrays_file=str(output_arrays.relative_to(ROOT)),
                arrays_sha256=sha(output_arrays),scope='Reused simulations, numerical enclosure diagnostics only.')
            save();summary.append({k:v for k,v in record.items() if k not in ['source_snapshot','local_cases','input_hashes']})
            print(ds,record['local_summary'],original_upper,updated_upper,flush=True)
        except Exception as error:
            record.update(complete=True,scientific_run_complete=False,error=repr(error));save();raise
    (out/'summary.json').write_text(json.dumps(dict(complete=True,records=summary,seconds=time.perf_counter()-start),indent=2)+'\n')


if __name__=='__main__':main()
