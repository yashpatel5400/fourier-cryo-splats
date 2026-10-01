#!/usr/bin/env python3
"""Declared post hoc replay of all frozen alignment trials, including true poses."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_continuous import cell_forward, cell_target_coefficients, continuous_certificate
from fourier_splats.uq_continuous_gaussian import continuous_gaussian
from fourier_splats.uq_continuous_quadrature import QuadratureObservationGram
from fourier_splats.uq_end_to_end import dephase_observations
from fourier_splats.uq_group_noise import cell_density_distance
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_refitting_diagnostics import dephase_adjoint_weights, realized_class_envelopes

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
TARGETS = [('center', [[0., 0., 0.]], [1.]), ('contrast', [[0., 0., .08], [0., 0., -.08]], [1., -1.])]
RADII = [0., .5, 1., 2.]
METHODS = ['fixed_folded', 'gaussian_tau2_fixed', 'gaussian_tau2_pose', 'gaussian_tau1_fixed', 'gaussian_tau1_pose']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(ds, args):
    start = time.perf_counter(); out = BASE/args.output/ds
    if out.exists():
        raise ValueError(f'Preserve prior reanalysis: {out}')
    out.mkdir(parents=True)
    original = BASE/'end-to-end-local-pose-v1'/ds
    cal = BASE/'local-alignment-calibration-v1'/ds
    summary = json.loads((original/'summary.json').read_text())
    gpath = cal/'generators.npz'
    if sha(gpath) != summary['calibration_hashes']['generators.npz']:
        raise ValueError('Generator hash changed')
    with np.load(gpath) as data:
        g = {k: data[k].copy() for k in data.files}
    noise, field = float(g['noise_std']), float(g['field_A'])
    rho, pilot, signal = g['truth'], g['pilot'], g['noiseless_signal']
    density_distance, distance_pad = cell_density_distance(pilot, 24, rho, 64)
    pilot_signal = cell_forward(g['k'], g['ctf'], pilot, 24, noise)
    gt = QuadratureObservationGram(g['k'], g['ctf'], noise, order=40, preconditioner_rank=1024)
    gt.nthreads = 1
    truths, pilots, true_fits, weights_archive = {}, {}, {}, {}
    max_center_error = 0.; rows = []; diagnoses = []; verified = []
    source_paths = [Path(__file__), ROOT/'src/fourier_splats/uq_refitting_diagnostics.py',
        ROOT/'research/uncertainty/REFITTING-REANALYSIS-PROTOCOL.md']
    metadata = dict(complete=False, dataset=ds, git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_paths},
        original_summary_sha256=sha(original/'summary.json'), generator_sha256=sha(gpath),
        scope='Post hoc saved-data diagnostic; realized envelopes use unknown true poses and are not deployable pose-set intervals.',
        density_distance=density_distance, density_distance_squared_roundoff_diagnostic_pad=distance_pad,
        radius_grid=RADII, errors=[])
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for target, centers, signs in TARGETS:
        truths[target] = float(cell_target_coefficients(64,centers,signs,.07)@rho)
        pilots[target] = float(cell_target_coefficients(24,centers,signs,.07)@pilot)
        audit = continuous_certificate(gt,centers,signs,.07,2.,alpha=.05,maxiter=100,rtol=.005)
        true_fits[target] = {'fixed_folded':audit}
        for tau in [2.,1.]:
            true_fits[target][f'gaussian_tau{tau:g}_fixed'] = continuous_gaussian(gt,centers,signs,.07,tau,alpha=.05,density_radius=2.)
        for name, fit in true_fits[target].items():
            weights_archive[target+'_'+name] = fit['weights']
            if not fit['converged']:
                metadata['errors'].append(f'True-pose fit nonconverged: {target}/{name}')
    np.savez_compressed(out/'true-pose-weights.npz', **weights_archive)
    (out/'true-pose-fits.json').write_text(json.dumps({t:{m:{k:v for k,v in f.items() if k!='weights'} for m,f in fits.items()} for t,fits in true_fits.items()},indent=2)+'\n')

    def add_row(rep, template, target, method, mode, w, y, expected, prediction, original_half, residual, old_center=None):
        nonlocal max_center_error
        sd = float(np.linalg.norm(w)); pilot_target = pilots[target]; truth = truths[target]
        center = float(pilot_target+w@(y-prediction))
        bias = float(pilot_target+w@(expected-prediction)-truth)
        error = center-truth
        if old_center is not None:
            discrepancy = abs(center-old_center); max_center_error=max(max_center_error,discrepancy)
            if discrepancy > 1e-7:
                raise ArithmeticError(f'Saved centre does not replay: {discrepancy}')
        row = dict(dataset=ds,replicate=rep,template=template,target=target,method=method,image_mode=mode,
            truth=truth,pilot=pilot_target,pilot_error=pilot_target-truth,center=center,error=error,
            noise_sd=sd,standardized_error=error/sd,realized_signal_bias=bias,
            standardized_signal_bias=bias/sd,projected_noise=center-truth-bias,
            standardized_projected_noise=(center-truth-bias)/sd,residual_norm=residual,
            original_half_width=original_half,original_covered=abs(error)<=original_half)
        for B in RADII:
            tag=f'B{B:g}'; half=float(bias_aware_half_width_stable(sd,B*residual,.05))
            row[tag+'_half_width']=half;row[tag+'_covered']=abs(error)<=half
            row[tag+'_contains_generator']=density_distance<=B
        rows.append(row)

    for entry in summary['records'][:args.limit]:
        rep=entry['replicate']; rp=original/f'replicate-{rep:03d}.json'; ap=rp.with_suffix('.npz')
        if sha(rp)!=entry['record_sha256']:
            raise ValueError('Trial JSON changed')
        old=json.loads(rp.read_text())
        if not old['complete'] or sha(ap)!=old['arrays_sha256']:
            raise ValueError('Incomplete or changed trial arrays')
        verified.append(dict(replicate=rep,json_sha256=sha(rp),arrays_sha256=sha(ap)))
        with np.load(ap) as data:
            arrays={k:data[k].copy() for k in data.files}
        for target,centers,signs in TARGETS:
            for method,fit in true_fits[target].items():
                residual=fit['bias']/2 if method=='fixed_folded' else fit['residual_norm']
                for mode,ykey in [('same_image','y_a'),('independent_image','y_b')]:
                    add_row(rep,'true_pose',target,method,mode,fit['weights'],arrays[ykey].ravel(),signal.ravel(),pilot_signal,fit['half_width'],residual)
        for template_record in old['templates']:
            template=template_record['template']; rotations=arrays[template+'_rotations']; shifts=arrays[template+'_shifts_A']
            khat=g['k']@rotations
            gf=QuadratureObservationGram(khat,g['ctf'],noise,order=40,preconditioner_rank=0);gf.nthreads=1
            prediction=cell_forward(khat,g['ctf'],pilot,24,noise)
            expected=dephase_observations(signal,g['q'],shifts,field)
            true_pilot=dephase_observations(pilot_signal.reshape(signal.shape),g['q'],shifts,field)
            observed={mode:dephase_observations(arrays[ykey],g['q'],shifts,field) for mode,ykey in [('same_image','y_a'),('independent_image','y_b')]}
            for fit_record in template_record['fits']:
                target=fit_record['target']; _,centers,signs=next(t for t in TARGETS if t[0]==target)
                prefix=template+'_'+target
                lookup={(r['method'],r['image_mode']):r for r in old['intervals'] if r['template']==template and r['target']==target}
                for method in METHODS:
                    key=prefix+('_audit_weights' if method=='fixed_folded' else '_'+method+'_weights')
                    w=arrays[key]
                    if method=='fixed_folded':
                        residual=fit_record['fixed_audit']['bias']/2
                        v=dephase_adjoint_weights(w,g['q'],shifts,field)
                        diag=realized_class_envelopes(gt,gf,w,v,centers,signs,.07,2.,float(w@(true_pilot-prediction)))
                        no_data=2*fit_record['fixed_audit']['target_norm']; mixed=fit_record['mixed_audit']
                        diag.update(dataset=ds,replicate=rep,template=template,target=target,
                            pose_inside_calibrated_ball=template_record['pose_inside_calibrated_ball'],
                            no_data_half_width=no_data,fixed_bias=fit_record['fixed_audit']['bias'],
                            original_deterministic_bias=fit_record['deterministic_pose_bias'],
                            original_mixed_half_width=mixed['half_width'],
                            first_order_only_half_width=mixed['density_bias']+mixed['subgaussian_tail'],
                            first_order_tail=mixed['subgaussian_tail'],
                            noise_sd=float(np.linalg.norm(w)))
                        diagnoses.append(diag)
                    else:
                        residual=next(f['residual_norm'] for f in fit_record['gaussian_fits'] if f['method']==method)
                    for mode,y in observed.items():
                        interval=lookup[method,mode]
                        add_row(rep,template,target,method,mode,w,y,expected,prediction,interval['raw_half_width'],residual,interval['raw_center'])
        if rep%10==0:
            print(ds,rep,'seconds',round(time.perf_counter()-start,2),'max_center_error',max_center_error,flush=True)
    for name,values in [('rows',rows),('realized-envelopes',diagnoses)]:
        with (out/(name+'.csv')).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(values[0]));writer.writeheader();writer.writerows(values)
    metadata.update(complete=True,replicates=len(verified),rows=len(rows),realized_envelopes=len(diagnoses),
        max_saved_center_discrepancy=max_center_error,verified_inputs=verified,seconds=time.perf_counter()-start,
        true_pose_fit_count=6,true_pose_fits_converged=not metadata['errors'],
        truths=truths,pilot_targets=pilots,outputs={p.name:sha(p) for p in out.iterdir() if p.name!='metadata.json'})
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(ds,'COMPLETE',len(rows),'rows',metadata['seconds'],'seconds',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--output',default='refitting-bias-reanalysis-v1')
    p.add_argument('--limit',type=int,default=200);a=p.parse_args()
    if a.limit!=200 and a.output=='refitting-bias-reanalysis-v1':
        raise ValueError('Development probes require a different output directory')
    for ds in a.datasets.split(','):run(ds,a)

if __name__=='__main__':main()
