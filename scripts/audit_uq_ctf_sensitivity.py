#!/usr/bin/env python3
"""Retain higher-band CTF/gain sensitivity widths for all prior fine-target fits."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.physics import ctf
from fourier_splats.uq_data import particle_geometry,VoxelReference
from fourier_splats.uq_continuous import cell_forward,cell_target_coefficients
from fourier_splats.uq_ctf_uncertainty import ctf_uncertainty_envelope,continuous_ctf_bias_bound
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot
from audit_uq_grid_refinement import model

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
MAPS={'10028':'2660','10049':'6487','10076':'8434'}
SETTINGS=[('nominal',{}),('defocus100A',{'defocus_radius_A':100.}),
          ('defocus500A',{'defocus_radius_A':500.}),
          ('mixed',{'defocus_radius_A':100.,'astigmatism_angle_radius_degrees':2.,
                    'phase_radius_degrees':1.,'relative_gain_radius':.05,'relative_B_radius_A2':10.})]


def main():
    out=BASE/'higher-band-ctf-sensitivity.json'
    if out.exists():raise RuntimeError('Preserve prior CTF outcomes')
    start=time.perf_counter();snapshot=source_snapshot(ROOT,Path(__file__),['scripts/audit_uq_grid_refinement.py'])
    result={'stage':'Post-review development, fixed-pose weights with bounded CTF sensitivities; known simulation noise',
            'source_snapshot':snapshot,'complete':False,'records':[],
            'assumptions':'Declared nuisance radii, not experimentally estimated confidence bounds; fixed voltage/Cs/amplitude contrast; no orientation error in these six designs.'}
    for dataset in MAPS:
        folder=BASE/('continuous-quadrature-high-band-probe' if dataset=='10028' else 'continuous-quadrature-high-band-additional')
        source=folder/f'{dataset}.json';prior=json.loads(source.read_text());cfg=prior['config']
        g=particle_geometry(ROOT,dataset,'inference_half0',radius=cfg['frequency_radius'],count=cfg['particles'],seed=cfg['seed'])
        parameters=np.load(ROOT/f'data/{dataset}/metadata.npz')['ctf'][g['indices']]
        checkpoint=np.load(BASE/'representation'/dataset/'real_particles-spacing-2.0.npz')
        op,pilot,_,noise=model(g,checkpoint,24);pilot=op.expand(pilot)
        rho=VoxelReference.from_mrc(ROOT/f'data/uncertainty/references/emd_{MAPS[dataset]}.map',box=64).volume.ravel()
        rho/=np.linalg.norm(rho);nominal=cell_forward(g['k'],g['ctf'],pilot,24,noise)
        for row in prior['targets']:
            name=row['target'];width=row['width_fraction_field'];fit=row['fit']
            centers=[[0,0,0]] if name=='center' else [[0,0,.08],[0,0,-.08]];signs=[1] if name=='center' else [1,-1]
            saved=np.load(folder/f'{dataset}-{name}-{width}-weights.npz');w=saved['weights']
            np.testing.assert_array_equal(saved['indices'],g['indices']);np.testing.assert_allclose(float(saved['noise_std']),noise)
            pilot_target=float(cell_target_coefficients(24,centers,signs,width)@pilot)
            truth=float(cell_target_coefficients(64,centers,signs,width)@rho)
            for label,setting in SETTINGS:
                envelope=ctf_uncertainty_envelope(g['q'][0]/g['field_A'],parameters,**setting)
                extra=continuous_ctf_bias_bound(w,noise,envelope['absolute_transfer_error'],3.)
                half=bias_aware_half_width_stable(fit['noise_sd'],fit['bias']+extra)
                no_data=2*fit['target_norm'];fallback=bool(half>=no_data);selected=min(half,no_data)
                checks=[]
                for sign in [-1.,1.]:
                    perturbed=parameters.copy();perturbed[:,2:4]+=sign*setting.get('defocus_radius_A',0.)
                    perturbed[:,4]+=sign*setting.get('astigmatism_angle_radius_degrees',0.)
                    perturbed[:,8]+=sign*setting.get('phase_radius_degrees',0.)
                    q2=np.sum((g['q'][0]/g['field_A'])**2,axis=1)
                    transfer=ctf(g['q'][0]/g['field_A'],perturbed)*(1+sign*setting.get('relative_gain_radius',0.))*np.exp(sign*setting.get('relative_B_radius_A2',0.)*q2[None]/4)
                    signal=cell_forward(g['k'],transfer,rho,64,noise)
                    mean=pilot_target if fallback else float(pilot_target+w@(signal-nominal))
                    sd=0. if fallback else fit['noise_sd'];bias=mean-truth
                    coverage=float(norm.cdf((selected-bias)/sd)-norm.cdf((-selected-bias)/sd)) if sd else float(abs(bias)<=selected)
                    power=float(norm.cdf((np.sign(truth)*mean-selected)/sd)) if sd else float(np.sign(truth)*mean>selected)
                    checks.append({'coherent_sign':sign,'actual_bias':bias,'analytic_coverage':coverage,'correct_sign_probability':power})
                result['records'].append({'dataset':dataset,'target':name,'sigma_A':width*g['field_A'],
                    'frequency_radius':cfg['frequency_radius'],'setting':label,'radii':setting,'ctf_bias_addition':extra,
                    'fixed_pose_density_bias':fit['bias'],'uncapped_half_width':half,'uses_no_data':fallback,
                    'selected_relative_half_width':selected/no_data,'reference_checks':checks,
                    'source_result_sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
                out.write_text(json.dumps(result,indent=2)+'\n')
                print(dataset,name,label,selected/no_data,flush=True)
    result.update(complete=True,seconds=time.perf_counter()-start);out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
