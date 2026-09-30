#!/usr/bin/env python3
"""Post-audit completed pose estimators with a residual/pose cross-term bound."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.stats import norm
from fourier_splats.uq_data import particle_geometry
from fourier_splats.uq_pose_operator import PolynomialPoseFieldOperator
from fourier_splats.uq_joint_bias import residual_pose_cross_bound,joint_density_pose_bias,sharp_cube_cubic_coefficients
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_provenance import source_snapshot

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--sharp-cubic',action='store_true');args=p.parse_args()
    source=ROOT/args.fit;prior=json.loads(source.read_text())
    if not prior.get('complete'):raise ValueError('Completed fit required')
    dataset=prior['dataset'];target=prior['target'];fit=prior['fit']
    out=BASE/('joint-bias-sharp-audit' if args.sharp_cubic else 'joint-bias-audit')/source.parent.name;out.mkdir(parents=True,exist_ok=True);path=out/source.name
    if path.exists():raise RuntimeError('Preserve previous audit')
    result={'stage':'Exploratory same-weight continuous cross-term audit; old spectral event reused, nuisance assumptions unchanged',
            'config':vars(args),'dataset':dataset,'target':target,'source_snapshot':source_snapshot(ROOT,Path(__file__)),
            'source_fit_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'complete':False}
    def save():path.write_text(json.dumps(result,indent=2)+'\n')
    save();start=time.perf_counter()
    try:
        g=particle_geometry(ROOT,dataset,'inference_half0',radius=5,count=128,seed=prior['source_geometry_config']['seed'])
        wp=source.with_name(source.stem+'-weights.npz');saved=np.load(wp);w=saved['weights'];noise=float(saved['noise_std'])
        result['source_weights_sha256']=hashlib.sha256(wp.read_bytes()).hexdigest()
        np.testing.assert_array_equal(saved['indices'],g['indices'])
        initial=np.load(BASE/'continuous-quadrature-optimized'/f'{dataset}-{target}-0.07-weights.npz')['weights']
        op=PolynomialPoseFieldOperator(g['k'],g['q'],g['ctf'],initial,noise,np.deg2rad(prior['rotation_radius_degrees']),
               prior['translation_radius_A']/g['field_A'],order=fit.get('pose_quadrature_order',32))
        scaling=op.establish_group_scaling();scale_sum=float(np.where(scaling['group_scales']>0,scaling['group_scales'],1.).sum())
        np.testing.assert_allclose(scale_sum,fit['initial_group_scale_sum'],rtol=1e-8)
        op.set_weights(w);centers=[[0,0,0]] if target=='center' else [[0,0,.08],[0,0,-.08]];signs=[1] if target=='center' else [1,-1]
        cross=residual_pose_cross_bound(op,w,noise,centers,signs,.07);vector=cross.pop('cross_vector')
        bound=joint_density_pose_bias(fit['density_bias']/2,fit['pose_polynomial_bias']/3,cross['norm_upper'],np.sqrt(scale_sum),2.,1.,fit['cubic_bias']/3)
        np.testing.assert_allclose(bound['triangle_bias_upper'],fit['bias'],rtol=1e-8)
        if args.sharp_cubic:
            wr=w.reshape(op.n,2*op.nq);amp=np.hypot(wr[:,:op.nq],wr[:,op.nq:])
            sharp=sharp_cube_cubic_coefficients(g['k'],g['q'],g['ctf']/noise,op.angle,op.shift,3.)
            cubic=float(np.sum(sharp*amp))
            result['cubic_refinement']={'original_bias':fit['cubic_bias'],'sharp_bias':cubic}
            bound=joint_density_pose_bias(fit['density_bias']/2,fit['pose_polynomial_bias']/3,cross['norm_upper'],np.sqrt(scale_sum),2.,1.,cubic/3)
        half=bias_aware_half_width_stable(fit['noise_sd'],bound['bias_upper'],fit['alpha_noise']);no_data=2*fit['target_norm']
        checks=[];diagnostic=BASE/'pose-optimized-diagnostics'/source.parent.name/source.name
        if diagnostic.exists():
            previous=json.loads(diagnostic.read_text());result['source_diagnostic_sha256']=hashlib.sha256(diagnostic.read_bytes()).hexdigest()
            for row in previous['reference_checks']:
                sd=fit['noise_sd'];h=min(half,no_data)
                # All development cases passed here have a non-vacuous original interval.
                if half>=no_data or previous['uses_no_data']:raise ValueError('Reference recentering needed for fallback case')
                coverage=norm.cdf((h-row['actual_bias'])/sd)-norm.cdf((-h-row['actual_bias'])/sd)
                power=norm.cdf((np.sign(row['true_target'])*row['expected_center']-h)/sd)
                checks.append({'scenario':row['scenario'],'analytic_coverage':float(coverage),'correct_sign_probability':float(power)})
            result['existing_feasible_over_new_upper']=max(x['feasible_bias_numerical'] for x in previous['candidates'])/bound['bias_upper']
        result.update(complete=True,seconds=time.perf_counter()-start,cross=cross,bound=bound,
            rotation_radius_degrees=prior['rotation_radius_degrees'],translation_radius_A=prior['translation_radius_A'],
            old_selected_relative_half_width=prior['selected_relative_half_width'],new_half_width=half,
            new_selected_relative_half_width=float(min(1.,half/no_data)),relative_width_change=float(half/fit['half_width']),
            alpha_noise=fit['alpha_noise'],alpha_numerical=fit['alpha_numerical'],reference_checks=checks,
            numerical_failure=bool(result.get('existing_feasible_over_new_upper',0.)>1.000001))
        np.savez(path.with_suffix('.npz'),residual_pose_cross=vector);save()
        print('DONE',source.parent.name,source.stem,result['new_selected_relative_half_width'],result['relative_width_change'],flush=True)
        if result['numerical_failure']:raise AssertionError('Existing feasible bias exceeds new upper; all results retained')
    except Exception as exc:
        result.update(error=repr(exc),seconds=time.perf_counter()-start);save();raise


if __name__=='__main__':main()
