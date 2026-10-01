#!/usr/bin/env python3
"""Retrospective saved-draw diagnostics; changes no gate or fitted proposal."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import logsumexp
from fourier_splats.uq_pose_importance import PoseProposal
from fourier_splats.uq_pose_catalog import CatalogPoseProposal

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--attempt',choices=['adaptive','catalog'],required=True)
    args=p.parse_args();attempt=args.attempt
    out=BASE/f'{attempt}-pose-autopsy-v1'
    if out.exists():raise ValueError('Preserve every autopsy attempt')
    for ds in ['10028','10049','10076']:
        assert json.loads((BASE/f'{attempt}-pose-integration-v1'/ds/'summary.json').read_text())['complete']
    out.mkdir();result=dict(complete=False,attempt=attempt,source_sha256=sha(Path(__file__)),
        protocol_sha256=sha(ROOT/'research/uncertainty/PLANTED-TRUTH-INTEGRATION-AUDIT.md'),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),datasets=[])
    def save():(out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    save()
    for ds in ['10028','10049','10076']:
        directory=BASE/f'{attempt}-pose-integration-v1'/ds
        j=json.loads((directory/'summary.json').read_text())
        old=json.loads((BASE/'adaptive-pose-integration-v1'/ds/'summary.json').read_text())
        image_path=BASE/'matched-information-ledger-v1'/ds/'replayed-images.npz'
        with np.load(image_path) as f:
            true_quat=Rotation.from_matrix(f['rotations']).as_quat()
            noise=f['noise']; means=f['means']
        catalog=np.load(directory/'catalogue-quaternions.npy') if attempt=='catalog' else None
        row=dict(dataset=ds,input_hashes={str(x.relative_to(ROOT)):sha(x) for x in [directory/'summary.json',image_path]},
                 images=[],stack_brackets=[],failure_strata=[])
        result['datasets'].append(row)
        for c,o in zip(j['cases'],old['cases']):
            assert (c['position'],c['state'],c['index'])==(o['position'],o['state'],o['index'])
            fit=o['proposal'];local=PoseProposal(fit['centers'],fit['covariances'],fit['weights'])
            index,state=c['index'],c['state']
            path=ROOT/c['arrays']['path'];assert sha(path)==c['arrays']['sha256']
            with np.load(path) as f:
                logq=f['log_proposal'];kernels=f['log_kernels'];y=f['observation']
                family=f['family'] if attempt=='catalog' else f['component']
                proposal=CatalogPoseProposal(local,catalog,f['catalogue_weights'],np.deg2rad(8)) if attempt=='catalog' else local
            np.testing.assert_allclose(y,means[0,index]-.25*state*means[1,index]+noise[index],atol=0,rtol=0)
            truth_logq=float(proposal.log_density(true_quat[index:index+1])[0])
            truth_kernel=-.5*float(np.sum(abs(noise[index])**2))
            truth_logweight=truth_kernel-truth_logq
            eigenvalues=np.array([m['hessian_eigenvalues'] for m in fit['modes']])
            clipped=(eigenvalues<1/np.deg2rad(30)**2)|(eigenvalues>1/np.deg2rad(.25)**2)
            entry=dict(position=c['position'],index=index,state=state,
                original_ratio_difference=c['between_bank_logratio_difference'],
                ratio_within_tolerance=bool(c['between_bank_logratio_difference']<=.01),
                clipped_eigenvalues=int(clipped.sum()),clipped_modes=int(np.any(clipped,axis=1).sum()),
                clipped_mode_weight=float(np.asarray(fit['weights'])@np.any(clipped,axis=1)),
                truth_log_proposal=truth_logq,truth_log_kernel=truth_kernel,banks=[])
            for bank in [0,1]:
                x=kernels[bank]-logq[bank,:,None];n=len(x)
                total=logsumexp(x,axis=0);weights=np.exp(x-total)
                integral_se=np.sqrt(np.maximum((n*np.sum(weights**2,axis=0)-1)/(n-1),0))
                ratio_se=float(np.sqrt(n/(n-1)*np.sum((weights[:,1]-weights[:,0])**2)))
                original=float(total[state]-np.log(n))
                planted=float(np.logaddexp(logsumexp(x[1:,state]),truth_logweight)-np.log(n))
                matching=weights[:,state];haar=family[bank]==0
                entry['banks'].append(dict(bank=bank,matching_model=state,
                    haar_weight_mass=float(matching[haar].sum()),
                    maximum_weight_from_haar=bool(haar[np.argmax(matching)]),
                    log_integral_delta_standard_errors=integral_se.tolist(),
                    logratio_delta_standard_error=ratio_se,
                    original_matching_log_integral=original,planted_matching_log_integral=planted,
                    planted_jump=planted-original,planted_share=float(np.exp(truth_logweight-planted)/n)))
            se=np.hypot(*[b['logratio_delta_standard_error'] for b in entry['banks']])
            entry['between_bank_ratio_delta_standard_error']=float(se)
            entry['standardized_bank_discrepancy']=entry['original_ratio_difference']/se if se>0 else None
            row['images'].append(entry)
        for state in [0,1]:
            for bank in [0,1]:
                subset=[x['banks'][bank] for x in row['images'] if x['state']==state]
                assert len(subset)==64
                lower=sum(x['original_matching_log_integral'] for x in subset)-np.log(40)
                upper=sum(x['planted_matching_log_integral'] for x in subset)+np.log(40)
                row['stack_brackets'].append(dict(state=state,bank=bank,particles=64,delta_lower=.025,delta_upper=.025,
                    lower=float(lower),upper=float(upper),width=float(upper-lower),
                    interpretation='Pointwise stochastic bracket under joint matched simulation; not simultaneous or deterministic'))
        for flag in [True,False]:
            subset=[x for x in row['images'] if x['ratio_within_tolerance']==flag]
            if not subset:continue
            banks=[b for x in subset for b in x['banks']]
            stats={}
            for name in ['haar_weight_mass','planted_jump','logratio_delta_standard_error','planted_share']:
                v=np.array([b[name] for b in banks]);stats[name]=dict(mean=float(v.mean()),
                    quantiles={str(q):float(np.quantile(v,q)) for q in [0,.1,.5,.9,.95,1]})
            row['failure_strata'].append(dict(ratio_within_tolerance=flag,images=len(subset),banks=len(banks),
                bank_statistics=stats,haar_maximum_fraction=float(np.mean([b['maximum_weight_from_haar'] for b in banks])),
                images_with_clipping=sum(x['clipped_modes']>0 for x in subset),
                mean_clipped_mode_weight=float(np.mean([x['clipped_mode_weight'] for x in subset])),
                median_standardized_discrepancy=float(np.median([x['standardized_bank_discrepancy'] for x in subset]))))
        save();print('AUTOPSY',attempt,ds,'complete',flush=True)
    result['complete']=True;save()


if __name__=='__main__':main()
