#!/usr/bin/env python3
"""Declared optimistic finite-view gate for paired-exposure map compatibility."""
import hashlib,json,time,subprocess
from pathlib import Path
import numpy as np
from fourier_splats.uq_power_cone import FinitePowerCone
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
FILES=['scripts/probe_uq_paired_power.py','src/fourier_splats/uq_power_cone.py','tests/test_power_cone.py',
    'research/uncertainty/paired-power-v1/FINITE-VIEW-PROBE.md','src/fourier_splats/uq_paired_power.py']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 for name in FILES:
  if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():
   raise ValueError('Commit probe protocol and implementation before outcomes')
 out=BASE/'paired-power-finite-view-v4'
 if out.exists():raise ValueError('Preserve earlier probe')
 out.mkdir();start=time.perf_counter();all_rows=[]
 metadata=dict(complete=False,git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
  source_sha256={name:sha(ROOT/name) for name in FILES},cases=[],input_hashes={},
  scope='Oracle expected-log growth on 64 sampled views; no continuous-pose certificate, no simulated coverage or experimental test.')
 def save(): (out/'summary.json').write_text(json.dumps(metadata,indent=2)+'\n')
 save()
 try:
  for ds in ['10028','10049','10076']:
   record_path=BASE/f'mixture-validation-preflight-v2/{ds}.json';record=json.loads(record_path.read_text())
   ap=ROOT/record['arrays_file']
   if sha(ap)!=record['arrays_sha256']:raise ValueError('Source Fourier means changed')
   metadata['input_hashes'][str(record_path.relative_to(ROOT))]=sha(record_path)
   metadata['input_hashes'][str(ap.relative_to(ROOT))]=sha(ap)
   with np.load(ap) as f:means=f['means'].copy()
   q=means.shape[-1]//2;powers=means[...,:q]**2+means[...,q:]**2
   solver=FinitePowerCone(64,q);arrays={}
   for profile in range(0,128,16):
    signal=powers[0,profile].mean(axis=0)
    for candidate,name in enumerate(['true_map','region_half_removed','region_removed','zero_signal']):
     for variance in [1.,2.,4.]:
      if time.perf_counter()-start>1800:raise TimeoutError('Declared 30-minute limit reached')
      tick=time.perf_counter()
      try:
       fitted=solver.solve(powers[candidate,profile],signal,variance)
      except Exception as error:
       metadata['cases'].append(dict(dataset=ds,profile=profile,candidate=name,variance_upper=variance,
           seconds=time.perf_counter()-tick,error=repr(error),status='failed'))
       save();print(ds,profile,name,variance,'FAILED',repr(error),flush=True);continue
      key=f'{profile}_{name}_v{variance:g}';arrays[key+'_weights']=fitted.pop('weights')
      arrays[key+'_raw_weights']=fitted.pop('raw_dimensionless_weights');arrays[key+'_multipliers']=fitted.pop('multipliers')
      case=dict(dataset=ds,profile=profile,original_index=int(record['original_indices'][profile]),candidate=name,
         seconds=time.perf_counter()-tick,**fitted)
      if name=='true_map' and case['expected_log_lower']>1e-6:
       raise ArithmeticError('Matched-map expected-log control failed')
      metadata['cases'].append(case);save()
      print(ds,profile,name,variance,case['status'],round(case['expected_log_lower'],7),round(case['expected_log_dual_upper'],7),flush=True)
   np.savez_compressed(out/f'{ds}-weights.npz',**arrays)
   metadata.setdefault('weights_sha256',{})[ds]=sha(out/f'{ds}-weights.npz');save()
  aggregates=[]
  for ds in ['10028','10049','10076']:
   for name in ['true_map','region_half_removed','region_removed','zero_signal']:
    for variance in [1.,2.,4.]:
     rows=[r for r in metadata['cases'] if r['dataset']==ds and r['candidate']==name and r['variance_upper']==variance]
     assert len(rows)==8
     aggregates.append(dict(dataset=ds,candidate=name,variance_upper=variance,profiles=8,
       expected_log_128_lower=16*sum(r.get('expected_log_lower',0.) for r in rows),
       expected_log_128_dual_upper=(16*sum(r['expected_log_dual_upper'] for r in rows) if all('error' not in r for r in rows) else None),
       failed_profiles=sum('error' in r for r in rows),
       maximum_profile_gap=max((r['gap'] for r in rows if 'gap' in r),default=None),
       nonoptimal_statuses=sum(r['status']!='optimal' for r in rows)))
  metadata.update(complete=True,aggregates=aggregates,seconds=time.perf_counter()-start);save()
 except Exception as exc:
  
  if 'arrays' in locals():
   partial=out/f'{ds}-partial-weights.npz';np.savez_compressed(partial,**arrays)
   metadata['partial_arrays_sha256']=sha(partial)
  metadata.update(error=repr(exc),seconds=time.perf_counter()-start);save();raise
 print('COMPLETE',len(metadata['cases']),metadata['seconds'],flush=True)

if __name__=='__main__':main()
