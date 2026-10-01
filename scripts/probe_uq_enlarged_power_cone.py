#!/usr/bin/env python3
"""Larger-view nonnegative power-cone witness diagnostic; no fresh observations."""
import hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from fourier_splats.uq_data import particle_geometry,VoxelReference
from fourier_splats.uq_continuous import cell_forward
from fourier_splats.uq_power_cone_witness import cone_witness,power_witness_log_upper
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'
FILES=['scripts/probe_uq_enlarged_power_cone.py','src/fourier_splats/uq_power_cone_witness.py',
 'src/fourier_splats/uq_power_cone.py','tests/test_power_cone_witness.py',
 'research/uncertainty/paired-power-v1/ENLARGED-CONE-PROTOCOL.md']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 for name in FILES:
  if subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)!=(ROOT/name).read_bytes():raise ValueError('Commit protocol and code first')
 out=BASE/'paired-power-enlarged-cone-v1'
 if out.exists():raise ValueError('Preserve previous probe')
 out.mkdir();start=time.perf_counter()
 meta=dict(complete=False,cases=[],source_sha256={n:sha(ROOT/n) for n in FILES},
  git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
  scope='Numerical nonnegative mixture witnesses and oracle expected-log upper diagnostics, not rejection probabilities.',input_hashes={})
 def save():(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')
 save()
 try:
  for ds,emd in [('10028','2660'),('10049','6487'),('10076','8434')]:
   rp=BASE/f'mixture-validation-preflight-v2/{ds}.json';record=json.loads(rp.read_text());ap=ROOT/record['arrays_file']
   assert sha(ap)==record['arrays_sha256']
   with np.load(ap) as f:old=f['templates'].copy()
   for name,h in record['input_hashes'].items():
    if sha(ROOT/name)!=h:raise ValueError('An original input changed: '+name)
   fp=BASE/f'pilot-selected-fixed-v1/{ds}-pilot_region_1.json';seed=json.loads(fp.read_text())['config']['seed']
   g=particle_geometry(ROOT,ds,'inference_half0',radius=12,count=128,seed=seed)
   np.testing.assert_array_equal(g['indices'],record['original_indices'])
   reference_path=ROOT/f'data/uncertainty/references/emd_{emd}.map'
   rho=VoxelReference.from_mrc(reference_path,box=64).volume.ravel();rho/=np.linalg.norm(rho)
   x=(np.arange(64)+.5)/64-.5;z,y,x=np.meshgrid(x,x,x,indexing='ij');xyz=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=1)
   center=np.asarray(record['target']['center_fraction_field']);bump=np.exp(-.5*np.sum((xyz-center)**2,axis=1)/(20/g['field_A'])**2)
   densities=[rho,rho*(1-.5*bump),rho*(1-bump)]
   rng=np.random.default_rng(261003+int(ds));rot=Rotation.random(4096,random_state=rng).as_matrix()
   plane=np.pad(g['q'][0],((0,0),(0,1)));newk=np.einsum('qi,nij->nqj',plane,rot)
   k=np.concatenate([g['k'][:64],newk]);n,nq=k.shape[:2];fourier=[];discrepancies=[]
   for candidate,density in enumerate(densities):
    values=cell_forward(k,np.ones((n,nq)),density,64,1.).reshape(n,2*nq)
    values=values[:,:nq]+1j*values[:,nq:];fourier.append(values)
    diff=float(np.max(abs(values[:64]-old[candidate])));discrepancies.append(diff)
    if diff>1e-8:raise ArithmeticError(f'Original template mismatch: {ds}/{candidate}/{diff}')
   fourier=np.stack(fourier);powers=abs(fourier)**2;signal=np.mean(abs(old[0])**2,axis=0)
   transfer=g['ctf'][::16]**2/record['noise_std']**2
   arrays=dict(fourier=fourier,rotations=rot,signal=signal,transfer_squared=transfer,original_k=g['k'][:64],q=g['q'][0])
   for candidate,name in enumerate(['true_map','region_half_removed','region_removed']):
    for added in [0,256,1024,4096]:
     if time.perf_counter()-start>1800:raise TimeoutError('Declared 30-minute budget reached')
     tick=time.perf_counter();key=f'{name}_{added}'
     try:
      result=cone_witness(powers[candidate,:64+added],signal)
      coefficients=result.pop('coefficients');approximation=result.pop('approximation')
      arrays[key+'_coefficients']=coefficients;arrays[key+'_approximation']=approximation
      profiles={str(v):[power_witness_log_upper(signal,approximation,tr,v) for tr in transfer] for v in [1.,2.,4.]}
      result.update(expected_log_upper_by_profile=profiles,
         expected_log_128_upper={v:16*sum(values) for v,values in profiles.items()})
      case=dict(dataset=ds,candidate=name,added_views=added,total_views=64+added,seconds=time.perf_counter()-tick,**result)
     except Exception as error:case=dict(dataset=ds,candidate=name,added_views=added,error=repr(error),seconds=time.perf_counter()-tick)
     meta['cases'].append(case);save()
     print(ds,name,added,'scaled residual',case.get('maximum_scaled_residual'),'log128',case.get('expected_log_128_upper'),'error',case.get('error'),flush=True)
   path=out/f'{ds}-arrays.npz';np.savez_compressed(path,**arrays)
   meta.setdefault('arrays',{})[ds]=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),original_template_max_discrepancies=discrepancies)
   for p in [rp,ap,fp,reference_path]:meta['input_hashes'][str(p.relative_to(ROOT))]=sha(p)
   save()
  meta.update(complete=True,seconds=time.perf_counter()-start);save()
 except Exception as exc:
  if 'arrays' in locals():
   path=out/f'{ds}-partial-arrays.npz';np.savez_compressed(path,**arrays);meta['partial_arrays_sha256']=sha(path)
  meta.update(error=repr(exc),seconds=time.perf_counter()-start);save();raise
 print('COMPLETE',len(meta['cases']),meta['seconds'],flush=True)

if __name__=='__main__':main()
