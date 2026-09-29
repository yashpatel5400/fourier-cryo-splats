#!/usr/bin/env python3
"""Development-only audit against independent voxel/NUFFT signal generators.

Fits both noise-free semisynthetic projections and real pilot particles, with
tuning separated by source groups. This is a dictionary/compute feasibility
audit; it makes no true-density confidence claim.
"""
import argparse
import json
import time
from pathlib import Path
import mrcfile
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.fft import fftshift,fftn,ifftn,ifftshift
from fourier_splats.uq_data import particle_geometry,particle_observations,VoxelReference
from fourier_splats.uq_physics import image_design,gaussian_pair_features,realify
from fourier_splats.fsc import fsc

ROOT=Path(__file__).resolve().parents[1]
MAPS={'10028':'2660','10049':'6487','10076':'8434'}


def dictionary(radius,spacing):
    g=np.arange(-np.ceil(radius/spacing),np.ceil(radius/spacing)+1)*spacing
    z,y,x=np.meshgrid(g,g,g,indexing='ij')
    k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    half=(k[:,2]>0)|((k[:,2]==0)&(k[:,1]>0))|((k[:,2]==0)&(k[:,1]==0)&(k[:,0]>=0))
    return k[half&(np.sum(k*k,axis=1)<=radius**2)]


def render(c,centers,sigma,box=64):
    g=np.arange(-box//2,box//2);z,y,x=np.meshgrid(g,g,g,indexing='ij')
    k=np.stack([x.ravel(),y.ravel(),z.ravel()],axis=-1)
    fourier=np.empty(len(k),complex)
    for start in range(0,len(k),1024):
        fourier[start:start+1024]=gaussian_pair_features(k[start:start+1024],centers,sigma)@c
    return fftshift(ifftn(ifftshift(fourier.reshape((box,)*3)))).real.astype(np.float32)


def run(dataset,count):
    out=ROOT/'results/uncertainty/development/representation'/dataset
    out.mkdir(parents=True,exist_ok=True)
    geometries={s:particle_geometry(ROOT,dataset,s,radius=6,count=count,seed=39) for s in ['pilot','tune','test']}
    ref=VoxelReference.from_mrc(ROOT/'data/uncertainty/references'/f'emd_{MAPS[dataset]}.map')
    syn={s:ref.fourier(g['k'],g['field_A'])*g['ctf'] for s,g in geometries.items()}
    real={s:particle_observations(ROOT,dataset,g) for s,g in geometries.items()}
    output={'stage':'development dictionary audit, no uncertainty result',
            'dataset':dataset,'particles_per_split':count,'frequencies':len(geometries['pilot']['q'][0]),
            'reference_EMDB':MAPS[dataset],'frequency_radius_bins':6,'variants':[]}
    for spacing,sigma in [(2.,1.),(1.5,.8),(1.,.65)]:
        started=time.perf_counter();centers=dictionary(8,spacing);p=2*len(centers)
        a={s:image_design(g['k'],g['ctf'],centers,sigma).reshape(-1,p) for s,g in geometries.items()}
        gram=a['pilot'].T@a['pilot'];diag=np.diag(gram);median=np.median(diag[diag>1e-9])
        for kind,observed in [('semisynthetic_noiseless',syn),('real_particles',real)]:
            scale=np.sqrt(np.mean(np.abs(observed['pilot'])**2))
            y={s:realify(v/scale).ravel() for s,v in observed.items()}
            rhs=a['pilot'].T@y['pilot'];tuned=[]
            for relative in [1e-6,1e-4,.01,1.]:
                c=cho_solve(cho_factor(gram+relative*median*np.eye(p)),rhs)
                tuning=float(np.mean((a['tune']@c-y['tune'])**2))
                tuned.append((tuning,relative,c))
            tuning,relative,c=min(tuned,key=lambda x:x[0])
            test_mse=float(np.mean((a['test']@c-y['test'])**2))
            volume=render(c,centers,sigma)
            with mrcfile.new(out/f'{kind}-spacing-{spacing}.mrc',overwrite=True) as m:
                m.set_data(volume);m.voxel_size=geometries['pilot']['field_A']/64
            np.savez(out/f'{kind}-spacing-{spacing}.npz',coefficients=c,centers=centers,sigma=sigma,scale=scale)
            record={'data_kind':kind,'spacing':spacing,'sigma':sigma,'coefficients':p,
                    'relative_ridge':relative,'pilot_scale':float(scale),'tune_mse':tuning,
                    'test_mse':test_mse,'test_relative_rmse':float(np.sqrt(test_mse/np.mean(y['test']**2))),
                    'cumulative_seconds':time.perf_counter()-started}
            if kind=='semisynthetic_noiseless':
                f1=fftshift(fftn(ifftshift(volume)))
                f2=fftshift(fftn(ifftshift(ref.volume)))
                curves=fsc(f1,f2,geometries['pilot']['field_A']/64)
                np.savetxt(out/f'{kind}-spacing-{spacing}-fsc.csv',curves,delimiter=',',header='shell,frequency_inverse_A,fsc,count',comments='')
                record['mean_fsc_shells_1_to_6']=float(np.nanmean(curves[(curves[:,0]>=1)&(curves[:,0]<=6),2]))
            output['variants'].append(record)
            (out/'metrics.json').write_text(json.dumps(output,indent=2)+'\n')
            print(dataset,kind,spacing,record,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--count',type=int,default=256)
    parser.add_argument('--datasets',default='10028,10049,10076');args=parser.parse_args()
    for dataset in args.datasets.split(','):run(dataset,args.count)
