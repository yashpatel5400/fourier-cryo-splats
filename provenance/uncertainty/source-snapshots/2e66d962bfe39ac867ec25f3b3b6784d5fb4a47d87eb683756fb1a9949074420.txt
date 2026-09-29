#!/usr/bin/env python3
"""Common-band real-particle prediction and corrected 3D FSC for neural fits.

All partitions are development partitions. Training losses/bandwidth differ;
this evaluates saved reconstructions on identical particles and Fourier samples.
Exposure-group bootstrap measures prediction variation, not density coverage.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import mrcfile
import torch
import yaml
from cryodrgn.models import load_decoder
from fourier_splats.physics import fft_volume_center
from fourier_splats.fsc import fsc,resolution
from run_experiment import observations,predict

ROOT=Path(__file__).resolve().parents[1]


def neural_prediction(folder,epoch,k,transfer):
    cfg=yaml.safe_load((folder/'config.yaml').read_text())
    cfg['model_args']['zdim']=0  # train_nn configs omit the loader's optional latent dimension.
    weights=folder/f'weights.{epoch}.pkl'
    model,_=load_decoder(cfg,str(weights),torch.device('cpu'));model.eval()
    points=k.reshape(-1,3)/(cfg['lattice_args']['D']-1)
    out=np.empty(len(points),np.complex64)
    with torch.no_grad():
        for start in range(0,len(points),8192):
            a=model.decode(torch.tensor(points[start:start+8192],dtype=torch.float32)).numpy()
            out[start:start+len(a)]=a[:,0]+1j*a[:,1]
    assert cfg['dataset_args']['norm'][0]==0
    out*=cfg['dataset_args']['norm'][1]
    # Check independently saved MRC conversion against direct network evaluation.
    with mrcfile.open(folder/f'reconstruct.{epoch}.mrc') as m:volume=m.data.copy()
    gridfft=fft_volume_center(volume);box=volume.shape[0]
    xyz=np.array([[3,2,1],[-4,1,-2],[0,0,3],[3,-2,0],[4,4,4]],dtype=np.int64)
    with torch.no_grad():
        v=model.decode(torch.tensor(xyz/box,dtype=torch.float32)).numpy()
    direct=(v[:,0]+1j*v[:,1])*cfg['dataset_args']['norm'][1]
    indexed=gridfft[xyz[:,2]+box//2,xyz[:,1]+box//2,xyz[:,0]+box//2]
    err=float(np.linalg.norm(indexed-direct)/max(np.linalg.norm(direct),1e-12))
    if err>1e-4:raise AssertionError(f'Network/map Fourier mismatch: {folder}, {epoch}, {err}')
    return out.reshape(k.shape[:-1])*transfer,gridfft,{'checkpoint_sha256':hashlib.sha256(weights.read_bytes()).hexdigest(),
        'map_direct_relative_error':err,'training_normalization':cfg['dataset_args']['norm']}


def metrics(pred,truth):
    error=np.mean(abs(pred-truth)**2,axis=1).astype(np.float64)
    power=np.mean(abs(truth)**2,axis=1).astype(np.float64)
    corr=float(np.vdot(pred,truth).real/np.sqrt(np.vdot(pred,pred).real*np.vdot(truth,truth).real))
    return {'nmse':float(error.sum()/power.sum()),'correlation':corr},error,power


def group_bootstrap(errors,power,groups,seed):
    labels,inv=np.unique(groups,return_inverse=True);rng=np.random.default_rng(seed)
    e=np.stack([np.bincount(inv,weights=x,minlength=len(labels)) for x in errors])
    p=np.bincount(inv,weights=power,minlength=len(labels))
    choices=rng.integers(len(labels),size=(2000,len(labels)))
    ratios=e[:,choices].sum(axis=-1)/p[choices].sum(axis=-1)
    result={}
    for idx,name in enumerate(['gaussian','voxel']):
        diff=ratios[idx]-ratios[2]
        result[name+'_minus_neural_nmse']={'estimate':float(e[idx].sum()/p.sum()-e[2].sum()/p.sum()),
             'percentile_95':np.quantile(diff,[.025,.975]).tolist()}
    return {'source_groups':len(labels),'resamples':2000,'seed':seed,'paired_differences':result,
            'interpretation':'exposure-group resampling of held-out prediction errors; not density uncertainty'}


def run(args,dataset):
    torch.set_num_threads(4)
    source=ROOT/'results/uncertainty/development/neural-reconstruction'/dataset
    classical=ROOT/'results/uncertainty/development/group-reconstruction'/dataset
    out=ROOT/'results/uncertainty/development/reconstruction-comparison'/dataset;out.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((ROOT/'research/uncertainty/splits'/f'{dataset}.csv').open()))
    groups=np.array([r['source_group'] for r in rows])
    k,obs,c,manifest,_=observations(ROOT/'data'/dataset,64,2048,42,True)
    parts=np.load(classical/'partitions.npz');training=np.concatenate([parts['half0'],parts['half1']])
    scale=np.sqrt(np.mean(abs(obs[training])**2))
    ids=np.concatenate([parts['validation'],parts['test']]);nv=len(parts['validation'])
    kk=k[ids];ct=c[ids];yy=obs[ids];curves=[]
    previous=json.loads((out/'metrics.json').read_text()) if (out/'metrics.json').exists() else {}
    outputs=previous.get('epochs',{})
    controls={};control_maps={};control_errors={}
    for name in ['gaussian','voxel']:
        coefficients=[np.load(classical/f'{name}-half{h}.npz') for h in [0,1]]
        re=sum(x['real'] for x in coefficients)/2;im=sum(x['imag'] for x in coefficients)/2
        p=predict(kk.reshape(-1,3),ct.ravel(),re,im,argparse.Namespace(box=64,sigma=.5,radius=2),name).reshape(yy.shape)*scale
        controls[name]={};control_errors[name]={}
        for split,sl in [('validation',slice(None,nv)),('test',slice(nv,None))]:
            controls[name][split],control_errors[name][split],_=metrics(p[sl],yy[sl])
        control_maps[name]=sum(x['fourier'] for x in coefficients)/2
    for epoch in map(int,args.epochs.split(',')):
        predictions=[];maps=[];checks=[]
        for half in [0,1]:
            p,vol,check=neural_prediction(source/f'half{half}',epoch,kk,ct)
            predictions.append(p);maps.append(vol);checks.append(check)
        pred=sum(predictions)/2
        curve=fsc(*maps,manifest['pixel_size_A']);curves.append(curve)
        np.savetxt(out/f'neural-epoch{epoch}-half-fsc.csv',curve,delimiter=',',header='shell,frequency_inverse_A,fsc,count',comments='')
        record={'epoch':epoch,'checkpoint_checks':checks,'conditional_half_fsc_resolution':resolution(curve),'prediction':{}}
        for split,sl in [('validation',slice(None,nv)),('test',slice(nv,None))]:
            record['prediction'][split],error,power=metrics(pred[sl],yy[sl])
            record[split+'_comparison']=group_bootstrap([control_errors[n][split] for n in ['gaussian','voxel']]+[error],power,groups[ids[sl]],609327)
            np.savez(out/f'epoch{epoch}-{split}-prediction.npz',indices=ids[sl],neural_error=error,
                     gaussian_error=control_errors['gaussian'][split],voxel_error=control_errors['voxel'][split],power=power)
        for name,volume in control_maps.items():
            cross=fsc(sum(maps)/2,volume,manifest['pixel_size_A'])
            np.savetxt(out/f'neural-epoch{epoch}-vs-{name}-fsc.csv',cross,delimiter=',',header='shell,frequency_inverse_A,fsc,count',comments='')
        outputs[str(epoch)]=record
        print(dataset,epoch,record['prediction'],record['conditional_half_fsc_resolution'],flush=True)
    # Replace provisional runner FSCs that used an image-only FFT before this audit.
    last_epoch=int(args.epochs.split(',')[-1])
    np.savetxt(source/'neural-half-fsc.csv',curves[-1],delimiter=',',header='shell,frequency_inverse_A,fsc,count',comments='')
    source_metrics=json.loads((source/'metrics.json').read_text())
    source_metrics['half_fsc_resolution']=outputs[str(last_epoch)]['conditional_half_fsc_resolution']
    source_metrics['fsc_transform']='full centered 3D FFT; cross-checked with direct neural Fourier predictions'
    source_metrics['fsc_epoch']=last_epoch
    (source/'metrics.json').write_text(json.dumps(source_metrics,indent=2)+'\n')
    result={'stage':'development; supplied consensus poses; no true density or gold-standard resolution',
            'dataset':dataset,'config':vars(args),'evaluation_runs':previous.get('evaluation_runs',[])+[vars(args)],
            'common_evaluation':'identical exposure-group partitions and 1410 Fourier pairs/image, radius 30 on 64 grid, window and per-image background normalization',
            'training_difference':'neural train_nn uses radius 32 and per-half normalization; Gaussian/voxel use radius 30 and ridge penalty; this is not an identical training objective',
            'classical_global_scale':float(scale),'controls':controls,'epochs':outputs}
    (out/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--datasets',default='10028,10049,10076');p.add_argument('--epochs',default='5,10,15,20');args=p.parse_args()
    for dataset in args.datasets.split(','):run(args,dataset)
