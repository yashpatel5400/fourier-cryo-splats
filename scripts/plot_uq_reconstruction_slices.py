#!/usr/bin/env python3
"""Show all locked reconstructions in three fixed planes; no view selection."""
import hashlib
import json
from pathlib import Path
import mrcfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.physics import fft_volume_center,volume_from_fourier

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'
METHODS=['reference','gaussian','voxel','neural','relion']
LABELS=['Registered reference','Gaussian, supplied pose','Voxel, supplied pose','Neural, supplied pose','RELION, fitted pose']


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    hashes={};cases=[]
    lockpath=ROOT/'research/uncertainty/confirmation/prediction-v1/locked-models.json'
    lock=json.loads(lockpath.read_text());hashes[str(lockpath.relative_to(ROOT))]=sha(lockpath)
    axis=np.arange(64)-32
    radius2=sum(x*x for x in np.meshgrid(axis,axis,axis,indexing='ij'))
    mask=(radius2<=30**2)&(radius2>0)
    for ds in ['10028','10049','10076']:
        stem=ROOT/f'paper/figures/reconstruction-slices-{ds}'
        if any(stem.with_suffix(s).exists() for s in ['.pdf','.png','.json']):
            raise RuntimeError('Preserve existing visualization')
        # Verify exactly the inputs used for the already reported FSCs.
        for method in METHODS[1:]:
            p=BASE/f'reference-registered-comparison-v1/{ds}-{method}.json'
            record=json.loads(p.read_text())
            if not record.get('scientific_run_complete'):raise ValueError('Incomplete reconstruction comparison')
            hashes[str(p.relative_to(ROOT))]=sha(p)
            for name,h in record['input_hashes'].items():
                if sha(ROOT/name)!=h:raise ValueError('Changed source: '+name)
                hashes[name]=h
        reg=json.loads((BASE/f'reference-registration-v1/{ds}.json').read_text())
        values={'reference':np.load(ROOT/reg['arrays_file'])['reference_registered']}
        for method in METHODS[1:]:
            if method=='relion':
                value=np.load(BASE/f'relion-evaluation-v2/{ds}/metrics-maps.npz')['relion']
            else:
                halves=[]
                for half in [0,1]:
                    if method=='neural':
                        p=BASE/f'neural-reconstruction/{ds}/half{half}/reconstruct.{lock["epochs"][ds]}.mrc'
                        with mrcfile.open(p) as m:x=m.data.astype(float)
                    else:
                        x=volume_from_fourier(np.load(BASE/f'group-reconstruction/{ds}/{method}-half{half}.npz')['fourier'])
                    halves.append(x)
                value=sum(halves)/2
            values[method]=value
        volumes=[];scales={}
        for method in METHODS:
            v=np.asarray(values[method],float)
            if v.shape!=(64,64,64) or not np.isfinite(v).all():raise ValueError('Invalid volume')
            v=volume_from_fourier(fft_volume_center(v)*mask)
            scale=float(np.sqrt(np.mean(v*v)))
            if scale<=0:raise ValueError('Empty display volume')
            scales[method]=scale;volumes.append(v/scale)
        limit=float(np.quantile(np.abs(np.stack(volumes)),.995))
        fig,axes=plt.subplots(3,5,figsize=(7.2,4.8),constrained_layout=True)
        extent=[-reg['field_A']/2,reg['field_A']/2]*2
        for row in range(3):
            for col,v in enumerate(volumes):
                ax=axes[row,col];im=ax.imshow(np.take(v,32,axis=row),origin='lower',
                    cmap='RdBu_r',vmin=-limit,vmax=limit,extent=extent,interpolation='nearest')
                ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(LABELS[col].replace(', ','\n'),fontsize=7)
                if col==0:ax.set_ylabel(f'Central plane, axis {row}',fontsize=7)
        colorbar=fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.65,pad=.01)
        colorbar.set_label('Zero-mean density / own RMS (display only)',fontsize=7)
        colorbar.ax.tick_params(labelsize=7)
        fig.suptitle(f'EMPIAR-{ds}: all three central planes; common radius-30 display band',fontsize=9)
        for suffix in ['.pdf','.png']:fig.savefig(stem.with_suffix(suffix),dpi=180)
        plt.close(fig)
        record=dict(complete=True,dataset=ds,input_hashes=dict(hashes),
            plotting_source_sha256=sha(Path(__file__)),display=dict(grid=64,radius=30,
                dc_removed=True,per_map_rms=scales,pooled_absolute_quantile=.995,color_limit=limit,
                planes=[{'axis':a,'index':32} for a in range(3)],field_A=reg['field_A']),
            outputs={s:sha(stem.with_suffix(s)) for s in ['.pdf','.png']},
            scope='Illustration of the existing frozen comparisons, not a new fit or numerical accuracy metric. Per-map display normalization removes amplitude differences. Reference registration is pilot-derived; the third reference represents only Class A.')
        stem.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n');cases.append(ds)
    print('Rendered all three fixed planes for reference and four methods on',cases)


if __name__=='__main__':main()
