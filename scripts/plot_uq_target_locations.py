#!/usr/bin/env python3
"""Physical context for all locked targets; no outcome selection or new fit."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    lockpath=ROOT/'research/uncertainty/pilot-selected-targets-v1/locked-targets.json'
    lock=json.loads(lockpath.read_text());hashes={str(lockpath.relative_to(ROOT)):sha(lockpath)}
    colors=['#0072B2','#D55E00','#009E73','#CC79A7']
    fig,axes=plt.subplots(1,3,figsize=(11.5,4.2),layout='constrained')
    coordinates=[]
    for ax,record in zip(axes,lock['datasets']):
        ds=record['dataset']
        metadata=ROOT/f'results/uncertainty/development/reference-registration-v1/{ds}.json'
        registration=json.loads(metadata.read_text())
        path=ROOT/registration['arrays_file']
        if sha(path)!=registration['arrays_sha256']:raise ValueError('Changed pilot array')
        hashes[str(metadata.relative_to(ROOT))]=sha(metadata)
        hashes[str(path.relative_to(ROOT))]=sha(path)
        np.testing.assert_allclose(registration['field_A'],record['field_A'],rtol=1e-6)
        plane=np.load(path)['pilot'].mean(axis=1)  # [z,x], integrating y.
        field=registration['field_A']
        extent=np.array([-.5,.5,-.5,.5])*field
        limit=float(np.quantile(np.maximum(plane,0),.99))
        ax.imshow(plane,origin='lower',extent=extent,cmap='Greys',vmin=0,vmax=limit,interpolation='nearest')
        for index,feature in enumerate(record['features']):
            x,y,z=feature['center_A'];label=f'R{index+1}' if index<3 else 'C'
            ax.add_patch(Circle((x,z),record['width_A'],fill=False,edgecolor=colors[index],lw=1.6))
            ax.plot(x,z,'o',color=colors[index],ms=3)
            # Fixed offsets separate the broad projected kernels' labels.
            offsets=[(22,16),(-38,14),(22,-20),(-34,-28)]
            ax.annotate(label,(x,z),xytext=offsets[index],textcoords='offset points',
                color=colors[index],weight='bold',fontsize=10,
                arrowprops=dict(arrowstyle='-',color=colors[index],lw=.8),
                bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.))
            coordinates.append(dict(dataset=ds,feature=feature['name'],center_A=[x,y,z],sigma_A=record['width_A']))
        suffix=' (heterogeneous)' if ds=='10076' else ''
        ax.set_title(f'EMPIAR-{ds}{suffix}\nOld independent Gaussian pilot',fontsize=10)
        ax.set_xlabel(r'$x$ ($\AA$)');ax.set_ylabel(r'$z$ ($\AA$)');ax.set_aspect('equal')
        ax.text(.03,.03,f'Band endpoint: {record["planned_radius12_band_endpoint_A"]:.2f} Å',
            transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
    fig.suptitle('Locked averaging features: projected centers and 20 Å standard-deviation contours',fontsize=11)
    out=ROOT/'paper/figures';out.mkdir(exist_ok=True)
    for suffix in ['png','pdf']:fig.savefig(out/f'locked-target-locations.{suffix}',dpi=220)
    plt.close(fig)
    record=dict(complete=True,source_sha256=sha(Path(__file__)),input_hashes=hashes,coordinates=coordinates,
        projection='Mean of the saved old independent Gaussian pilot 64-cube along y; native x,z coordinates and full field. Pilot rendering uses its declared radius-0.35 support.',
        display='Per-map grayscale from zero to the 99th percentile of positive projection values; display clipping only.',
        interpretation='Old pilot reconstruction context, not a truth map or uncertainty map. Circles show projected Gaussian one-SD contours, not hard support or confidence regions.',
        outputs={suffix:sha(out/f'locked-target-locations.{suffix}') for suffix in ['png','pdf']})
    (out/'locked-target-locations.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':main()
