#!/usr/bin/env python3
"""External high-signal coordinate audit; both 90-degree sign controls retained."""
import csv
import hashlib
import json
from pathlib import Path
import mrcfile
import numpy as np
import starfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/uncertainty/development/cahra-inplane-convention-v1'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def centered_quarter_turn(image,k):
    # Even-size DFT coordinates have their origin at N/2. np.rot90 instead
    # rotates about (N-1)/2; the one-pixel circular shift corrects the origin.
    return np.roll(np.rot90(image,k=k),1,axis=0 if k==1 else 1)


def main():
    if OUT.exists():raise ValueError('Preserve all external convention checks')
    OUT.mkdir(parents=True)
    rows=[];results=[]
    for state in ['closed','open']:
        folder=ROOT/'data/cahra-pose-validation'/state
        df=starfile.read(folder/'relion_dataset.star')['particles']
        with mrcfile.open(folder/'000.mrcs',permissive=False) as f:images=f.data.astype(float)
        assert images.shape==(216,256,256)
        assert np.all(df[['rlnOriginXAngst','rlnOriginYAngst']].values==0)
        lookup={tuple(r[c] for c in ['rlnAngleRot','rlnAngleTilt','rlnAnglePsi']):i for i,r in df.iterrows()}
        assert len(lookup)==216
        for i,r in df.iterrows():
            image_index=int(r['rlnImageName'].split('@')[0])-1;assert image_index==i
            key=(r['rlnAngleRot'],r['rlnAngleTilt'],r['rlnAnglePsi']+90)
            if key not in lookup:continue
            j=lookup[key];a=images[i];b=images[j]
            # No fitted amplitude, translation, filtering or per-image selection.
            for k in [-1,1]:
                prediction=centered_quarter_turn(a,k)
                rel=float(np.linalg.norm(prediction-b)/np.linalg.norm(b))
                x=prediction.ravel()-prediction.mean();y=b.ravel()-b.mean()
                corr=float(x@y/np.sqrt((x@x)*(y@y)))
                rows.append(dict(state=state,source=i,target=j,rot=r['rlnAngleRot'],tilt=r['rlnAngleTilt'],psi=r['rlnAnglePsi'],
                    numpy_quarter_turn=k,relative_l2_error=rel,centered_correlation=corr))
        selected=[x for x in rows if x['state']==state];assert len(selected)==216
        results.append(dict(state=state,paired_images=108,inputs={str(p.relative_to(ROOT)):sha(p) for p in [folder/'000.mrcs',folder/'relion_dataset.star']},
            sign_controls=[dict(numpy_quarter_turn=k,
                median_relative_l2_error=float(np.median([v['relative_l2_error'] for v in selected if v['numpy_quarter_turn']==k])),
                maximum_relative_l2_error=max(v['relative_l2_error'] for v in selected if v['numpy_quarter_turn']==k),
                median_centered_correlation=float(np.median([v['centered_correlation'] for v in selected if v['numpy_quarter_turn']==k]))) for k in [-1,1]]))
    with (OUT/'pairs.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result=dict(complete=True,scope='External high-signal 2D in-plane convention description; no template fitting, pose recovery or noisy population claim',
        source_hash=sha(Path(__file__)),datasets=results)
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
