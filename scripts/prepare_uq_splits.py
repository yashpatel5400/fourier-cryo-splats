#!/usr/bin/env python3
"""Assign whole source micrographs/films to development and inference splits.

The published consensus poses remain conditional inputs; grouping does not
retroactively make those poses independently estimated gold-standard poses.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import requests
import starfile

ROOT=Path(__file__).resolve().parents[1]
SEED=20260929
LABELS=['pilot','tune','inference_half0','inference_half1','test']
EDGES=[.20,.30,.55,.80,1.]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def groups_for(dataset):
    p=ROOT/'background/cryodrgn_empiar'/('empiar'+dataset)/'inputs'
    cs=np.load(next(p.glob('*.cs')))
    selected=np.load(ROOT/'data'/dataset/'indices.npy')
    if dataset=='10076':
        path=ROOT/'data/10076/Frealign9Parameter_0_r1.par'
        url='https://ftp.ebi.ac.uk/empiar/world_availability/10076/data/Frealign9Parameter_0_r1.par'
        if not path.exists():
            r=requests.get(url,timeout=90);r.raise_for_status();path.write_bytes(r.content)
        par=np.loadtxt(path,comments='C')
        assert len(par)==len(cs)
        assert np.array_equal(par[:,0],np.arange(1,len(par)+1))
        source_indices=cs['blob/idx'][selected].astype(int)
        groups=par[source_indices,7].astype(int).astype(str)
        # Published STAR micrograph numbers are sequential particle numbers.
        # FILM in the Frealign file carries actual repeated exposure groups.
        star=starfile.read(p/'Parameters.star')
        assert np.allclose(par[:,8],star['rlnDefocusU'])
        assert np.allclose(par[:,9],star['rlnDefocusV'])
        detail={'source':url,'sha256':sha(path),'field':'Frealign FILM, column 8',
                'all_source_groups':int(len(np.unique(par[:,7]))),
                'identity_validation':'One-based particle index and all defocus U/V values agree with deposited STAR file.'}
    else:
        path=p/('shiny_2sets.star' if dataset=='10028' else 'allimg.star')
        frame=starfile.read(path)
        def key(image):
            index,filename=image.split('@',1)
            parts=Path(filename).parts[-2:] if dataset=='10028' else [Path(filename).name]
            return int(index)-1,'/'.join(parts)
        mapping={key(name):str(group) for name,group in zip(frame['rlnImageName'],frame['rlnMicrographName'])}
        assert len(mapping)==len(frame)
        groups=[]
        for row in cs[selected]:
            filename=row['blob/path'].decode()
            parts=Path(filename).parts[-2:] if dataset=='10028' else [Path(filename).name]
            groups.append(mapping[(int(row['blob/idx']),'/'.join(parts))])
        groups=np.array(groups)
        detail={'source':str(path.relative_to(ROOT)),'sha256':sha(path),
                'field':'rlnMicrographName','all_source_groups':int(frame['rlnMicrographName'].nunique()),
                'identity_validation':'Joined by original source stack path and zero-based image index, not by row order.'}
    return selected,groups,detail


def main():
    out=ROOT/'research/uncertainty/splits';out.mkdir(parents=True,exist_ok=True)
    manifest={'seed':SEED,'policy':'SHA256(seed/dataset/source-group), deterministic uniform threshold assignment',
              'split_fractions':dict(zip(LABELS,np.diff([0]+EDGES))),
              'limitation':'Published full-consensus poses and CTFs are conditional inputs, not independently estimated gold-standard halves.',
              'datasets':{}}
    for dataset in ['10028','10049','10076']:
        selected,groups,detail=groups_for(dataset)
        def assign(g):
            digest=hashlib.sha256(f'{SEED}/{dataset}/{g}'.encode()).digest()
            fraction=int.from_bytes(digest[:8],'big')/2**64
            return LABELS[np.searchsorted(EDGES,fraction)]
        assignments={g:assign(g) for g in np.unique(groups)}
        split=np.array([assignments[g] for g in groups])
        path=out/(dataset+'.csv')
        with path.open('w') as f:
            writer=csv.writer(f,lineterminator='\n');writer.writerow(['output_index','metadata_source_index','source_group','split'])
            writer.writerows(zip(range(len(selected)),selected,groups,split))
        detail.update(selected_particles=len(selected),selected_groups=len(assignments),split_file_sha256=sha(path),
                      counts={label:{'particles':int(np.sum(split==label)),
                                     'groups':len(np.unique(groups[split==label]))} for label in LABELS})
        manifest['datasets'][dataset]=detail
        print(dataset,detail['selected_groups'],detail['counts'],flush=True)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
