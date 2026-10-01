#!/usr/bin/env python3
"""Post hoc acquisition-level description of saved corner spectra, without CIs."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
OUT = BASE/'acquisition-background-description-v1'


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    if OUT.exists():
        raise ValueError('Preserve descriptive attempts')
    OUT.mkdir(parents=True)
    result = dict(complete=False, scope='Post hoc descriptive acquisition heterogeneity; no pure-noise identification, independence, confidence or new held-out claim', datasets=[])
    for ds in ['10028', '10049', '10076']:
        split = ROOT/'research/uncertainty/splits'/(ds+'.csv')
        array = BASE/'background-spectrum-inventory-v2'/(ds+'-moments.npz')
        with split.open() as f:
            rows = list(csv.DictReader(f))
        labels = np.array([r['source_group'] for r in rows])
        folds = np.array([r['split'] for r in rows])
        groups, inverse, counts = np.unique(labels, return_inverse=True, return_counts=True)
        with np.load(array) as f:
            energy = f['particle_patch_energy'].mean(axis=1)
            raw_power = f['dct_power']
            power = raw_power / np.diag(f['dct_reference_covariance'])[None, :]
        group_energy = np.bincount(inverse, weights=energy)/counts
        group_power = np.stack([np.bincount(inverse, weights=power[:, j])/counts for j in range(63)], axis=1)
        global_mean = energy.mean()
        between = np.sum(counts*(group_energy-global_mean)**2)
        within = np.sum((energy-group_energy[inverse])**2)
        total = np.sum((energy-global_mean)**2)
        np.testing.assert_allclose(between+within, total, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(np.average(group_power, axis=0, weights=counts), power.mean(axis=0), rtol=1e-13)
        weighted = power.mean(axis=0); weighted /= weighted.mean()
        equal_group = group_power.mean(axis=0); equal_group /= equal_group.mean()
        pilot = power[folds=='pilot'].mean(axis=0)
        strata=[]
        for fold in ['pilot','tune','inference_half0','inference_half1','test']:
            mask=folds==fold; p=power[mask].mean(axis=0); ratio=p/pilot; shape=ratio/ratio.mean()
            strata.append(dict(split=fold,particles=int(mask.sum()),source_groups=len(np.unique(labels[mask])),
                energy_mean=float(energy[mask].mean()),relative_to_pilot_coordinate_power=ratio.tolist(),
                pilot_corrected_shape_range=[float(shape.min()),float(shape.max())]))
        record=dict(dataset=ds,source_groups=len(groups),particles=len(rows),
            inputs={str(p.relative_to(ROOT)):sha(p) for p in [split,array]},
            group_mean_energy_quantiles={str(q):float(np.quantile(group_energy,q)) for q in [0,.05,.25,.5,.75,.95,1]},
            descriptive_between_group_energy_fraction=float(between/total),
            anova_identity_absolute_error=float(abs(between+within-total)),
            particle_weighted_shape=weighted.tolist(),equal_group_weighted_shape=equal_group.tolist(),
            equal_group_to_particle_shape_ratio_range=[float((equal_group/weighted).min()),float((equal_group/weighted).max())],splits=strata)
        np.savez_compressed(OUT/(ds+'-group-moments.npz'),source_groups=groups,particle_counts=counts,
            group_mean_energy=group_energy,group_mean_reference_adjusted_dct_power=group_power,
            particle_weighted_shape=weighted,equal_group_weighted_shape=equal_group)
        result['datasets'].append(record)
        print(ds,record['descriptive_between_group_energy_fraction'],record['group_mean_energy_quantiles'],flush=True)
    result['complete']=True
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
