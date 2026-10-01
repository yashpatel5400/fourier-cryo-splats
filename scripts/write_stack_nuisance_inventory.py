#!/usr/bin/env python3
"""Report every declared metadata cohort; plot recorded view histograms."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development/stack-nuisance-inventory-v1'
s=json.loads((BASE/'summary.json').read_text())
check=json.loads((BASE/'independent-check.json').read_text())
assert s['complete'] and check['complete']
lines=['# Recorded imaging nuisances on the three stacks','','1 October 2026 UTC. All three declared metadata inventories finish. This responds to review 4 R18/R19 at the descriptive level; it does not establish the true viewing law, physical amplitude bounds or noise spectra. Source poses are estimates and can reflect symmetry conventions, shared processing and particle selection. The [frozen protocol](STACK-NUISANCE-INVENTORY-PROTOCOL.md) states all grids and subsets.','','## Recorded viewing directions','','For 128 equal-area bins of unoriented plane normals, the peak ratios and cross-half chi-square statistics are well above those compatible with a recorded distribution near uniform. Both are resolution-dependent summaries. A binned maximum is not a continuous-density upper bound. Half-set agreement does not certify pose accuracy or independent processing.','','| Stack/cohort | Particles | Peak ratio | Cross-half chi-square | Half-set TV |','|---|---:|---:|---:|---:|']
for d in s['datasets']:
 for label,r in d['subsets'].items():
  a=r['angular']['normal-8'];lines.append(f"| {d['dataset']} / {label} | {r['particles']:,} | {a['maximum_empirical_ratio']:.4f} | {a['cross_half_chi_square']:.4f} | {a['half_total_variation']:.4f} |")
lines+=['','All 32/128/288 normal-bin and 64/512/1,728 full-rotation-bin counts, including both recorded halves, are retained in the JSON. On all-metadata cohorts, normal-bin peak ratios range from 2.62 to 6.07 (10028), 4.87 to 7.33 (10049), and 2.06 to 3.13 (10076) as bin count increases. The uniform finite-count expectation for plug-in chi-square is recorded separately; the observed deviations are much larger. The published filter changes 10076 appreciably, illustrating that source selection matters.','','## CTFs and recorded scale fields','','| Stack, all metadata | Defocus mean, 5th–95th percentile (micrometres) | Alpha, 5th/50th/95th percentile | Alpha variation associated with view bins |','|---|---|---|---:|']
for d in s['datasets']:
 r=d['subsets']['all_metadata'];f=r['scalar_fields'];df=f['derived/defocus_mean_A']['quantiles'];a=f['alignments3D/alpha']['quantiles'];v=r['association']['alignments3D/alpha']['between_bin_variance_fraction']
 lines.append(f"| {d['dataset']} | {df['0.05']/1e4:.3f}–{df['0.95']/1e4:.3f} | {a['0.05']:.3f} / {a['0.5']:.3f} / {a['0.95']:.3f} | "+('constant recorded field' if v is None else f'{100*v:.2f}%')+' |')
lines+=['','All CTF scale entries equal one and all alignment weight entries equal zero. Those constants are metadata values, not physical confidence statements. Alpha varies on 10028 and 10076, but its exact legacy processing interpretation has not been independently calibrated. The 10076 view association falls from 19.44% to 0.74% after the published filter. No causal claim or valid amplitude interval follows. The constant alpha on 10049 cannot establish absence of physical amplitude variation.','','The nontrivial per-particle defocus spread is incompatible with interpreting the earlier one-CTF simulator as an experimental forward model for the full cohort. These files contain scalar residual/image powers, not frequency-resolved pure-noise covariance estimates. They therefore cannot justify the exactly known white-noise assumption. Existing image-based calibration still has the previously reported signal-contamination and dependence limits.','','## Verification and consequence','','An independent direct-atan2 / `histogramdd` implementation reproduces all 90 full/half histogram count vectors exactly, and verifies every input hash and all scalar finite counts/extrema. The main runner checks every stored rotation against the transposed cryoSPARC rotation-vector convention and checks orthogonality. Neither verifier supplies outward-rounded numerical certification, all-quantile replay, or a statistical guarantee for latent poses.','','The earlier kappa=1.1 Haar comparison is a controlled simulation assumption, not a nuisance level supported by these metadata. This diagnosis supports freezing that branch and comparing latent-pose likelihood information before selecting a replacement. It does not prove that cryo-EM uncertainty inference is impossible or that a more flexible viewing model will work.']
assert len(check['checks'])==90
(ROOT/'research/uncertainty/STACK-NUISANCE-INVENTORY-RESULTS.md').write_text('\n'.join(lines)+'\n')
fig,axs=plt.subplots(2,3,figsize=(10,5.8),layout='constrained')
for col,d in enumerate(s['datasets']):
 for row,label in enumerate(['all_metadata','published_filter']):
  actual=label if label in d['subsets'] else 'all_metadata'
  a=d['subsets'][actual]['angular']['normal-12'];density=np.asarray(a['counts']).reshape(24,12).T*a['bins']/a['particles']
  im=axs[row,col].imshow(density,origin='lower',extent=[0,360,0,1],aspect='auto',vmin=0,vmax=8,cmap='viridis')
  axs[row,col].set_title(d['dataset']+' — '+('all metadata' if row==0 else ('published filter' if actual==label else 'no published filter')))
  axs[row,col].set_xlabel('Normal azimuth (degrees)');axs[row,col].set_ylabel('|normal z|')
fig.colorbar(im,ax=axs,label='Recorded density / uniform',shrink=.84)
fig.suptitle('Estimated viewing directions; not a bound on the true viewing law')
fig.savefig(BASE/'recorded-view-distributions.png',dpi=160)
fig.savefig(BASE/'recorded-view-distributions.pdf')
print('wrote report and figure')
