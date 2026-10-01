#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/uncertainty/development/background-spectrum-inventory-v2'


def main():
    j=json.loads((OUT/'summary.json').read_text());assert j['complete']
    text=['# Descriptive background spectra of the three particle stacks','','1 October 2026 UTC. All 8,192 downloaded particles per stack are included, with four 8×8 corner patches per particle. These are the already processed extracted particles, not new raw movie acquisitions. No particle or corner is removed. The [protocol](BACKGROUND-SPECTRUM-INVENTORY-PROTOCOL.md) was committed before outcomes.','','The known Fourier-crop operation, including removal of the Nyquist row and column, is incorporated in a white-noise reference covariance. Patch means are removed and DCT coefficients transformed using that reference covariance. The real data are **not** empirically whitened. Consequently a white-noise reference would have equal expected normalized powers; the observed ratios below describe departures from that reference.','','| EMPIAR | Patches | Scale | Normalized coordinate-power range | Second-moment eigenvalue range | Off-diagonal norm fraction | Median marginal kurtosis |','|---|---:|---:|---:|---:|---:|---:|']
    fig,axes=plt.subplots(2,3,figsize=(11,6),constrained_layout=True)
    half=[];corners=[]
    for i,d in enumerate(j['datasets']):
        groups=d['groups'];allrow=groups[0];ds=d['dataset']
        text.append(f"| {ds} | {allrow['patches']} | {allrow['reference_whitened_scale']:.6g} | {allrow['diagonal_min']:.3f}–{allrow['diagonal_max']:.3f} | {allrow['eigenvalue_min']:.3f}–{allrow['eigenvalue_max']:.3f} | {allrow['off_diagonal_frobenius_fraction']:.3f} | {allrow['marginal_kurtosis_median']:.3f} |")
        for r in groups:
            if r['corner']=='pooled':
                half.append(f"| {ds} | {r['subset']} | {r['reference_whitened_scale']:.6g} | {r['diagonal_min']:.3f}–{r['diagonal_max']:.3f} | {r['off_diagonal_frobenius_fraction']:.3f} |")
            elif r['subset']=='all':
                corners.append(f"| {ds} | {r['corner']} | {r['reference_whitened_scale']:.6g} | {r['diagonal_min']:.3f}–{r['diagonal_max']:.3f} |")
        power=np.asarray(allrow['dct_power'])/np.asarray(allrow['dct_reference_power']);power/=power.mean()
        tile=np.r_[np.nan,power].reshape(8,8)
        im=axes[0,i].imshow(tile,vmin=0,vmax=3.1,cmap='magma',origin='lower')
        axes[0,i].set(title=f'EMPIAR {ds}',xlabel='DCT x mode',ylabel='DCT y mode')
        for label,style in [('all','-'),('source_half0','--'),('source_half1',':')]:
            r=next(v for v in groups if v['subset']==label and v['corner']=='pooled')
            diagonal=np.diag(np.array(r['second_moment']))/r['reference_whitened_scale']
            axes[1,i].plot(np.arange(1,64),diagonal,ls=style,label=label)
        axes[1,i].axhline(1,color='black',lw=.8);axes[1,i].set(xlabel='Non-DC coordinate',ylabel='Normalized power',ylim=(0,3.2));axes[1,i].grid(alpha=.2)
    fig.colorbar(im,ax=axes[0],shrink=.8,label='Relative DCT power')
    fig.legend(*axes[1,0].get_legend_handles_labels(),loc='outside lower center',ncol=3)
    for ext in ['png','pdf']:fig.savefig(OUT/f'corner-spectra.{ext}',dpi=180)
    plt.close(fig)
    text.extend(['','The nonflat patterns persist in both recorded halves. The largest coordinate-power ratios span approximately 3.8–5.7-fold within a stack after correcting the known preprocessing. These are descriptive departures, not p-values. The column labelled covariance eigenvalues uses the second-moment matrix; the centered covariance is also retained. Mean-vector energy is below 0.01% in each pooled half, so that distinction is small here. Scales are in the downloaded images’ intensity units and are not directly comparable to the unit-noise simulator without its whitening convention.','','## Recorded-half stability','','| EMPIAR | Cohort | Scale | Power range | Off-diagonal fraction |','|---|---|---:|---:|---:|',*half,'','## Separate corners','','Corner order is top-left, top-right, bottom-left, bottom-right in array coordinates.','','| EMPIAR | Corner | Scale | Power range |','|---|---|---:|---:|',*corners,'','## Verification and interpretation','','The explicit 63-frequency Fourier-basis covariance agrees with the analytic reference within 3.78e-15. Thirty-six selected patch transforms reproduce using an independent cosine matrix and linear solve within 1e-12. The first attempt stopped before reading particle outcomes because a 2e-15 reference-check tolerance was smaller than ordinary accumulation error. Its log and original source commit are retained; the second attempt uses 1e-13 and changes no data or procedure.','','Background can contain particle tails, other particles, ice, and preprocessing artifacts. These patches are correlated within each particle. Existing metadata joins recover 229 / 137 / 351 source acquisition groups; the earlier inference that identities were unavailable was incorrect (see the [dated correction](BACKGROUND-SOURCE-GROUP-CORRECTION.md)). This analysis pooled particles and did not use those groups for inference. The recorded halves are not independent noise exposures. The analysis therefore cannot identify pure-noise covariance inside a molecule, establish stationarity, assign valid acquisition-level error bars, or certify the sub-percent noise precision demanded by the earlier moment scores. It does show that treating these archived image coordinates as exactly known white noise needs a measured and validated whitening model. A future likelihood study must retain this limitation instead of transferring unit-noise simulation guarantees directly to experiment.','','Full second moments, centered covariances, marginal moments and source hashes are in `results/uncertainty/development/background-spectrum-inventory-v2/summary.json`; per-particle patch energies and powers are saved in the three NPZ files.'])
    (ROOT/'research/uncertainty/BACKGROUND-SPECTRUM-INVENTORY-RESULTS.md').write_text('\n'.join(text)+'\n')


if __name__=='__main__':main()
