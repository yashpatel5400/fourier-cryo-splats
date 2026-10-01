#!/usr/bin/env python3
"""Complete tables for third-order diagnostics; no selected success table."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def main():
    r=json.loads((BASE/'bispectrum-hull-v1/summary.json').read_text())
    p=json.loads((BASE/'bispectrum-adversarial-pose-v1/summary.json').read_text())
    verification=json.loads((ROOT/'provenance/uncertainty/bispectrum-independent-verification.json').read_text())
    assert r['complete'] and p['complete']
    assert verification['complete']
    lines=['# Translation-invariant moment diagnostics: complete outcomes','',
     'This is post-review method development, not experimental calibration. '
     'All 24 declared cells finish in %.3f seconds: twelve exact true-map controls '
     'and twelve removal fits. Every removal fit reaches the 2,500-iteration limit; '
     'none meets the 1e-5 distance-gap criterion. Feasible mixtures and separating '
     'directions bound the finite-hull distance independently of convergence.'%r['seconds'],
     'The moments use the first saved CTF/noise profile and an oracle alternative '
     'averaged over 4,160 orientations. Power plus 512 complex bispectra has '
     '1,244 real coordinates. Amplitude assumptions differ from the earlier '
     'unrestricted-amplitude bilinear e-value model.','',
     '| Stack | Features | Amplitude | Distance lower | Distance upper | Separation after other 10,000 views | After continuous searches | Finite-catalog sufficient n |',
     '|---|---|---|---:|---:|---:|---:|---:|']
    rows=[]
    for c in r['cases']:
        if c['candidate']!='region_removed':continue
        amp='1' if len(c['amplitudes'])==1 else '[.9,1.1]'
        key=c['features']+('_fixed' if amp=='1' else '_range09_11')
        adversary=next(a for a in p['cases'] if a['dataset']==c['dataset'] and a['key']==key)
        remaining=adversary.get('remaining_mean_separation')
        n=c['finite_catalog_cantelli_sufficient_particles']
        text='no direction' if remaining is None else f'{remaining:.6g}'
        lines.append(f"| {c['dataset']} | {c['features']} | {amp} | {c['distance_lower']:.6g} | {c['distance_upper']:.6g} | {c['combined_profiled_separation']:.6g} | {text} | {f'{n:.6g}' if n is not None else 'none'} |")
        rows.append(dict(dataset=c['dataset'],features=c['features'],amplitude=amp,
            distance_lower=c['distance_lower'],distance_upper=c['distance_upper'],status=c['status'],
            finite_union_separation=c['combined_profiled_separation'],adversarial_separation=remaining,
            finite_catalog_sufficient_particles=n,maximum_null_noise_variance=c['maximum_null_noise_variance'],
            alternative_variance=c['alternative_variance']))
    searched=[c for c in p['cases'] if c.get('complete')]
    total=sum(len(c['starts']) for c in searched);successful=sum(x['success'] for c in searched for x in c['starts'])
    maximum=max(c['direct_cell_discrepancy'] for c in searched)
    lines.extend(['',
     '**The last column is not a continuous-pose sample-size guarantee.** It is '
     'the declared two-Cantelli sufficient bound for the enlarged finite rotation '
     'catalog, continuously profiled over the amplitude interval, before the '
     'adversarial searches. It is not a necessary sample size, lower complexity '
     'bound or empirical power measurement. The negative 10028 margins supersede '
     'any tempting favorable interpretation of finite-catalog diagnostics.','',
     '## Continuous orientations and independent calculations','',
     f'All {total} declared searches complete ({successful} optimizers report success) in '
     f"{p['seconds']:.3f} seconds. Four all-zero 10076 directions are explicitly skipped. "
     'Each selected pose is recomputed by the physical cell Fourier operator. '
     f'The eight worst selected views also agree with independent direct cell sums within {maximum:.3g}. '
     'The interpolation guide never serves as a continuous certificate.','',
     'A separate polynomial-moment enumerator, independent of the real-derivative '
     f"implementation, checks all {verification['variance_comparisons']} declared conditional "
     'variance comparisons across eight nonzero contrasts, three means and two '
     f"amplitudes. The maximum absolute discrepancy is {verification['maximum_variance_difference']:.3g}. "
     'All twelve saved removal mixtures and distance brackets also replay. '
     'These checks establish numerical agreement under the assumed noise law, '
     'not its physical calibration.','',
     'Both ranged-amplitude contrasts on 10028 have null orientations whose '
     'expected score exceeds the alternative expectation. The fixed-amplitude '
     '10028 contrasts and all four 10049 contrasts retain positive margins at '
     'the searched poses, which supplies no upper bound over all orientations. '
     'On 10049, the ranged bispectrum margin is about .0296 after search, compared '
     'with .0174 for power only. This is an oracle directional comparison with '
     'different noise variance and nonzero fitting gaps, not demonstrated testing '
     'power or a new reconstruction result. No positive separator was found for '
     '10076; its nonzero upper residuals and fitting limits preclude a claim '
     'that the full third-moment family cannot distinguish the structures.','',
     '## Decision and open requirements','',
     'Do not promote any saved direction to a calibrated continuous-view test. '
     '10049 warrants a genuinely global continuous-orientation bound or an '
     'explicitly narrower justified model before noise/power trials. The first '
     'profile, assumed independent circular Gaussian noise, map normalization '
     'and amplitude interval still require physical validation. The present '
     'Cantelli calculation needs no exponential moment; a cubic Gaussian '
     'statistic cannot simply replace a bilinear term in the earlier e-value. '
     'The exact Hermite variance, convex projection and moment invariants are '
     'classical mathematics, not established novelty. All three full Fable '
     'reviews remain rejections; no fourth full acceptance review occurred.',''])
    (ROOT/'research/uncertainty/paired-power-v1/BISPECTRUM-RESULTS.md').write_text('\n'.join(lines))
    out=BASE/'bispectrum-hull-v1/results.csv'
    with out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print('Wrote complete twelve-removal-case report and CSV')


if __name__=='__main__':main()
