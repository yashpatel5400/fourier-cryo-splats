#!/usr/bin/env python3
"""Write complete method-development tables without selecting favorable cells."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/uncertainty/development'


def read(name):return json.loads((BASE/name/'summary.json').read_text())


def main():
    full=read('paired-covariance-full-frequency-v1');shifted=read('paired-covariance-shifted-v1')
    adversarial=read('paired-covariance-adversarial-10049-v1');corrected=read('paired-statistics-bound-corrections-v1')
    assert all(r['complete'] for r in [full,shifted,adversarial,corrected])
    lines=['# Paired-exposure method gates after review 3','',
      'These are completed oracle information and discretization checks, not experimental '
      'uncertainty validation. The authentic focused Fable audit checks the Gaussian algebra '
      'but identifies missing discrimination bounds, compression loss, and continuous-pose '
      'constraints. The three full ICML reviews still reject. No fourth full review has occurred.','',
      '## Diagonal-power correction and amplitude sensitivity','',
      'All 36 enlarged-cone witnesses replay. Correcting the stationary-root shortfall with '
      'a full-domain concave tangent and a rounding diagnostic changes the 10076 full-removal '
      '128-particle upper from the raw approximately 7.55e-28 to 6.411e-10 at v=1. '
      'Both calculations use ordinary floating point. The near-membership finding survives; '
      'the smaller number is not an arithmetic certificate.','',
      'The 4,160-view feasible mixture mass ranges from .964499 to 1.154280 for half removal '
      'and 1.043575 to 1.146993 for full removal. Half removal therefore has a unit-mass '
      'mixture by convex interpolation; its failure need not use gain inflation. Full removal '
      'needs at least sqrt(1.043575), about 1.0216 amplitude, in this finite catalog. '
      'Those extrema are not continuous-orientation extrema. Diagonal power is defeated '
      'numerically on 10076; 10049 remains open rather than disproved by a small 128-particle expectation.','',
      '## Full-frequency and fresh-view checks','',
      'The 18 compressed cases complete, and every symmetric feature matrix has full '
      'numerical column rank at relative singular-value threshold 1e-8. Rank-32 compression '
      'retains .819/.973/.828 of the true-minus-removal covariance Frobenius norm. '
      'A negative compressed screen is not a negative result for the full matrix family. '
      'The next six-case screen uses all 440 real coordinates and a nonnegative covariance '
      'mixture to bound the best possible bilinear expected growth from above. All three '
      'removal fits converge.','',
      '| Stack | Finite-view achieved growth/particle | Spectral upper | Fraction fresh views violating direction | Achieved growth after repair |',
      '|---|---:|---:|---:|---:|']
    for c in full['cases']:
      if c['candidate']!='region_removed':continue
      g=c['growth']['1.0'];after=c['growth_after_union_repair']['1.0']
      lines.append(f"| {c['dataset']} | {g['expected_log_lower']:.6g} | {g['upper']:.6g} | {c['fraction_fresh_violations_above_1e_minus10']:.4f} | {after['expected_log_lower']:.6g} |")
    lines.extend(['','These initial directions fail on 26.3–44.0% of 10,000 fresh rotations. '
      'The surviving finite union is still not a continuous null. The stored normal '
      'approximation to product-test power is a development calculation, not measured power.','',
      '## Common shifts and adversarial poses','',
      'All twelve larger-catalog cases finish; all twelve NNLS solves report convergence. '
      'The null has 14,160 rotation/shift pairs. Alternative moments integrate the Gaussian '
      'common shift exactly, including both conjugate and unconjugated complex moments. '
      'Each direction is checked on a third set of 10,000 rotations and independent shifts.','',
      '| Stack | Shift SD (64-box pixels) | Catalog lower | Spectral upper | Fresh-union lower | After adversarial repair |',
      '|---|---:|---:|---:|---:|---:|'])
    csvrows=[]
    for c in shifted['cases']:
      sigma=c['shift_sd_pixels'];g=c['fit']['growth']['1.0'];after=c['growth_after_union_repair']['1.0']['expected_log_lower']
      adversary=next((a for a in adversarial['cases'] if c['dataset']=='10049' and a['shift_sd_pixels']==sigma),None)
      last=adversary['growth_after_verified_repair']['1.0']['expected_log_lower'] if adversary else None
      lines.append(f"| {c['dataset']} | {sigma:g} | {g['expected_log_lower']:.6g} | {g['upper']:.6g} | {after:.6g} | {last if last is not None else 'not searched'} |")
      csvrows.append(dict(dataset=c['dataset'],shift_sd_pixels=sigma,catalog_lower=g['expected_log_lower'],
          spectral_upper=g['upper'],fresh_union_lower=after,adversarial_repaired_lower=last))
    lines.extend(['','The most promising stack, 10049, receives 80 declared local searches '
      'over continuous rotations and translations, with every selected pose verified by '
      'the continuous-cell forward operator. All optimizers report success. The four '
      'maximum verified normalized violations are .06710, .04813, .02494 and .007502 '
      'at shift SDs 0,.5,1,2. Interpolation-guide discrepancies are at most 2.95e-5, '
      'far below these violations. The search retains 42 positive verified poses. '
      'The repaired achieved growth becomes zero in the first three cases and '
      '4.132e-5 per particle in the last. That final repair is still only sampled.','',
      'With unrestricted amplitude, a positive null exponent can be amplified without '
      'bound. These are concrete counterexamples to promoting the saved directions '
      'to a continuous-pose e-value. They do not invalidate the conditional Gaussian '
      'proposition, prove the whole bilinear family powerless, or rule out other '
      'cryo-EM validation methods. All bounds, lower values and failures remain available.','',
      '## Decision','',
      'Do not spend a large noise-simulation or GPU budget evaluating these directions '
      'as calibrated continuous-pose tests. Their key null constraint is unmet. '
      'Any continuation needs an actual continuous-orbit bound or a justified narrower '
      'nuisance model, with usefulness demonstrated before experimental calibration. '
      'Mean agreement between paired frames, covariance bounds, particle dependence, '
      'and an experimentally calibrated density class remain open.',''])
    (ROOT/'research/uncertainty/paired-power-v1/METHOD-GATES-RESULTS.md').write_text('\n'.join(lines))
    output=BASE/'paired-method-gate-summary-v1';output.mkdir(exist_ok=True)
    with (output/'shifted-covariance.csv').open('w',newline='') as f:
      writer=csv.DictWriter(f,fieldnames=list(csvrows[0]));writer.writeheader();writer.writerows(csvrows)


if __name__=='__main__':main()
