#!/usr/bin/env python3
"""Algebraic decomposition of the saved test, with unidentified bias labelled."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'matched-event-decomposition-v1'
    if out.exists():raise ValueError('Preserve diagnostic attempts')
    out.mkdir()
    rows=[];hashes={}
    for ds in ['10028','10049','10076']:
        p=BASE/'candidate-fisher-score-v1'/ds/'summary.json'
        hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        j=json.loads(p.read_text());assert j['complete']
        for c in j['cases']:
            cal=c['calibration'];counts=c['counts'];n=j['heldout_views']
            fixed=counts[0][1]/n;alt=counts[2][1]/n;global_grid=max(counts[0])/n
            for b in cal['bounds']:
                upper=cal['joint_mean_interval'][1];group=cal['grouped_mean']
                actual=b['view_variance']
                row=dict(dataset=ds,key=c['key'],kappa=b['kappa'],
                    heldout_null_at_one=fixed,heldout_alternative_quarter=alt,
                    heldout_global_amplitude_three_point_max=global_grid,
                    observed_signal_probability_gap=alt-fixed,
                    grouped_envelope_mean=group,
                    grouped_minus_heldout_null=group-fixed,
                    finite_calibration_mean_increment=upper-group,
                    final_view_cap_increment=actual-upper,
                    final_view_variance_bound=actual,observed_bound_minus_alternative=actual-alt,
                    individual_envelope_mean=cal['individual_mean'],
                    product_estimate=c['decomposition']['product_estimate'],
                    product_sampling_term=c['decomposition']['square_root_term'],
                    product_range_term=c['decomposition']['range_term'],
                    product_distance_subtraction=c['decomposition']['distance_subtraction'],
                    variance_bound=cal['view_variance_upper'],
                    cvar_dkw=b['cvar_dkw'],cvar_split=b['cvar_split'])
                assert abs(fixed+row['grouped_minus_heldout_null']+row['finite_calibration_mean_increment']+row['final_view_cap_increment']-actual)<1e-14
                rows.append(row)
    with (out/'decomposition.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'summary.json').write_text(json.dumps(dict(complete=True,source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),inputs=hashes,rows=len(rows),
        scope='Exact algebraic decomposition of archived floating-point quantities. The grouped-minus-heldout term combines physical nuisance, finite inner replication, cell suprema and different-sample error; none is separately identified.'),indent=2)+'\n')
    lines=['# Saved event-probability decomposition','','1 October 2026 UTC. This is an algebraic diagnosis of the existing Fisher-study test, not a new score, calibration or power run. All 60 score/cap combinations are retained in the CSV. The table shows κ=1.1, which the experimental metadata does not support as a measured cap.','','| EMPIAR | Score | Observed 25% probability gap | Grouped envelope − held-out null | Outer mean increment | View-cap increment | Final bound − observed alternative |','|---|---|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['kappa']==1.1:
            lines.append('| '+r['dataset']+' | '+r['key']+' | '+' | '.join(f'{r[k]:.6g}' for k in ['observed_signal_probability_gap','grouped_minus_heldout_null','finite_calibration_mean_increment','final_view_cap_increment','observed_bound_minus_alternative'])+' |')
    lines.extend(['','The three increments plus the held-out null rate equal the archived final bound to numerical precision. The last increment is the net result of the original minimum and clipping; it is not independently certified physical viewing uncertainty. The positive final gaps explain why the original test is insensitive at this alternative under the specified simulator. Sampling uncertainty in the held-out probabilities remains as reported in the original study.','','**The second column of increments is not an estimate of finite-L bias.** It combines view-dependent amplitude, the cellwise event supremum, maximization after finite noise replication, and differences between calibration/test Monte Carlo draws. The old L=32 versus L=128 runs also used different views and noise, so their difference cannot isolate any one component. The complete CSV separately retains the global three-amplitude grid maximum, individual-noise envelope, mean radius, paired-product sampling/range terms and CVaR comparators. The three-amplitude maximum does not optimize a continuous amplitude law.','','This completes the identifiable saved-data portion of the requested decomposition. A separate matched nested calculation would be needed to estimate inner-replication bias. That calculation cannot repair the demonstrated noise sensitivity, poor retained separation or unsupported experimental viewing model; the moment-test branch remains frozen.','','Inputs and source hashes are in `results/uncertainty/development/matched-event-decomposition-v1/summary.json`.'])
    (ROOT/'research/uncertainty/SAVED-EVENT-DECOMPOSITION.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
