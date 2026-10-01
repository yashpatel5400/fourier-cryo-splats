#!/usr/bin/env python3
"""All-feature, fixed-weight class/noise sensitivity; no data-driven refit."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.optimize import brentq
from fourier_splats.uq_intervals import bias_aware_half_width_stable
from fourier_splats.uq_joint_bias import joint_density_pose_bias, scaled_pose_radius
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
PROTOCOL = ROOT/'research/uncertainty/BREAKDOWN-RADIUS-PROTOCOL.md'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def increasing_root(function):
    if function(0.) >= 0:
        return dict(value=0., already_nonpositive_margin_at_zero=True)
    hi = 1.
    while function(hi) < 0 and hi < 1e10:
        hi *= 2
    if function(hi) < 0:
        return dict(value=None, lower_bound=hi, already_nonpositive_margin_at_zero=False)
    return dict(value=float(brentq(function, 0., hi, xtol=1e-12, rtol=1e-12)),
        already_nonpositive_margin_at_zero=False)


def main():
    for path in [Path(__file__), PROTOCOL]:
        if subprocess.check_output(['git', 'show', f'HEAD:{path.relative_to(ROOT)}'], cwd=ROOT) != path.read_bytes():
            raise ValueError('Commit reporting rule before calculation')
    out = BASE/'breakdown-radii-review2-v1'
    if out.exists():
        raise RuntimeError('Preserve previous calculation')
    out.mkdir()
    report = dict(complete=False, input_hashes={}, records=[],
        source_snapshot=source_snapshot(ROOT, Path(__file__), [str(PROTOCOL.relative_to(ROOT))]),
        scope='Post-outcome conditional sensitivity of all fixed estimators; no physical class calibration, no new inference, and no choice among the reported procedures.')
    def read(p, expected=None):
        digest = sha(p)
        if expected is not None and digest != expected:
            raise ValueError(f'Changed input: {p}')
        report['input_hashes'][str(p.relative_to(ROOT))] = digest
        data = json.loads(p.read_text())
        if not data.get('complete') or data.get('error'):
            raise ValueError(f'Incomplete input: {p}')
        return data
    def save():
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    try:
        for ds in ['10028', '10049', '10076']:
            raw = read(ROOT/f'results/uncertainty/confirmation/noise-calibration-v1/{ds}.json')
            centered = read(BASE/f'centered-noise-calibration-v1/{ds}.json')
            for feature in raw['features']:
                name = feature['target']; cr = next(r for r in centered['records'] if r['target'] == name)
                for old in feature['records']:
                    source = read(ROOT/old['bias_source'], old['bias_source_sha256'])
                    if old['pose_class'] == 'fixed_pose':
                        fit = source['targets'][0]['fit']
                        slope, intercept = fit['bias']/2, 0.
                    else:
                        original = read(ROOT/source['config']['audit'], source['source_audit_sha256'])
                        B, P = source['density_radius'], source['pilot_norm_bound']
                        h = original['density_bias']/B
                        field = original['pose_polynomial_bias']/(B+P)
                        cross = original['cross']['norm_upper']
                        L = scaled_pose_radius(np.asarray(original['pose_scaling']['group_scales']))
                        remainder = source['selected_cubic_bias']/(B+P)
                        def bias(radius):
                            return joint_density_pose_bias(h, field, cross, L, radius, P, remainder,
                                original['bound']['pilot_cross_upper'])['bias_upper']
                        intercept = bias(0.); slope = bias(1.)-intercept
                        # Verify the affine representation at a third radius too.
                        np.testing.assert_allclose(bias(3.7), intercept+3.7*slope, rtol=1e-12, atol=1e-12)
                    np.testing.assert_allclose(intercept+2*slope, old['bias_upper'], rtol=1e-9, atol=1e-9)
                    cs = next(r for r in cr['intervals'] if r['pose_class'] == old['pose_class'])
                    for noise_label, row, sd in [('raw', old, feature['variance_calibration']['noise_sd_upper']),
                            ('centered', cs, cr['calibration']['centered']['noise_sd_upper'])]:
                        center = row['raw_observed_center']; alpha = raw['alpha_noise']
                        half = bias_aware_half_width_stable(sd, intercept+2*slope, alpha)
                        np.testing.assert_allclose(half, row['raw_half_width'], rtol=1e-9, atol=1e-9)
                        density = increasing_root(lambda B: bias_aware_half_width_stable(sd, intercept+B*slope, alpha)-abs(center))
                        noise = increasing_root(lambda s: bias_aware_half_width_stable(s, intercept+2*slope, alpha)-abs(center))
                        report['records'].append(dict(dataset=ds, target=name, pose_class=old['pose_class'],
                            calibration=noise_label, bias_intercept=intercept, bias_slope=slope,
                            declared_density_radius=2., raw_center=center, noise_sd_upper=sd,
                            raw_half_width=half, raw_half_width_over_abs_center=half/abs(center) if center else None,
                            selected_center=row['interval_center'], selected_half_width=row['interval_half_width'],
                            selected_half_width_over_abs_center=row['interval_half_width']/abs(row['interval_center']) if row['interval_center'] else None,
                            uses_no_data=row['uses_no_data'], declared_excludes_zero=row['excludes_zero'],
                            registered_reference_inside=cs['registered_reference_inside'] if noise_label == 'centered' else None,
                            density_radius_breakdown=density, noise_sd_breakdown=noise,
                            noise_sd_breakdown_over_reported=noise['value']/sd if noise['value'] is not None else None,
                            width_replay_absolute_error=abs(half-row['raw_half_width'])))
                    save()
        if len(report['records']) != 96:
            raise AssertionError('All 96 rows required')
        report['complete'] = True; save()
        rows = []
        for r in report['records']:
            row = {k: v for k, v in r.items() if not isinstance(v, dict)}
            for key in ['density_radius_breakdown', 'noise_sd_breakdown']:
                row[key] = r[key]['value']
                row[key+'_already_fails_at_zero'] = r[key]['already_nonpositive_margin_at_zero']
            rows.append(row)
        with (out/'summary.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        print('All', len(rows), 'rows replayed and reported', flush=True)
    except Exception as exc:
        report['error'] = repr(exc); save(); raise


if __name__ == '__main__':
    main()
