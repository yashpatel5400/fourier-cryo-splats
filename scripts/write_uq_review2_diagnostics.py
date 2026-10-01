#!/usr/bin/env python3
"""Render completed reporting diagnostics from saved, hash-verified results."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'results/uncertainty/development'
NOTES = ROOT/'research/uncertainty'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    hashes = {}
    def read(path):
        data = json.loads(path.read_text())
        if not data.get('complete') or data.get('error'):
            raise ValueError(f'Incomplete result {path}')
        hashes[str(path.relative_to(ROOT))] = sha(path)
        return data
    lines = ['# Registered replay of all original dictionary examples', '',
        '1 October 2026 UTC. All twelve estimators exactly reproduce the archived original widths, biases and analytic coverages (maximum stored replay discrepancy zero). Original-versus-two-stage reference resampling differs by at most 2.15e-7 relatively. Registration is the previously frozen pilot-only transform; no estimator is refitted to improve coverage. These are conditional finite-voxel simulations, not experimental truth labels.', '',
        '| Stack | Target | Width / field | Original coverage | Registered coverage | Registered audited coverage |',
        '| --- | --- | ---: | ---: | ---: | ---: |']
    table = []
    for ds in ['10028', '10049', '10076']:
        path = BASE/f'registered-dictionary-replay-v1/{ds}.json'; data = read(path)
        ap = path.with_name(f'{ds}-arrays.npz')
        if sha(ap) != data['arrays_sha256']:
            raise ValueError('Replay arrays changed')
        hashes[str(ap.relative_to(ROOT))] = sha(ap)
        for row in data['targets']:
            def coverage(frame, method):
                return next(c['analytic_coverage'] for c in row['checks'] if c['frame'] == frame and c['method'] == method)
            old, new, audit = coverage('original', 'restricted_dictionary'), coverage('registered', 'restricted_dictionary'), coverage('registered', 'audited_dictionary')
            lines.append(f"| {ds} | {row['target']} | {row['width_fraction_field']:.2f} | {old:.6g} | {new:.6g} | {audit:.6g} |")
            table.append(dict(dataset=ds, target=row['target'], width_fraction_field=row['width_fraction_field'],
                original_coverage=old, registered_coverage=new, registered_audited_coverage=audit,
                historical_replay_max_error=row['historical_replay_max_error']))
    lines += ['', 'The highlighted 10028 central .07 example changes from 4.18037e-5 to 3.73035e-14 coverage. Other targets improve or worsen; none are omitted. The ambient-audited intervals cover these two specific generators with probability numerically one, which is conservative and does not establish that an experimental density belongs to the class. The 10076 generator remains one Class A map, not heterogeneous consensus truth.', '',
        'Declared protocol: [REGISTERED-DICTIONARY-PROTOCOL.md](REGISTERED-DICTIONARY-PROTOCOL.md). Saved weights and both reference generators accompany the three geometry records.']
    (NOTES/'REGISTERED-DICTIONARY-RESULTS.md').write_text('\n'.join(lines)+'\n')
    with (BASE/'registered-dictionary-replay-v1/summary.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0])); writer.writeheader(); writer.writerows(table)

    breakdown = read(BASE/'breakdown-radii-review2-v1/summary.json')
    lines = ['# Density and noise breakdowns of the existing experimental intervals', '',
        '1 October 2026 UTC. All 96 fixed-weight combinations replay the archived B=2 raw half-widths within the declared tolerance. The calculation varies the assumed density radius while holding the observed center, pilot, poses, image weights and noise calibration fixed. It does not choose a new radius or validate the physical class.', '',
        'The table gives the centered-noise procedure at fixed poses. B* is the boundary where the raw data-based interval first includes zero as B increases. A zero B* means it already includes zero at B=0. The original declared radius remains B=2. Width/center uses the raw interval, avoiding the different no-data center.', '',
        '| Stack | Feature | Half-width / absolute center | B* | Critical SD / reported SD |',
        '| --- | --- | ---: | ---: | ---: |']
    for row in breakdown['records']:
        if row['calibration'] == 'centered' and row['pose_class'] == 'fixed_pose':
            lines.append(f"| {row['dataset']} | {row['target']} | {row['raw_half_width_over_abs_center']:.3f} | {row['density_radius_breakdown']['value']:.3f} | {row['noise_sd_breakdown_over_reported']:.3f} |")
    lines += ['', 'At the one-degree pose class, every feature has a zero noise-SD breakdown: its recorded bias bound already exceeds the observed center, even with measurement SD zero. At B=0, neither homogeneous stack has a one-degree exclusion under either calibration; centered 10076 regions 1 and 3 have only tiny density-radius thresholds (.054 and .016). This diagnoses the existing fixed-weight audit, not intrinsic impossibility or a physical error radius. The separately optimized cubic estimators are not part of these twelve frozen estimators.', '',
        'The complete CSV retains raw and centered calibration, fixed/shift/one-/two-degree classes, selected and raw widths, fallbacks, replay errors and thresholds. Large B* values can reflect a small fixed-pose adjoint residual; they do not validate unknown-pose inference. The five centered 10076 registered-reference disagreements remain unchanged.', '',
        'Declared protocol: [BREAKDOWN-RADIUS-PROTOCOL.md](BREAKDOWN-RADIUS-PROTOCOL.md).']
    (NOTES/'BREAKDOWN-RADIUS-RESULTS.md').write_text('\n'.join(lines)+'\n')

    phase = read(BASE/'phase-split-control-v1/summary.json')
    for case in phase['records']:
        ap = ROOT/case['arrays_file']
        if sha(ap) != case['arrays_sha256']:
            raise ValueError('Phase arrays changed')
        hashes[str(ap.relative_to(ROOT))] = sha(ap)
    labels = ['Known phase', 'Same-image aligned', 'Independent aligned: nominal SD',
              'Independent aligned: sample SD', 'Invariant cross-power bound']
    colors = ['#555555', '#d55e00', '#0072b2', '#56b4e9', '#009e73']
    fig, axes = plt.subplots(2, 4, figsize=(9, 4.8), sharex=True)
    for col, sigma in enumerate([2., 1., .5, .25]):
        cases = [r for r in phase['records'] if r['sigma'] == sigma]
        ns = [r['particles'] for r in cases]
        for mi, (label, color) in enumerate(zip(labels, colors)):
            records = [r['methods'][mi] for r in cases]
            means = np.array([r['coverage'] for r in records])
            bounds = np.array([r['exact_binomial_95'] for r in records])
            axes[0, col].errorbar(ns, means, yerr=np.vstack([means-bounds[:, 0], bounds[:, 1]-means]),
                label=label, color=color, marker='o', ms=3, lw=1, capsize=1)
            axes[1, col].plot(ns, [r['median_half_width'] for r in records], color=color, marker='o', ms=3, lw=1)
        axes[0, col].axhline(.95, color='black', linestyle=':', lw=.8)
        axes[0, col].set_title(f'Amplitude / noise = {1/sigma:g}', fontsize=9)
        axes[0, col].set_ylim(-.035, 1.04)
        axes[1, col].set_yscale('log'); axes[1, col].set_xlabel('Independent particles', fontsize=8)
        for ax in axes[:, col]:
            ax.set_xscale('log', base=4); ax.set_xticks([16, 256, 4096], ['16', '256', '4096'])
            ax.tick_params(labelsize=7); ax.grid(alpha=.2)
    axes[0, 0].set_ylabel('Coverage of true amplitude', fontsize=8)
    axes[1, 0].set_ylabel('Median half-width', fontsize=8)
    handles, names = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, names, loc='lower center', ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=[0, .17, 1, 1])
    fig.savefig(ROOT/'paper/figures/phase-split-control.pdf')
    fig.savefig(ROOT/'paper/figures/phase-split-control.png', dpi=160)
    plt.close(fig)
    lines = ['# Completed phase-alignment independence control', '',
        '1 October 2026 UTC. Twenty declared cases, 2,000 independently simulated datasets per case, five paired procedures, and all 40,000 sets of interval endpoints are retained. This is a one-frequency phase model with known Gaussian coordinate noise, not a cryo-EM reconstruction benchmark. Numerical tests check the product-normal MGF/moments and confidence-set inversion.', '',
        '| Amplitude / noise | Particles | Oracle coverage | Same-image Student | Independent nominal | Independent Student | Invariant bound |',
        '| ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for case in phase['records']:
        rows = case['methods']
        lines.append(f"| {case['amplitude_snr']:g} | {case['particles']} | "+' | '.join(f"{r['coverage']:.4f}" for r in rows)+' |')
    lines += ['', 'At 4,096 particles, every naive aligned interval has zero observed coverage in all four SNR cases, including the independent-image intervals with sample variance. The oracle remains near .95 and the invariant finite-sample bound is conservative. These outcomes demonstrate the stated attenuation/noise-selection counterexample; they do not establish a new inverse method or a universal failure of independent data splitting.', '',
        'The invariant bound uses a different, phase-invariant estimand before taking a square root, and requires the known equal-variance Gaussian model. It does not calibrate experimental noise, CTFs or physical heterogeneity. Empty sets and their counts are retained in the JSON; binomial intervals and widths are available for every case.', '',
        'Protocol: [PHASE-SPLIT-CONTROL-PROTOCOL.md](PHASE-SPLIT-CONTROL-PROTOCOL.md). Proof: [PHASE-SPLIT-DERIVATION.md](PHASE-SPLIT-DERIVATION.md).']
    (NOTES/'PHASE-SPLIT-CONTROL-RESULTS.md').write_text('\n'.join(lines)+'\n')
    (ROOT/'provenance/uncertainty/review2-reporting-diagnostics.json').write_text(json.dumps(dict(
        complete=True, script_sha256=sha(Path(__file__)), input_hashes=hashes,
        figures={str(p.relative_to(ROOT)): sha(p) for p in [ROOT/'paper/figures/phase-split-control.pdf', ROOT/'paper/figures/phase-split-control.png']}), indent=2)+'\n')


if __name__ == '__main__':
    main()
