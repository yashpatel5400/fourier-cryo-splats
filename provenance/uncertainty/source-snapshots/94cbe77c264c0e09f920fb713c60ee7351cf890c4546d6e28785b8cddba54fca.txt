#!/usr/bin/env python3
"""Frozen phase-only calibration diagnostic, all paired outcomes retained."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from scipy.stats import norm, t
from fourier_splats.uq_phase_control import cross_power_amplitude_interval
from fourier_splats.uq_end_to_end import binomial_interval
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/'research/uncertainty/PHASE-SPLIT-CONTROL-PROTOCOL.md'
METHODS = ['known_phase_oracle', 'same_image_student', 'independent_image_nominal',
           'independent_image_student', 'invariant_cross_power_bound']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    for path in [Path(__file__), PROTOCOL, ROOT/'research/uncertainty/PHASE-SPLIT-DERIVATION.md',
            ROOT/'src/fourier_splats/uq_phase_control.py', ROOT/'tests/test_phase_control.py']:
        if subprocess.check_output(['git', 'show', f'HEAD:{path.relative_to(ROOT)}'], cwd=ROOT) != path.read_bytes():
            raise ValueError('Commit diagnostic protocol and source before outcomes')
    out = ROOT/'results/uncertainty/development/phase-split-control-v1'
    if out.exists():
        raise RuntimeError('Preserve earlier control')
    out.mkdir(); start = time.perf_counter()
    report = dict(complete=False, planned_replicates_per_case=2000, methods=METHODS, records=[],
        source_snapshot=source_snapshot(ROOT, Path(__file__), [str(PROTOCOL.relative_to(ROOT)),
            'research/uncertainty/PHASE-SPLIT-DERIVATION.md', 'tests/test_phase_control.py']),
        scope='One-frequency phase-alignment diagnostic with equal known Gaussian coordinate noise. Not a full cryo-EM experiment, novel reconstruction method, or experimental calibration.')
    def save():
        (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    save()
    try:
        for si, sigma in enumerate([2., 1., .5, .25]):
            for n in [16, 64, 256, 1024, 4096]:
                tick = time.perf_counter(); estimates = np.empty((2000, 5)); bounds = np.full((2000, 5, 2), np.nan)
                rng = np.random.default_rng(np.random.SeedSequence([261003, si, n]))
                for begin in range(0, 2000, 100):
                    stop = begin+100; shape = (100, n)
                    theta = rng.uniform(-np.pi, np.pi, size=shape)
                    signal = np.exp(1j*theta)
                    a = signal+sigma*(rng.normal(size=shape)+1j*rng.normal(size=shape))
                    b = signal+sigma*(rng.normal(size=shape)+1j*rng.normal(size=shape))
                    phase = np.angle(a)
                    oracle = (b*np.exp(-1j*theta)).real
                    same = abs(a); separate = (b*np.exp(-1j*phase)).real
                    cross = (a.conj()*b).real.mean(axis=1)
                    means = [oracle.mean(axis=1), same.mean(axis=1), separate.mean(axis=1)]
                    halves = [np.full(100, norm.isf(.025)*sigma/np.sqrt(n)),
                              t.isf(.025, n-1)*same.std(axis=1, ddof=1)/np.sqrt(n),
                              np.full(100, norm.isf(.025)*sigma/np.sqrt(n)),
                              t.isf(.025, n-1)*separate.std(axis=1, ddof=1)/np.sqrt(n)]
                    for mi, mean in enumerate([means[0], means[1], means[2], means[2]]):
                        estimates[begin:stop, mi] = mean
                        bounds[begin:stop, mi, 0] = mean-halves[mi]
                        bounds[begin:stop, mi, 1] = mean+halves[mi]
                    estimates[begin:stop, 4] = np.sqrt(np.maximum(0., cross))
                    for offset, value in enumerate(cross):
                        interval = cross_power_amplitude_interval(float(value), n, sigma)
                        if not interval['empty']:
                            bounds[begin+offset, 4] = [interval['lower'], interval['upper']]
                path = out/f'sigma{sigma:g}-n{n}.npz'
                np.savez_compressed(path, estimates=estimates, endpoints=bounds)
                case = dict(sigma=sigma, amplitude_snr=1/sigma, particles=n, replicate_count=2000,
                    arrays_file=str(path.relative_to(ROOT)), arrays_sha256=sha(path), methods=[])
                for mi, method in enumerate(METHODS):
                    lower, upper = bounds[:, mi, 0], bounds[:, mi, 1]
                    finite = np.isfinite(lower)&np.isfinite(upper)
                    covered = finite&(lower <= 1)&(1 <= upper); count = int(covered.sum())
                    width = (upper[finite]-lower[finite])/2
                    case['methods'].append(dict(method=method, covered=count, coverage=count/2000,
                        exact_binomial_95=binomial_interval(count, 2000), empty_count=int((~finite).sum()),
                        mean_estimate=float(estimates[:, mi].mean()), signed_bias=float(estimates[:, mi].mean()-1),
                        sampling_sd=float(estimates[:, mi].std(ddof=1)),
                        finite_width_count=int(finite.sum()), median_half_width=float(np.median(width)) if len(width) else None,
                        mean_half_width=float(width.mean()) if len(width) else None))
                case['seconds'] = time.perf_counter()-tick
                report['records'].append(case); save()
                print(sigma, n, [(r['method'], r['coverage']) for r in case['methods']], flush=True)
        report.update(complete=True, seconds=time.perf_counter()-start); save()
    except Exception as exc:
        report['error'] = repr(exc); save(); raise


if __name__ == '__main__':
    main()
