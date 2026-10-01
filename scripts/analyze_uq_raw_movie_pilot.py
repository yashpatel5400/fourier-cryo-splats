#!/usr/bin/env python3
"""Declared one-movie acquisition diagnostics; no reconstruction/noise calibration."""
import hashlib
import json
from pathlib import Path
import subprocess
import time
import mrcfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fourier_splats.uq_movie_diagnostics import correlation, radial_power, block_average
from fourier_splats.uq_provenance import source_snapshot

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/'research/uncertainty/RAW-MOVIE-PILOT-PROTOCOL.md'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    files = [Path(__file__), PROTOCOL, ROOT/'src/fourier_splats/uq_movie_diagnostics.py',
             ROOT/'tests/test_movie_diagnostics.py']
    for p in files:
        if subprocess.check_output(['git', 'show', f'HEAD:{p.relative_to(ROOT)}'], cwd=ROOT) != p.read_bytes():
            raise ValueError('Commit acquisition-diagnostic code before movie analysis')
    record_path = ROOT/'provenance/uncertainty/raw-movie-pilot-v1.json'
    source = json.loads(record_path.read_text())
    if not source.get('complete'):
        raise ValueError('Complete verified movie required')
    path = ROOT/source['path']
    if path.stat().st_size != source['expected_bytes'] or sha(path) != source['sha256']:
        raise ValueError('Movie changed')
    out = ROOT/'results/uncertainty/development/raw-movie-pilot-v1'
    if out.exists():
        raise RuntimeError('Preserve existing movie diagnostics')
    out.mkdir(parents=True); tick = time.perf_counter()
    summary = dict(complete=False, input_hashes={str(path.relative_to(ROOT)): source['sha256'],
        str(record_path.relative_to(ROOT)): sha(record_path)},
        source_snapshot=source_snapshot(ROOT, Path(__file__), [str(p.relative_to(ROOT)) for p in files[1:]]),
        config=dict(patch_size=512, patch_grid=[8, 8], frame_pairs=[[i, i+1] for i in range(0, 16, 2)],
            difference_scale='(second-first)/sqrt(2)', periodogram='Demeaned, unwindowed, all patches',
            display_block_factor=16, display_limits='Pooled 1st/99th percentiles of odd/even sum displays'),
        scope='Acquisition development only. Difference fields contain signal/motion/dose/detector effects; not independent pure-noise samples or experimental density validation.',
        source_identity_limit=source['identity_scope'])
    sp = out/'summary.json'
    def save():
        sp.write_text(json.dumps(summary, indent=2)+'\n')
    save()
    try:
        with mrcfile.mmap(path, mode='r', permissive=False) as m:
            data = m.data
            if data.shape != (16, 4096, 4096) or data.dtype != np.dtype('float32'):
                raise ValueError('Unexpected acquisition shape/type')
            summary['header'] = dict(shape=list(data.shape), mode=int(m.header.mode),
                nsymbt=int(m.header.nsymbt), voxel_size=[float(m.voxel_size[a]) for a in ['x', 'y', 'z']])
            frames, displays = [], []
            for i, frame in enumerate(data):
                x = np.asarray(frame, dtype=np.float64)
                if not np.isfinite(x).all():
                    raise ValueError(f'Nonfinite frame {i}')
                frames.append(dict(frame_zero_based=i, mean=float(x.mean()), sd=float(x.std()),
                    quantile_levels=[0, .01, .25, .5, .75, .99, 1],
                    quantiles=np.quantile(x, [0, .01, .25, .5, .75, .99, 1]).tolist()))
                displays.append(block_average(x, 16))
            summary['frames'] = frames; save()
            spectra, stats, temporal, modes = [], [], [], None
            for py in range(8):
                for px in range(8):
                    f = np.asarray(data[:, py*512:(py+1)*512, px*512:(px+1)*512], dtype=np.float64)
                    diff = (f[1::2]-f[::2])/np.sqrt(2)
                    for pair, x in enumerate(diff):
                        p = radial_power(x); modes = p['modes']; spectra.append(p['shell_energy'])
                        stats.append(dict(patch_y=py, patch_x=px, pair=pair, mean=float(x.mean()),
                            sd=float(x.std()), horizontal_correlation=correlation(x[:, :-1], x[:, 1:]),
                            vertical_correlation=correlation(x[:-1], x[1:]),
                            parseval_absolute_error=abs(float(p['shell_energy'].sum())-p['variance'])))
                    for a in range(8):
                        for b in range(a+1, 8):
                            temporal.append(dict(patch_y=py, patch_x=px, pair_a=a, pair_b=b,
                                correlation=correlation(diff[a], diff[b])))
            displays = np.stack(displays); odd = displays[::2].sum(axis=0); even = displays[1::2].sum(axis=0)
            summary.update(patch_differences=stats, disjoint_pair_correlations=temporal,
                unaligned_sum_display_correlation=correlation(odd, even),
                maximum_parseval_error=max(r['parseval_absolute_error'] for r in stats))
            ap = out/'diagnostics.npz'
            np.savez_compressed(ap, frame_displays=displays, unaligned_odd_sum_display=odd,
                unaligned_even_sum_display=even, shell_energy=np.stack(spectra), shell_mode_count=modes)
            summary['arrays_sha256'] = sha(ap)
        fig, ax = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
        limits = np.quantile(np.r_[odd.ravel(), even.ravel()], [.01, .99])
        for a, x, label in zip(ax[0, :2], [odd, even], ['Unaligned odd-frame sum', 'Unaligned even-frame sum']):
            a.imshow(x, cmap='gray', vmin=limits[0], vmax=limits[1]); a.set_title(label); a.set_axis_off()
        ax[0, 2].plot(range(1, 17), [r['mean'] for r in frames], 'o-'); ax[0, 2].set(xlabel='Frame', ylabel='Mean pixel value', title='All acquisition frames')
        colors = plt.cm.viridis(np.linspace(0, 1, 8))
        power = np.stack(spectra).reshape(64, 8, -1)
        freq = np.arange(power.shape[-1])/512
        for i, color in enumerate(colors):
            mean = power[:, i].mean(axis=0)/modes
            ax[1, 0].semilogy(freq[1:], mean[1:], color=color, label=f'{2*i+1}–{2*i+2}')
        ax[1, 0].set(xlabel='Radial frequency (cycles/pixel)', ylabel='Power per Fourier mode', title='Disjoint frame differences'); ax[1, 0].legend(ncol=2, fontsize=7)
        for key, label in [('horizontal_correlation', 'Horizontal'), ('vertical_correlation', 'Vertical')]:
            ax[1, 1].hist([r[key] for r in stats if r[key] is not None], bins=30, alpha=.5, label=label)
        ax[1, 1].set(xlabel='Neighbor correlation', ylabel='Patch/pair count'); ax[1, 1].legend()
        ax[1, 2].hist([r['correlation'] for r in temporal if r['correlation'] is not None], bins=40)
        ax[1, 2].set(xlabel='Correlation between disjoint differences', ylabel='Patch/pair count')
        fig.suptitle('One raw movie: acquisition diagnostics, not a noise-independence test', fontsize=12)
        for suffix in ['pdf', 'png']:
            fig.savefig(out/f'acquisition-diagnostics.{suffix}', dpi=150)
        plt.close(fig)
        summary.update(complete=True, seconds=time.perf_counter()-tick,
            figure_hashes={name: sha(out/name) for name in ['acquisition-diagnostics.pdf', 'acquisition-diagnostics.png']})
        save(); print('complete', summary['seconds'], 'seconds;', len(stats), 'patch differences;', len(temporal), 'cross-pair comparisons', flush=True)
    except Exception as exc:
        summary.update(error=repr(exc), seconds=time.perf_counter()-tick); save(); raise


if __name__ == '__main__':
    main()
