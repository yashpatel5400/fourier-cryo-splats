#!/usr/bin/env python3
"""Retain all thirty comparisons and the unchanged primary gate."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def main():
    out=BASE/'population-materiality-report-v1'
    if out.exists():raise ValueError('Preserve reports')
    entries=[]
    for ds in ['10028','10049','10076']:
        d=BASE/'population-materiality-v1'/ds
        j=json.loads((d/'summary.json').read_text());v=json.loads((d/'independent-check.json').read_text())
        assert j['complete'] and v['complete'];entries.append((ds,j,v))
    out.mkdir();rows=[];table=[];amplitude=[];checks=[]
    for ds,j,v in entries:
        primary=next(g for g in j['grids'] if g['bins_per_axis']==8)['pairs']
        harmonic=[r['amplitude_tangent_harmonic'] for r in primary]
        table.append(f"| {ds} | {max(abs(r['bias']) for r in primary):.8f} | {min(harmonic):.6f}–{max(harmonic):.6f} | {'pass' if j['primary_screen_passed'] else 'fail'} |")
        h=j['haar'];s=h['energy']['mean'];perp=h['amplitude_residual_full']['mean']
        amplitude.append(f"| {ds} | {s:.6f} | {perp:.6f} | {100*(1-perp/s):.2f}% | {h['inverse_energy']['mean']:.6f} | {h['inverse_amplitude_residual_full']['mean']:.6f} |")
        for g in j['grids']:
            for p in g['pairs']:rows.append(dict(dataset=ds,grid=g['bins_per_axis'],**p))
        checks.append(f"- {ds}: all 65,536 information/projection values and every metadata histogram replay; {len(v['expected_scores'])} expected population scores checked by direct density-ratio integration with a separate 160-node quadrature. Maximum score residual {max(abs(x['score_160']) for x in v['expected_scores']):.3g}; maximum amplitude-residual discrepancy {max(v['projected_residual_maximum_errors']):.3g}; twelve direct cell sums across the first, boundary and last training batches differ by at most {max(x['maximum_error'] for x in j['physical_checks']):.3g}.")
    with (out/'all-pairs.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    passed=sum(j['primary_screen_passed'] for _,j,_ in entries)
    manifest=dict(complete=True,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  comparisons=len(rows),stacks_passed=passed,required_stacks=2,gate_passed=passed>=2)
    (out/'summary.json').write_text(json.dumps(manifest,indent=2)+'\n')
    lines=['# Population materiality screen: result','',
        '1 October 2026 UTC. **The primary gate fails on all three stacks.** All thirty prelisted pair/assignment/grid comparisons complete. The frozen criterion requires at least two stacks to have modeled absolute population bias at least .02 and harmonic amplitude-tangent surrogate at most 8.8125, for the same pair on the primary 8-bin-per-axis grid. No sensitivity grid substitutes for a primary outcome.',
        '', '| EMPIAR | Largest absolute primary bias | Primary harmonic surrogate range | Stack gate |',
        '|---|---:|---:|---|',*table,'',
        'The failure is materiality, not the harmonic threshold. The largest primary bias, .003493 on 10076, is about 0.35 percentage points versus the required 2 percentage points. No sensitivity-grid bias reaches .02 either. Under the stopping rule, we do not build the proposed population-feature method on these compact regions. This does not show that pooled likelihood is generally adequate, that latent poses are harmless, or that other conformations have negligible population bias.',
        '', '## What was calculated','',
        'A is each previously fitted normalized map; B deletes its fixed 20-Angstrom peak region. The population is .75. The calculation knows the pose, one CTF and white noise, and uses 220 Fourier frequencies. We compare the two published source halves on all stacks and published filtering versus its complement on 10028/10076, in both assignments to A/B. These metadata groups are not biological states. Their estimated rotations define a piecewise-Haar proxy through fixed 4/8/12-bin ZYZ grids. No actual biological class-specific law or experimental population truth is available in these files.',
        '', 'The exact-in-model population-score root uses one-dimensional Gaussian integration; its odds relation is implicit. The weak-signal approximation and every signed outcome are retained in the CSV. The training draws are reused development data, not an independent confirmation set. Within-cell Monte Carlo errors in each JSON condition on the fixed histograms and omit pose/model/metadata uncertainty.',
        '', '## Amplitude separation diagnostic','',
        '| EMPIAR | Haar mean region energy | Full-map tangent residual | Arithmetic reduction | Haar inverse-energy mean | Haar projected harmonic surrogate |',
        '|---|---:|---:|---:|---:|---:|',*amplitude,'',
        'Real-amplitude projection removes about 5%, 37% and 8% of arithmetic mean separation, respectively. This is a local Gaussian-mean diagnostic. Its harmonic substitution is not the variance theorem for a population estimator with unknown amplitudes; that theorem in the algebra note requires amplitude one and known poses. A small projected harmonic value neither certifies latent-pose identifiability nor provides an uncertainty interval.',
        '', '## Verification','',*checks,'',
        'Four targeted tests compare the formulas with direct density integration, separately optimized expected likelihood and explicit real-coordinate least squares. The later-batch projection checks also address the external reviewer’s concern about orientation regeneration beyond its first sample. No new images, regions, seeds, extra rotations or relaxed thresholds were introduced.',
        '', 'See the [frozen protocol](POPULATION-MATERIALITY-PROTOCOL.md), [algebra audit](POPULATION-SCREEN-ALGEBRA.md) and [unchanged focused consultation](reviews/post-round04-method-consultation/critique.md). This negative development result supplies neither a new uncertainty method nor an ICML acceptance claim.']
    (ROOT/'research/uncertainty/POPULATION-MATERIALITY-RESULTS.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
