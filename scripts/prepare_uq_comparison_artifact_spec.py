#!/usr/bin/env python3
"""Collect complete comparison evidence, explicitly including numerical failures."""
import hashlib,json
from pathlib import Path
from review_uq_candidate import snapshot_references
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/uncertainty/development'


def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    owner=ROOT/'provenance/uncertainty/review2-comparison-artifact-collection.json'
    specpath=ROOT/'provenance/uncertainty/review2-comparison-artifact-spec-v1.json'
    if owner.exists() or specpath.exists():raise RuntimeError('Preserve previous collection')
    summary=BASE/'continuous-gaussian-review2-summary-v2/summary.json'
    d=json.loads(summary.read_text())
    if not d.get('complete') or d['fits']!=48 or d['conditional_scenarios']!=672:
        raise ValueError('All baseline fits/scenarios required')
    for p,h in d['input_hashes'].items():
        if sha(ROOT/p)!=h:raise ValueError('Changed baseline input '+p)
    files=set()
    for study in ['continuous-gaussian-review2-v1','continuous-gaussian-review2-v2',
                  'continuous-gaussian-review2-summary-v2','phase-split-control-v1',
                  'registered-dictionary-replay-v1','breakdown-radii-review2-v1']:
        for p in (BASE/study).rglob('*'):
            if p.is_file() and p.suffix in ['.npz','.json','.csv']:files.add(p)
    for name in ['continuous-gaussian-warm-start-runtime','continuous-gaussian-diagonal-preflight',
                 'continuous-gaussian-rank8192-preflight','continuous-gaussian-solve-gap-first-case',
                 'cached-quadrature-runtime']:
        files.update((BASE/'audit-regressions').glob(name+'.*'))
    files.update((ROOT/'provenance/uncertainty').glob('continuous-gaussian-v1-stop*.json'))
    for p in list(files):
        if p.suffix=='.json':
            for name in snapshot_references(json.loads(p.read_text())):
                q=ROOT/name
                if sha(q)!=q.stem:raise ValueError('Changed source snapshot')
                files.add(q)
    scope=('This complete collection includes all 48 converged matched continuous-prior fits, all 672 prescribed conditional cases, '
        'the interrupted incomplete v1 with four nonconverged fits, all numerical-only probes, the 40,000-dataset phase control, '
        'all 12 registered dictionary replays, and all 96 breakdown cases. Collection completeness is not success of every scientific/numerical attempt. '
        'The v1 incomplete flags and all unfavorable outcomes remain unchanged. Original experimental inputs are obtained from the public data URLs and prior release dependencies.')
    owner.write_text(json.dumps(dict(complete=True,scope=scope,baseline_summary_sha256=sha(summary),
        member_count_before_collection_record=len(files)),indent=2)+'\n')
    files.add(owner);oh=sha(owner);on=str(owner.relative_to(ROOT))
    spec=dict(version='v0.7.0-dev-comparisons',scope=scope,
        files=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p),owner=on,owner_sha256=oh) for p in sorted(files)],
        required_prior_arrays=['Prior v0.6.2/v0.6.1/v0.6.0 bundles and original EMPIAR/EMDB data described in README and per-study protocols.'],
        reproduction='Use immutable source snapshots and each declared protocol; complete final reports verify source/input/array hashes. The incomplete v1 is preserved historical numerical evidence and is not silently retried.')
    specpath.write_text(json.dumps(spec,indent=2)+'\n');print(len(files),'members',flush=True)


if __name__=='__main__':main()
