#!/usr/bin/env python3
"""Compact, all-case presentation of the completed local-refinement study."""
from pathlib import Path
import hashlib
import json
from fourier_splats.uq_end_to_end_summary import METHODS

ROOT=Path(__file__).resolve().parents[1]
LABELS=['Fixed folded','Fixed sum',r'Bounded pose$^\dagger$',r'Mixed 0$^\dagger$',r'Mixed .1$^\dagger$',r'GP $B$, fixed',r'GP $B$, pose',r'GP $B/2$, fixed',r'GP $B/2$, pose']

def main():
    datasets=['10028','10049','10076']; data={}; hashes={}
    for ds in datasets:
        p=ROOT/f'results/uncertainty/development/end-to-end-local-pose-summary-v1/{ds}.json'
        d=json.loads(p.read_text())
        if not d['complete'] or d['attempted_replicates']!=200 or len(d['groups'])!=72:
            raise ValueError('All frozen outcomes required')
        data[ds]=d;hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    groups=[r for d in data.values() for r in d['groups']]
    # Verify the descriptive caption against every procedure/control cell.
    assert all(r['raw_covered']['successes']==200 and r['covered']['successes']==200 for r in groups)
    lines=[r'\begin{table}[t]',r'\centering\scriptsize',
        r'\caption{Local-refinement results: range of median raw half-width/no-data ratios across both targets and templates. All 216 cells have raw and selected coverage 200/200 (individual exact 95\% Monte Carlo interval [.9817, 1]). $^\dagger$Every trial selects the no-data fallback, with relative width 1 and no sign exclusion. Widths are identical for the paired image controls. This tests the declared recipe, not all possible pose-aware estimators.}',
        r'\label{tab:refittingwidths}',r'\begin{tabular}{lrrr}',r'\toprule',r'Procedure & 10028 & 10049 & 10076 \\',r'\midrule']
    records=[]
    for method,label in zip(METHODS,LABELS):
        vals=[]
        for ds in datasets:
            rows=[r for r in data[ds]['groups'] if r['method']==method and r['image_mode']=='independent_image']
            assert len(rows)==4
            for r in rows:
                paired=next(x for x in data[ds]['groups'] if all(x[k]==r[k] for k in ['method','template','target']) and x['image_mode']=='same_image')
                assert paired['raw_relative_half_width']==r['raw_relative_half_width']
                if method in ['deterministic_pose','mixed_common0','mixed_common01']:
                    assert r['fallback']['successes']==200 and r['correct_sign_exclusion']['successes']==0
            widths=[r['raw_relative_half_width']['median'] for r in rows]
            lo,hi=min(widths),max(widths)
            records.append(dict(dataset=ds,method=method,min_median_relative_width=lo,max_median_relative_width=hi))
            vals.append(f'{lo:.0f}--{hi:.0f}' if hi>1 else f'{lo:.3f}--{hi:.3f}')
        lines.append(label+' & '+' & '.join(vals)+r' \\')
    lines += [r'\bottomrule',r'\end{tabular}',r'\end{table}']
    (ROOT/'paper/tables/end-to-end-widths.tex').write_text('\n'.join(lines)+'\n')
    (ROOT/'provenance/uncertainty/end-to-end-main-presentation.json').write_text(json.dumps(dict(input_hashes=hashes,records=records,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Post-outcome all-case presentation only; no fitted or selected interval changes'),indent=2)+'\n')

if __name__=='__main__': main()
