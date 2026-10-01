#!/usr/bin/env python3
"""Package the incremental post-review-3 replay and verify every archive member."""
from pathlib import Path
import hashlib,json,subprocess,tarfile
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 files=set()
 for directory in ['results/uncertainty/development/refitting-bias-reanalysis-v1','results/uncertainty/development/refitting-bias-summary-v1']:
  files.update(p for p in (ROOT/directory).rglob('*') if p.is_file())
 for name in ['scripts/analyze_uq_refitting_bias.py','scripts/summarize_uq_refitting_bias.py',
  'scripts/verify_uq_refitting_envelopes.py','scripts/write_uq_revision4_tables.py','scripts/plot_uq_phase_control_revision4.py',
  'scripts/verify_uq_release_artifacts.py','scripts/package_uq_revision4.py','scripts/build_paper.sh',
  'src/fourier_splats/uq_refitting_diagnostics.py','tests/test_refitting_diagnostics.py',
  'research/uncertainty/REFITTING-REANALYSIS-PROTOCOL.md','research/uncertainty/REFITTING-BIAS-REANALYSIS-RESULTS.md',
  'research/uncertainty/REPRODUCE-REVISION4.md','research/uncertainty/REVIEW3-SOURCE-CORRECTIONS.md',
  'research/uncertainty/reviews/round-03/review.md','research/uncertainty/reviews/response-to-round-03-development.md',
  'provenance/uncertainty/refitting-realized-envelope-verification.json','provenance/uncertainty/revision4-paper-verification.json',
  'logs/uncertainty/revision4-focused-tests.log','logs/uncertainty/refitting-envelope-verification.log']:
  files.add(ROOT/name)
 # All manuscript source/assets are in the Git tag; include its current inputs.
 for name in ['paper/main.tex','paper/focused-main.tex','paper/focused-theory.tex','paper/focused-experiments.tex',
  'paper/focused-survey.tex','paper/figures/refitting-alignment-bias.pdf','paper/figures/refitting-alignment-bias.png',
  'paper/tables/refitting-bias-detail.tex','paper/tables/continuous-gaussian-class-coverage.tex']:
  files.add(ROOT/name)
 path=ROOT/'output/artifacts/v0.7.1-dev-refitting-diagnosis.tar.gz'
 manifest=path.with_name(path.name[:-7]+'-manifest.json')
 if path.exists() or manifest.exists():raise ValueError('Preserve existing package')
 rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
 with tarfile.open(path,'w:gz') as archive:
  for p in sorted(files):archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
 expected={r['path']:r for r in rows};seen=set()
 with tarfile.open(path,'r|gz') as archive:
  for m in archive:
   assert m.isfile() and m.name in expected and m.name not in seen
   with archive.extractfile(m) as f:h=hashlib.file_digest(f,'sha256').hexdigest()
   assert h==expected[m.name]['sha256'] and m.size==expected[m.name]['bytes'];seen.add(m.name)
 assert seen==set(expected)
 record=dict(complete=True,stream_verification_complete=True,bytes=path.stat().st_size,sha256=sha(path),
  source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),files=rows,
  dependencies=['Tagged repository source/environment','v0.7.0-dev refitting bundles for all three stacks',
    'v0.7.0-dev comparison bundle and its stated older dependencies for contextual tables'],
  scope='Incremental post hoc replay artifacts, not new simulation evidence or a new independent acceptance verdict.')
 manifest.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='files'}));print('members',len(rows))
if __name__=='__main__':main()
