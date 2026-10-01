# Building the manuscript

The style files are the official ICML 2026 archive, downloaded from:
https://media.icml.cc/Conferences/ICML2026/Styles/icml2026.zip

The manuscript uses the provided `preprint` option, so it does not imply an
ICML submission or acceptance. Authorship/contact details remain unassigned.
The original `initial-draft.tex` records the proposal before real-data fitting.
`reconstruction-v0.1.0.tex` and its PDF preserve the earlier reconstruction
feasibility manuscript. The current `main.tex` is the uncertainty development
draft and explicitly lists unfinished validation.

From `paper/`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The final deliverable is copied to `output/pdf/fourier-cryo-splats.pdf`.
Numbers and figures are generated from saved experiment outputs. Rasterized
pages are inspected before release; compilation alone is not visual QA.
The historical generators `scripts/write_uq_results.py` and
`scripts/summarize_uq_coverage.py` reproduce earlier artifacts. Do not run them
blindly on the current paper: the original dictionary table was superseded by
the all-case registered replay, and the old generator would overwrite that
correction.

The complete restructuring is being developed as `focused-main.tex`, with
`focused-theory.tex`, `focused-experiments.tex` and `focused-survey.tex`. Build it
with the same LaTeX/BibTeX sequence, substituting `focused-main` for `main`.
It remains a working draft until the frozen 600-dataset local-pose study and
48-fit matched continuous-prior comparison have complete, verified summaries.
`write_uq_review2_diagnostics.py` reproduces the completed phase-control figure
and registered-replay/breakdown reports. `summarize_uq_end_to_end_local.py`
refuses final summaries until every prescribed replicate has been attempted.
No incomplete working PDF is promoted to the release deliverable.
