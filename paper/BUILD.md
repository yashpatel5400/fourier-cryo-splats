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

The authoritative `main.tex` entry point inputs `focused-main.tex`, with
`focused-theory.tex`, `focused-experiments.tex` and `focused-survey.tex`.
The frozen 600-dataset local-pose study and 48-fit continuous-prior comparison
are complete. All unfavorable outcomes and no-data fallbacks are reported.
`paper/uncertainty-v0.6.2-dev/` preserves the historical 55-page checkpoint;
its complete source dependencies remain in the v0.6.2-dev Git tag.
`write_uq_review2_diagnostics.py` reproduces the phase-control figure and
registered-replay/breakdown reports. `summarize_uq_end_to_end_local.py`
requires every prescribed attempt. `write_uq_end_to_end_figures.py` reports
all cells; `summarize_uq_refitting_main.py` creates the compact all-case
main-text width table. These reporting scripts do not refit the study.
