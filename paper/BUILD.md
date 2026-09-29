# Building the manuscript

The style files are the official ICML 2026 archive, downloaded from:
https://media.icml.cc/Conferences/ICML2026/Styles/icml2026.zip

The manuscript uses the provided `preprint` option, so it does not imply an
ICML submission or acceptance. Authorship/contact details remain unassigned.
The original `initial-draft.tex` records the proposal before real-data fitting.

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
