# Building the diagnostic manuscript

The authoritative entry point is `paper/main.tex`, which inputs
`paper/diagnostic-main.tex` and then `paper/diagnostic-appendix.tex`. The old
`focused-main.tex` is the reviewed historical source, not the current entry.
The v0.7.9-dev source/QA archive contains the required figure, tables, bibliography
and ICML 2026 style files. TeX Live 2025 was used on the Mac.

From the repository root, with `pdflatex` and `bibtex` on PATH:

```sh
bash scripts/build_paper.sh
```

This compiles the entry point, resolves the bibliography and cross-references,
and copies the result to `output/pdf/fourier-cryo-splats.pdf`. The checked PDF
has sixteen pages. A rebuild can differ bytewise because of TeX timestamps;
the manifest identifies the exact published file. There are no overfull boxes
or unresolved citations/references in the checked build. Two harmless hyperref
warnings concern a title line break and the style's empty affiliation anchor.

The three tables are generated from recorded results, not hand-entered:

```sh
python scripts/write_diagnostic_paper_tables.py
```

Their input summaries and independent-check records are included in this
archive. Regenerating the numerical results requires the v0.7.8/v0.7.7 archives
and their dependencies; see `REPRODUCE-POST-REVIEW4-GATES.md`. The manuscript
rewrite itself introduces no new numerical run or selected favorable cohort.

The QA record includes the published PDF hash, exact source hashes, page count,
build checks and all sixteen page-render hashes. Every rendered page was
visually inspected after correcting a long equation and code identifier. The
proof review clarified realification of complex observations and distinguished
the ACG angular scale parameter from its actual angular standard deviation.
These checks concern correctness and presentation, not novelty or acceptance.
