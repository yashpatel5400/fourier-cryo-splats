# Reproduce the saved-array diagnosis after full review 3

This is post hoc reanalysis, not a new frozen simulation. The original results
remain in v0.7.0-dev; the revised manuscript and incremental diagnostic bundle
are in v0.7.1-dev. Full Fable review 3 rejects the earlier 26-page manuscript.
The current 25-page revision has not received a fourth full review.

Start from the tagged repository and its Python environment. Download and
verify the three v0.7.0-dev refitting bundles first; they contain the original
calibration generators, all 600 trial arrays/records and their source snapshots.
The v0.7.1-dev increment includes new true-pose weights, every diagnostic CSV,
all summaries, numerical verification and the exact reporting code. It does
not duplicate the three older numerical bundles or the raw movie.

The source at `187d7cc` computed the post hoc replay. It verifies original hashes
and refuses to replace an existing output. To recompute in the same checkout,
use a separate output name, or extract inputs into an isolated checkout without
the increment's already-computed output directory:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  scripts/analyze_uq_refitting_bias.py --output refitting-bias-reanalysis-local
```

The default output `refitting-bias-reanalysis-v1` is what the reporting scripts
read. In an isolated recomputation use that default; then run:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  scripts/summarize_uq_refitting_bias.py
PYTHONPATH=src .venv/bin/python scripts/write_uq_revision4_tables.py
PYTHONPATH=src .venv/bin/python scripts/plot_uq_phase_control_revision4.py
bash scripts/build_paper.sh
```

The Gaussian class table also needs the v0.7.0-dev comparisons bundle and its
stated earlier dependencies. The phase redraw only changes tick formatting of
an existing result. The independent envelope checker uses 36 preselected
higher-order checks and three full dense sinc integrations; it refuses to
replace its existing provenance output:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  scripts/verify_uq_refitting_envelopes.py
```

Individual binomial intervals are not simultaneous across procedures. Reduced
radii below 2 exclude all three complete generators, and true-pose envelopes
are oracle diagnostics rather than deployable intervals. Passing these
replays does not establish scientific novelty, experimental calibration or
conference acceptance. Separate paired-power prototypes are not claims of
this revised manuscript.
