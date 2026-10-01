# Author-code population baseline replay

1 October 2026 UTC. This is retrospective reproduction of released summaries,
not a new calibration/coverage experiment, blind test or CAHRA application.
Freeze this protocol and runner before computing population estimates.

Use `aevans1/counting_particles_paper` at commit
`b9099e1b7fb3f03207f94d89b153839bfcd6a7c4` (Evans et al., 2026,
DOI 10.1038/s42003-026-09859-6). Preserve the clone and original outputs.
Its README says the likelihoods were produced by older RECOVAR pipelines;
we do not regenerate or independently validate those likelihoods. The source
LICENSE is GPLv3 while its package metadata says MIT: retain this discrepancy
and do not copy the implementation into our package.

Run the unmodified author's `multiplicative_gradient` and
`deconvolve_assignments` functions on all ten synthetic spike datasets and
all thirteen released experimental log-likelihood arrays (including the
additional noisy-image/ground-truth-likelihood control). Hard assignment and
one-step soft assignment are retained. Use the original synthetic tolerance
1e-8, experimental tolerance 1e-3, maximum 10,000 iterations and default
JAX 32-bit arithmetic. Do not silently relax the tolerance or hide iteration
limits. Save the complete original stdout and recompute the returned-weight
simplex gradient gap in float64. These are current dependencies, not a
reconstruction of the authors' original software environment.

Independently optimize the two-state likelihood in float64 using its monotone
one-dimensional derivative, including boundary optima. Compare the returned
weights and per-image objective values, and compare deconvolution against
its analytic constrained two-state least-squares solution. Preserve all
discrepancies, nonconverged cases and existing published outputs. Synthetic
truth labels may be summarized; the experimental constructed 80/20 reference
is not treated as biological truth. No confidence interval coverage is inferred
from one saved dataset per condition.

The isolated local environment lives under `tmp/counting-baseline-venv`.
Record every input hash, package versions, source identity and runtime. Save
our outputs in a new directory outside the author clone. This baseline replay
tests computational reproducibility, not novelty, nuisance robustness, or
experimental validation of our Gaussian reconstruction method.
