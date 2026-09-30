# Response to the follow-up mathematical audit

The authentic follow-up closes B1–B4 and B6–B8, and confirms the positive-scale
selection argument (B12). It does not assess ICML acceptance. The full round-1
rejection and its unresolved usefulness/calibration objections remain in force.

**B9, legacy fallback reconstruction:** fixed by asserting that the operator's
broadcast radii are constant, then passing scalars to the older nonlinear
projection helper. The archived studies all use uniform radii. Nonuniform
legacy recovery raises an explicit error rather than silently using one radius.
`scripts/validate_uq_legacy_fallback.py` forces the missing-raw-mean/fallback
metadata path on the archived 10049 two-degree fit and reconstructs all three
reference means. Its maximum absolute difference from the original raw means
is **0.0**. The copied metadata is labeled a regression fixture; no scientific
record, image, pose or weight was modified. Evidence is in
`results/uncertainty/development/audit-regressions/legacy-fallback-reconstruction.json`.

**B10, dummy pilot center:** removed the zero default. A missing center is
represented by `None`, allowed only while the affine estimator is selected.
Selecting fallback without an explicit pilot target raises an error; a test
checks that failure as well as the non-fallback branch.

**B11, reporting population:** the summary now scans all result directories for
the audit stage, including sharp, product and custom output folders. It also
enumerates case records in the defined pose-optimization/exchange source grid,
marking unattempted and unfinished fits. At this snapshot there are 22 audit
attempts, all complete, across the variants; the main comparison table retains
its stated ten-row joint/unchanged-remainder scope. The source population has
21 fit records, eleven not yet post-audited. These counts are a dated inventory,
not evidence of calibration or an assertion that all optimization is finished.

After these fixes the full suite passes **89 tests in 7.52 seconds**. Ordinary
floating-point arithmetic remains numerically checked rather than validated
by interval arithmetic. The new tests cannot establish the experimental noise,
pose or density assumptions.
