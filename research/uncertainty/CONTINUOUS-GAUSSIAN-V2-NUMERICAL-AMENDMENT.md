# Numerical amendment to the matched continuous Gaussian comparison

1 October 2026 UTC. The original [48-fit scientific protocol](CONTINUOUS-GAUSSIAN-PROTOCOL.md) is unchanged: all twelve locked features, both prior scales, fixed/local-Gaussian poses, both reference frames and every declared stress signal. No scientific outcome is used to choose a feature, prior, pose scenario or density class. No old result is overwritten.

The original rank-1024 implementation repeatedly reaches the 2,000-iteration limit. Its first fixed-pose fit has relative residual 1.038e-5 and a numerical variance-gap upper bound about 9% of its reported variance. Warm starts, a rank-2048 residual-diagonal preconditioner and plan caching alone do not solve convergence. Preserve every probe, including its unfavorable result.

A numerical-only rank-8192 test reduces the unrepresented Gram diagonal to .00244 or less, versus ridge .25. Starting from the first failed fit, its residual reaches 9.50e-12 in three iterations after 315 seconds of setup. This selection uses only linear-system residuals, not new particles, reference coverage or interval utility. It motivates a separate complete rerun.

Version 2 uses the same continuous operator, order-80 quadrature, target scales, initial zero weights, CG tolerance 1e-10 and limit 2,000, with rank-8192 preconditioning and cached FINUFFT plans, one CPU thread. All 48 cases must be attempted and all solver flags retained. Save preconditioner diagnostics, source snapshots, metadata/split hashes and unchanged input/reference hashes. Ordinary floating-point accuracy remains distinct from validated error enclosures.

The slow original run is stopped as a numerically inadequate development implementation, with a separate interruption record preserving its last completed records and hashes. The v1 directory remains incomplete and must not be presented as a completed 48-case study. Version 2 is the required complete comparison. The 600-dataset end-to-end study keeps its own already frozen radius-five solver and is unaffected.
