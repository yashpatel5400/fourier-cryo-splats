# Frozen validation of the integrated-moment pose remainder, v2

29 September 2026. Freeze before this study's geometry fits, signal generation,
pose draws or noise draws. This is a new study; continuous-v1 remains unchanged.

## Motivation and prior information

The development 1/2/5-degree sweep made the spatial-maximum cubic remainder
uninformative in all 18 settings. Integrating the known Fourier column over the
unit cube and retaining its actual estimator weights reduced this remainder by
6.79--14.74 times across all 30 broad-target development settings. All six
one-degree settings then had relative widths 0.300--0.532. These outcomes and
partial ongoing continuous-v1 outcomes were known before this freeze. The new
method is a response to those results, not an original prespecified hypothesis.

The new theorem and implementation are recorded in
`research/uncertainty/CONTINUOUS-MOMENT-REMAINDER.md` and
`src/fourier_splats/uq_continuous_moments.py`. It replaces the cubic remainder
only; first/second derivative fields, their quadrature padding, fixed-pose
weight optimization, and the independent pilot remain unchanged.

## Fixed design

Use all three EMPIAR acquisition populations, four 128-particle geometry
replicates each, radius-five Fourier half planes, B=2 and unit-L2 pilots/truths.
Use the same deterministic permutation of the 4096 additional-image metadata
pool as v1, but take positions 512 through 1023 in consecutive blocks of 128.
These particles are disjoint from v1's positions 0 through 511. Exposure groups
can recur, and none of these subsets is a new independent biological specimen.
Published poses/CTFs remain conditional. No new experimental density labels
are introduced.

Prespecify both center and axial-contrast targets at width 0.07 of the physical
field. Use rotation budgets 1 and 2 degrees, with the existing shift budget
0.01/24 of the field. These limits are sensitivity assumptions, not estimates
of experimental alignment accuracy. Geometry seed is 609451 + accession, with
offset 512. Use new signal/pose/noise base seeds 609521/609531/609541 with the
same deterministic formulas as v1. Keep 10000 scalar Gaussian draws per case.

Use Gauss order 40, preconditioner rank 1024, 100 outer steps and 0.5% gap for
fixed-pose weights; use order 32 for polynomial-field pose Grams. Retain every
case, including nonconvergence and no-data fallback. Recheck the continuous
adjoint residual with analytic sinc integrals independently of quadrature.

## Signals, comparators and outcomes

As in v1, use four unit-L2 constant-cell generators (deposited reference,
ellipsoid union, modulated shell, signed multiscale); five scenarios each
(nominal, random boundary, coherent boundary, pose radius exceeded fourfold,
noise scale doubled); and four exact continuous density-ball boundaries
(both signs at nominal and coherent poses). The last two generator scenarios
deliberately violate assumptions and remain separately reported.

Compare six procedures: noise only, 24-cell fixed-pose audit, continuous
fixed-pose audit, continuous pose audit with the spatial-maximum remainder,
the integrated-moment pose audit, and no data. All use identical fixed weights;
each switches both center and width to the independent no-data rule if needed.
No empirical noise scale, support or density radius is learned in this study.

There are 48 audit settings and 1152 signal/scenario records: 768 in class and
384 deliberate violations. Report all exact Gaussian coverage values, paired
Monte Carlo checks and marginal binomial intervals, width/fallback/power
summaries, old/new remainder ratios, fit gaps, feasible bias comparisons and
runtime. Minimum in-class analytic coverage >=0.95 (tolerance 1e-8) is the
correctness criterion. The usefulness comparison retains all one- and two-
degree settings; report improvements and failures equally. Neither conservatism
nor passing this simulation establishes experimental calibration or novelty.

## Integrity

Commit this protocol, driver, implementation and hash manifest before running
the study. A smoke run uses only eight old development particles and 100 noise
draws with the unrelated 900xxx seeds. Archive exact executed source bytes.
Numerical bugs require a dated amendment and retained outputs; scientific
changes require a new labeled study. The simulation has no outcome-based
stopping rule. No independent reviewer has seen or approved this method yet.
