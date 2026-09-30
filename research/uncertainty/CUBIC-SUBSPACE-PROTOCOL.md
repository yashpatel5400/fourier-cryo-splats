# Reduced convex cubic weight design

Development protocol, 30 September 2026. The original cubic fit has a 0.940
full-space surrogate gap and low power; the separate coordinate follow-up is
still running. This probe addresses optimizer behavior, not experimental noise,
pose or density-class calibration. It does not replace either prior fit.

## Design and classical construction

Use the same 10049 first pilot-selected region, 128 particles, radius twelve,
20 A Gaussian target, density radius two, pilot norm one, one-degree/0.5 A joint
pose balls, fixed noise and inherited cubic column scales. All data and these
earlier outcomes are development inputs. Form thirteen candidate weight columns:
the original fixed-pose estimator restricted to each of twelve unit-width
frequency shells, plus the completed original cubic optimized estimator.
Orthonormalize this span with a declared relative SVD threshold of 1e-12. Keep
the singular values and projection residual. Do not use the running coordinate
fit's weights or its outcomes to choose this span.

For weights w=Q a, compute the full cross terms of the density Gram, pilot
polynomial norm and each derivative-remainder Gram in this small coordinate
space. Bound each complex-coordinate amplitude by the sum over basis columns
of |a_j| times that column's amplitude. This triangle inequality is used only
for the quadrature-integration and embedding-residual pads, not the principal
density or derivative terms. Each reduced Gram uses recorded ordinary-floating-
point PSD padding; this is not validated interval arithmetic.

The remaining expensive term is the spectral norm of the linear cubic field
F(Q a). For any unit vectors u,v, its norm is at least |v' F(Q a) u|. Therefore
each queried singular pair gives two valid linear supporting constraints in a.
Solve a small second-order-cone master with these constraints, then query the
new weights and add three directions. This is a classical spectral supporting-
plane method, not a new statistical or convex-optimization principle. Primary
historical source: Kelley (1960), *The Cutting-Plane Method for Solving Convex
Programs*, DOI https://doi.org/10.1137/0108053 (publisher bibliographic record
checked; no claim to have read the full historical paper).

In exact arithmetic, the master relaxes the reduced majorant problem and each
new direction preserves that relaxation. In implementation, a conic primal
objective and approximate Ritz value are numerical diagnostics, not certified
lower/upper values. Any approximate final weights remain eligible for the
original independent coverage audit. A subspace result is not a certificate of
full-space optimality, and changing the fitting space does not shrink the
inference density class.

## Bounded experiment

Initialize at the completed original cubic estimator. Use the unsmoothed norm
objective; the previous final audits already evaluate that unsmoothed bound.
Use density quadrature 80, pose quadrature 64 and the smaller 7,040-dimensional
pose-column Gram for three Ritz directions (tolerance 1e-4, subspace 17, at most
150 eigensolver iterations). Use at most twenty spectral queries. Stop earlier
only when the numerical restricted guide/master gap is at most 0.001. Keep
every query, conic status, gap, coefficient vector and rejected/nonselected
candidate. Conic solver: Clarabel, absolute/relative/feasibility tolerances 1e-8,
at most 500 iterations. Partial or failed eigensolves terminate and preserve the
attempt; do not invent missing leading modes.

Optimizer seed 710011; independent final-certificate seed 710001. Check the
existing result seed inventory before creating the case. The final upper audit
uses four Gaussian probes and forty power iterations at numerical failure
budget 1e-6/12. Total error budget is 0.05/12, as in the original conditional
comparison. Retain the original full-space dual diagnostic, cross-term audit,
direct and PSD quartic remainders, and the three already specified nonlinear
reference scenarios. No reference value selects the estimator.

Run on two CPU threads after the coordinate follow-up terminates to avoid
adding another concurrent heavy fit. This scheduling condition does not select
on its scientific outcome. Report runtime, peak resident memory, subspace size,
fresh spectral upper, full-space surrogate gap, interval width and sign power.
The empirical criterion is improved audited precision and power. Numerical
convergence, a smaller guide value or a passing small test alone is insufficient.

The source and this protocol must be committed before empirical execution.
All prior frozen files remain unchanged. The original and new intervals are
separate alternatives; no unadjusted data-dependent minimum across attempts
is presented as a simultaneous experimental confidence statement.
