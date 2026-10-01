# Enlarged viewing-mixture diagnostic

All 36 declared cells completed in 52 seconds. The original Fourier templates
reproduce to at most 5.1e-13 absolute discrepancy. The stored nonnegative
coefficients, Fourier arrays, rotations and source hashes permit replay without
rerunning the linear program. These are numerical witnesses, not validated
interval-arithmetic equalities.

At 4,160 candidate views, the full-removal EMPIAR-10076 candidate reproduces
the original map's uniform-64-view mean power with maximum scaled residual
1.25e-14. Its mixture mass is 1.092598159; thus probabilities proportional to
the coefficients and a common amplitude of the square root of that mass match
the power to numerical precision. The half-removal residual is 1.12e-14 with
mass 1.061743333. The elementary cone lemma in the protocol implies zero
positive expected log growth for the *diagonal linear cross-power family*
when the cone equality is exact. Numerical separable dual values here are
below 8e-28 for 128 particles at unit variance bound.

| Stack | Half-removal scaled residual | Full-removal scaled residual | Full-removal 128-particle log-growth upper diagnostic, v=1 |
|---|---:|---:|---:|
| 10028 | .002384 | .008816 | .13834 |
| 10049 | .172224 | .288448 | 1.97399 |
| 10076 | 1.12e-14 | 1.25e-14 | 7.55e-28 |

The other two stacks do not have equality witnesses in this catalog. Their
displayed upper values use the particular minimax-residual mixture, not an
optimized dual. In particular, the earlier 64-view optimized dual for 10049
is tighter; enlarging the null cannot improve the best achievable growth.
Every intermediate catalog and all variance bounds remain in `summary.json`.
The true-map controls have numerical equality witnesses on all three stacks.

This is a reason to stop development of the diagonal statistic as the principal
method. It does not prove that full particle-image distributions agree,
that a statistical test can never reject, or that reconstruction is impossible.
An expectation bound is not a rejection-power calculation. Arbitrary viewing
mixtures and unrestricted amplitude are material parts of this failure.

A new feasibility check may retain off-diagonal Fourier correlations, with
the added cost of modeling common translations. It must distinguish a sampled
view relaxation from continuous-pose validity and avoid presenting a moment
metric or standard Gaussian moment formula as new general theory.
