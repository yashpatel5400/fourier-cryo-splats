# Post hoc diagnosis of the frozen local-alignment study

1 October 2026 UTC. Motivated by full Fable review 3, R13/R14. This is a
**post hoc reanalysis**, not an additional independent validation study. The
original 600 datasets and their released outputs will not be modified. No
new random observations are generated. All three stacks, both targets, both
alignment templates, all 200 trials and both image controls are retained.

Before calculating the following diagnostics, the analysis plan is:

1. Verify each generator, trial JSON and trial NPZ against its archived hash.
   Recompute the saved centres from the observations, weights and poses.
2. For each of the five distinct estimators (audit and four Gaussian variants),
   report mean error, empirical standard error, RMSE against the fixed pilot,
   standardized error using the weight norm, and paired same-minus-independent
   image shifts. Decompose error into its realized signal contribution and
   projected measurement noise. Only the independent-image noise has the
   claimed conditional standard normal law after conditioning on alignment.
3. Compute actual coverage of noise-only intervals and folded-normal intervals
   at B in {0, .5, 1, 2}, with **unchanged weights**. Measure the exact continuous
   distance between the 24-cell pilot and 64-cell truth, so out-of-class
   reduced-radius results are identifiable. These radii are diagnostic choices
   motivated by the review, not selected for their observed coverage.
4. Add a true-pose control: one deterministic audit and two fixed-pose Gaussian
   fits per target/stack using the original operator, noise, pilot and solver
   settings. Replay both saved images. Retain any solve failures. This is an
   added control on reused data, not part of the original frozen comparison.
5. For every audit fit calculate the continuous class-worst-case bias at the
   **realized true and estimated poses**, without a Taylor expansion. If D
   dephases by the estimated shift, A is the true operator and Ahat the fitted
   operator, set delta = w' (DA-Ahat) rho0. The exact envelope is
   |delta| + B ||ell-A*D*w||. The pose-only envelope is
   |delta| + B ||(A*D*-Ahat*)w||. These use the same continuous cube class.
   Finite Fourier sums are integrated with order-40 Gauss--Legendre quadrature
   and its analytic remainder; roundoff/NUFFT tolerances remain numerical
   diagnostics, not a validated arithmetic enclosure. Independently check
   signs and norms against dense sinc integration on small problems.
6. Compare realized envelopes with the old fixed, first-order-without-remainder,
   and full pose bounds, separately recording calibrated-ball inclusion.
   A realized envelope is not a uniform pose-set bound, a global minimax lower
   bound, or a confidence interval available without true poses. No coverage
   claim follows by substituting it into a data-adaptive procedure.

Individual binomial intervals describe Monte Carlo uncertainty, without a
simultaneous multiple-comparison claim. Mean standard errors use independent
simulation replicates. All rows and numerical discrepancies will be retained.
No redesigned simulation starts before these diagnostics and a working
alignment gate inform a new protocol.
