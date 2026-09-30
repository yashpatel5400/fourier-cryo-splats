# Completed local Gaussian pose baseline

The protocol and source were published in commit 9df4aac before fitting.
All 48 target/prior/pose-model solves completed across the three stacks, with
384 analytic reference scenarios. The latter are conditional calculations,
not independent experiments. All 24 fixed-pose control widths exactly repeat
the original weak-prior outputs. The maximum variance-identity discrepancy is
2.923e-13. Setup plus fitting takes 8.45/4.58/4.71 s per stack; cumulative
process peak resident memory reaches 1.27 GB. Exact hashes, Jacobians and
weight arrays are recorded in `pilot-selected-fourier-pose-v1`.

Marginalizing the supplied local pose prior widens target intervals by only
0.0066–0.2875 percent. For coordinate prior SD one, continuous-reference
coverage ranges from 0.9861 to 0.9999 across prescribed perturbations; the
minimum sign power remains 0.04276. Both broad priors, all targets and all
scenarios are retained. The summary separates matched prior-predictive
calibration from fixed-reference probabilities.

Small inflation describes this fixed-pilot Jacobian model. It must not be
interpreted as showing that nonlinear pose uncertainty is negligible. In the
one-degree checks, linearizing the trilinear pilot already leaves errors of
0.71–2.41 percent of signal norm, or 14.1–27.3 percent of the pose-induced
change. At two degrees, signal-relative errors reach 2.00–6.71 percent.
Density/pose products are also omitted, and no image-derived pose confidence
radius or independently calibrated experimental noise law is supplied.

This is a classical Gaussian marginalization control, not a new statistical
method, an ab-initio reconstruction, or a reproduction of Rangan et al.'s
volumetric Hessian pipeline. Favorable reference behavior is retained without
claiming uniform coverage over the continuous density/pose class.
