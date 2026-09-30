# Bounded joint density/pose cubic design

30 September 2026. Declared after the original, coordinate and reduced cubic
fits, registered-reference and centered/projected-noise outcomes, and while
the separate adaptive triangle-objective enrichment is still running. This
study changes the design objective, not the density/pose/noise assumptions.
It is post-outcome numerical development, not fresh confirmation.

Keep the already selected 10049 region 1, 20 A Gaussian target SD, 128
radius-12 particles, one-degree/0.5 A joint pose class, density radius two,
pilot norm one, original cubic scales and supplied simulation noise. Use the
same thirteen-dimensional span of the twelve original fixed-weight shells
and completed original cubic estimator. Initial weights are the original
cubic weights, regardless of the outcome of the concurrent enrichment.

Use the joint robust-norm guide described in JOINT-TRUST-REGION-THEORY.md.
Construct nominal fields and the pointwise Gaussian target at order 64 only
for design. Each oracle uses 24 fully reorthogonalized Krylov columns, starting
from the pose cross vector and a fixed independent random vector. Solve the
small trust-region problem with explicit hard-case handling; compare its
witness to all prior cuts. Add the selected witness as a quadratic-norm cone
constraint. The restricted CLARABEL master keeps the original noise,
known-pilot and quartic-remainder norm terms, with the same basis-dependent
remainder majorant as the earlier reduced study. Initial zero-pose cut retained.

Allow twelve evaluations and 5,400 design seconds checked at evaluation
boundaries, with a .001 numerical restricted guide-gap stopping rule only
after at least one proposal is evaluated. Select by the joint guide with all
past cut values retained, not a reference/pixel criterion. Save every improving
checkpoint and all evaluations. Such a guide gap is not a continuous or
full-space convergence certificate. Independent small SDP and analytic tests
must pass before this empirical run. Failed attempts are retained.

Use new design seed 953011 and final spectral seed 953001, inventoried against
all development JSON seed fields before starting. After selection, compute
the unchanged order-80 continuous density bound, order-64 cross/quartic bounds,
and a fresh four-probe/forty-step spectral upper event at delta=1e-6/12.
The residual-controlled joint audit uses up to 40 Krylov columns, a fixed
61-point log-shift grid from -20 to 10, and bounded local scalar refinement.
Every selected scalar upper remains valid on the same spectral event in
real arithmetic; no extra weight selection after the audit is permitted.
Use the minimum of the new and old joint upper bounds on that event.

Two FINUFFT threads, isolated CPU runtime, 10,800-second process alarm
(subject to native interrupt latency). Wait for the active enrichment fit to
terminate before starting this heavy job. Retain its outcome separately.
Report all timing/memory, complete or failed status, all old/new bound terms,
and nominal/coherent/random-boundary generator checks in both original and
pilot-registered frames, with the original reference-check seed. No reference
enters design. No real experimental interval is implied by this protocol;
that application, if performed later, must retain its calibration conditions.
