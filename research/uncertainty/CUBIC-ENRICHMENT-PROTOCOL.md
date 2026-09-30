# Adaptive subspace enrichment for the unchanged cubic design problem

30 September 2026. Declared after the original, coordinate-metric and fixed
subspace fits, registered-reference checks, and centered-noise development.
This is a numerical optimization study, not new experimental calibration.

Keep the already selected 10049 region 1, 20 A target SD, 128 radius-12
particles, one-degree/0.5 A joint pose class, density radius two, pilot norm
one, original positive cubic scales, and supplied simulation noise. Start
from the original optimized weights and the same thirteen-dimensional span
used by the completed subspace study. Preserve all earlier fits, including
their large full-space gaps and unfavorable experimental intervals.

The prior restricted guide gap was small while the full-space audited gap
remained .940. Enlarge the span using the original unsmoothed objective's
full-coordinate subgradient and a fixed nominal-Gram-preconditioned version
of that direction. Orthonormalize twice against the existing span and retain
only directions with relative orthogonal norm above 1e-10. The preconditioner
ridge is z times the initial majorized density-residual norm divided by B
times the initial weight norm. It is fixed before enrichment; the preconditioner
only guides optimization, not coverage. No reference or experimental image
value chooses a direction.

First evaluate the inherited point, then allow three enrichment rounds, with
at most six spectral evaluations and restricted conic solves per round.
Always evaluate at least one conic proposal in a new span. Use the existing
Ritz support oracle (three modes, tolerance 1e-4, subspace 17, maxiter 150)
and CLARABEL master. A restricted numerical guide gap of .001 can end a
round but is not a full-space convergence certificate. Selection across
different bases uses the original unsmoothed guide objective, not the
basis-dependent integration-pad relaxation. Retain every evaluated point
and every improving checkpoint. A nonconverged guide remains an admissible
candidate for a final bound; it is not an optimized-solution claim.

Design seed 952011 and fresh final certificate seed 952001 must not have
appeared previously in development result seed fields. Allow 7,200 seconds
for design, checked at evaluation boundaries, followed by the existing four-
probe/forty-step independent spectral audit at delta=1e-6/12. Use order 80
for the nominal density Gram and 64 for pose fields. A process alarm at
10,800 seconds bounds the whole attempt subject to native-call interrupt
latency; preserve a timeout as failure, never as a completed certificate.
Use two FINUFFT threads in the isolated CPU environment. Retain peak memory,
all guide/dual gaps, and the original final joint-density/pose bias refinement.

Before the empirical run, compare a small enriched solve to a separately
assembled full dense conic problem and verify preservation/enlargement of
the span. The initial CLARABEL dense reference attempt failed and is retained;
the established independent SCS reference succeeds. This is not an empirical
case change. Neither passing tests nor an improved guide value guarantees
improvement in the final audited interval.

After selection, report both original and registered reference-generator
checks for nominal/coherent/random-boundary poses with their old seed and
normalization. Neither reference participates in design. No new real-image
interval or observed-data weight selection is part of this protocol.
