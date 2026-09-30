# Original cubic-weight optimization outcome

Completed 30 September 2026 under the separately committed single-case
protocol (2123d33). EMPIAR-10049, pilot region 1, Gaussian sigma 20 A,
128 particles/radius 12, one degree and 0.5 A in the declared joint pose ball.
This is a development result with supplied class/noise/nuisance assumptions.

The 30-iteration cap is reached with `optimizer_success=false`. All 36
objective evaluations and improving weight checkpoints remain archived,
including rejected large-objective trials. The selected approximate objective
is 2.26735; it is not an upper certificate. The fresh independent spectral
audit gives eigenvalue upper 0.0651904 versus Rayleigh lower 0.0597610.
The audited unsmoothed triangle objective is 2.32211, with continuous dual
lower 0.138869 and relative gap 0.940197. This is not near convergence.

The final joint-refined half-width is 2.18595 versus no-data 12.15870, relative
width 0.179785. The preceding fixed-weight cubic relative width was 0.252312.
The three declared approximate-reference scenarios have minimum correct-sign
probability 0.00653666, despite conditional coverage rounded to one. Width
improvement does not establish useful inference. The original quadratic,
fixed cubic and optimized cubic procedures remain separate alternatives;
their unadjusted minimum is not a simultaneous confidence procedure.

Wall time is 7,191.39 seconds, including 5,907.88 seconds of optimization and
1,265.02 seconds of spectral certification; peak resident memory is
2,593,570,816 bytes. Source and array hashes are in the immutable case record.

The separately declared final-weight checks complete in 22.43 seconds. Changing
NUFFT tolerance from 1e-12 to 1e-14 changes the nominal Gram action by
7.31e-13 relatively and its quadratic form by 1.43e-12, 0.170 times the existing
roundoff guard. One cubic forward/adjoint pair changes by 4.04e-13/5.64e-13
relatively. The fresh spectral upper/Ritz ratio is 1.09033 and Rayleigh/Ritz
ratio 0.999521. These are empirical sensitivity checks, not a validated global
operator-error bound or a replacement certificate.

The gap exceeds the prespecified 0.005 trigger, so the already published
coordinate-metric follow-up runs from the original fixed-pose weights with
the same iteration budget and a fresh certification seed. Its outcome is
pending. This record does not claim that a different coordinate system will
improve convergence or scientific utility.
