# Final-weight tolerance and guide diagnostic

Declared 30 September 2026 in response to focused Fable findings W1/W4/W5,
while the single cubic optimization is running. Preserve its original procedure,
weights, final spectral probe, intervals and source snapshot without alteration.
After that run completes, evaluate the same order-80 nominal Gram at FINUFFT
relative tolerances 1e-12 and 1e-14. Record the change in its quadratic form and
Gram action. Show the effect of adding the observed changes as heuristic density
and dual-defect sensitivity margins; do not label them validated bounds.

At the same final weights, compare one cubic forward/adjoint pair at the two
tolerances using independent diagnostic seed 650129, normalized Gaussian column
and spatial vectors and the same order-64 quadrature. These checks cannot bound
all operator directions or validate the randomized spectral event in floating
point. The original scalar cancellation guard is roundoff-only; it does not
cover NUFFT approximation, which remains an explicit numerical limitation.

Compare the selected evaluation's largest squared Ritz singular value with the
final certificate's Rayleigh lower and probabilistic upper. No extra probability
claim or probe selection is made. Report eta, entropy and roundoff contributions
that can obstruct a tight lower bound. They are not asserted to decompose the
whole optimization gap. Do not retune or select a different estimator from these
checks. Numerical failures are retained.
