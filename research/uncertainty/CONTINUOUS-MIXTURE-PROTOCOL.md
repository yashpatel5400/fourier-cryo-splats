# First full-orientation likelihood-bracket computation

30 September 2026. Declare before generating these observations or fitting
these continuous mixtures. This follows the discrete-view screen but does
not reuse its observations or extend its detection claim by relabeling a grid.

For each of EMPIAR-10028, 10049 and 10076, take the same 128 acquisition CTFs
and radius-twelve half-plane frequencies as its prior known-map information
diagnostic. Use the previously fitted Fourier Gaussian pilot as the specified
map, normalized by its analytic whole-space density L2 Gram norm. Do not
clip, mask or convert that map to voxel cells. The selected Gaussian map is
both the fixed candidate and the generator in this compute-only diagnostic.
It is not a claim that the experimental stack is homogeneous or Gaussian.

Generate 128 independent Haar rotations and Gaussian noise with seed
`910001 + int(dataset)` and the supplied scale from the information diagnostic.
Translations are zero, CTFs exact and whitening isotropic. No new experimental
pixels are accessed. Simulated orientations remain hidden from the optimizer.

Use the full Euler cover [0,2pi] x [0,pi] x [0,2pi], initially 4 x 2 x 4 boxes.
Use the combined Taylor/coordinate likelihood envelopes in the published
derivation. Refit the finite support mixture every 256 splits, with at most
100 EM iterations and its recorded gap. Allow 8,192 splits and 900 seconds
per stack (the time limit is checked between completed cover updates/refits).
The global continuous target gap is one log-likelihood unit. A limit is an
unresolved likelihood bracket, not a converged solution. The whole batch has
at most three such stage budgets, plus bounded setup and final checks.

After the optimization, also evaluate a uniform mixture over the 128 true
latent rotations as an oracle feasible likelihood diagnostic. Do not feed
those rotations to the optimizer or choose runs from that comparison. The
continuous upper bound must contain this likelihood as well as its own best
finite support lower value, up to reported floating-point tolerances.

Retain every stack and its stopping status, likelihood-bracket trajectory,
Euler bisection history, final and best covering partitions, exact center
kernels, mixture anchors/weights, input map and geometry, all generated
observations/rotations, and source hashes. The bound can be replayed from
the retained best cover without rerunning the optimizer. A separate global
numerical enclosure is not established by a successful replay.

This study measures denominator tightness and computational cost. It has no
predictive numerator, no e-value or rejection outcome, and no claim of
experimental density coverage. Before advancing to a structural test, useful
continuous brackets, an independently learned normalized predictor, and
nuisance/model-error validation remain necessary. If the bracket is too wide,
retain that outcome and diagnose the dominant relaxation rather than present
finite-grid likelihoods as the continuous maximum.
