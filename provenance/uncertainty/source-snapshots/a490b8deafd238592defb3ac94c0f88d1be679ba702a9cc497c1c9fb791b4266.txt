# Pilot-only map alignment details

Implementation choices made before evaluating the RELION maps. The fixed
three-stack protocol already specifies 48 principal-axis/hand starts and local
rigid searches. This note makes the numerical choices reproducible; they are
not a new statistical result.

Use the radial Fourier multiplier `exp(-(40 A)^4 |frequency|^4)` and Fourier
resample to a 32-cube. Positive-part mass provides each filtered map's centroid
and principal axes; the correlation objective itself uses the signed maps,
after subtracting their means. Enumerate all 48 signed axis permutations
(24 per handedness). For each, optimize an extra three-component rotation
vector and a translation with Powell's method, at most 300 evaluations/40
iterations. Rotation-vector components lie in [-pi/3,pi/3], translation
components in [-6,6] working-grid pixels. Tolerances are `xtol=1e-4` and
`ftol=1e-6`. Retain the original point if termination gives a worse correlation.
Keep every start, evaluation count and convergence flag; select solely by
filtered pilot correlation. Neither local convergence nor the best correlation
proves a unique or globally optimal frame.

The transform maps output array indices (z,y,x) to input indices. Scale its
offset to the native cube and apply the same trilinear/zero-exterior transform
to both half maps. Native half FSC is computed before this interpolation.
Cross-method/reference FSC after alignment retains interpolation and gauge-
selection caveats. No deposited map participates in this selection.

The independent check renders an asymmetric multi-Gaussian phantom twice from
analytic coordinates, with a known rotation, translation, and each handedness.
It does not generate test inputs with the interpolation being tested. This is
an implementation check, not an experimental alignment validation.

The comparator checkpoints remain the previously frozen neural epochs
60/50/100 for 10028/10049/10076, and the previously frozen Gaussian/voxel half
fits. Their hashes are checked against the existing prediction-study model lock.
No new epoch selection is performed. The native reference fields match the
particle fields; the evaluator verifies this before resampling. Unconverged
last paired refinement maps may be evaluated descriptively, with that status
retained; an initializer without paired refinement maps yields a skipped record.
