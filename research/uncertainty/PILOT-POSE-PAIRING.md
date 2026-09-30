# Using the known pilot in the polynomial pose term

This prospective post-audit refinement keeps the same continuous density class,
pose set, weights and spectral event. It does not change the density radius or
declare a new experimental calibration. The first trial is the existing 10049
two-degree full-weight averaged estimator, with sharp cubic constants, before
considering other completed cases. Both successful and unfavorable outcomes
will be retained in a separate output directory.

Write the scalar bias as `<u,h-Fv-e> - <rho0,Fv+e>`, with `||u||<=B`,
`||rho0||<=P`, `||v||<=L`, `||Fv||<=f`, `||e||<=r`, and `c>=||F*h||`.
The previous bound uses `|<rho0,Fv>|<=P f`. The pilot is actually known, so
with `p>=||F*rho0||` we also have `|<rho0,Fv>|<=L p`. Therefore

    |bias| <= B sqrt(||h||²+2 L c+f²) + min(P f,L p) + (B+P)r.

The triangle alternative similarly becomes
`B(||h||+f)+min(P f,L p)+(B+P)r`. Their minimum is valid on the same spectral
event. The proof is Cauchy–Schwarz on the **known** pilot pairing; it is not a
new general statistical theorem. The remaining density–pose interaction and
cubic remainder still cover the entire original continuous class.

For the archived constant-cell pilot, Fourier moments of degree at most two
are analytic. On each cell write `x=c+z`; the Fourier moment factors into
the transform of coefficient-weighted center monomials and uniform-cell
moments. With cell width `1/m`, the latter are

    E exp(2 pi i k z) = sinc(k/m),
    E z exp(2 pi i k z) = sinc'(k/m)/(2 pi i m),
    E z² exp(2 pi i k z) = -sinc''(k/m)/(4 pi² m²).

The normalized constant-cell basis contributes `m^(-3/2)`, and a half-cell
phase shifts the centered FFT grid to actual cell centers. Ten NUFFT transforms
provide all center monomials needed by the existing polynomial pose columns.
No quadrature of the discontinuous pilot, and no restriction of the unknown
density to those cells, enters this calculation. Ordinary floating point,
NUFFT and near-zero sinc-derivative evaluation remain numerical operations,
not interval-arithmetic guarantees.

Independent tests compare all moments with composite quadrature inside every
cell, check the zeroth moment against the existing exact-cell forward model,
compare pilot pose pairings against nonlinear finite differences, and test the
resulting bias inequality on unrelated Hilbert-space instances. These tests
are implementation checks, not evidence of experimental coverage or utility.

The first ten completed post-audits are also checked at their actual 24³ pilot
scale by direct exponential sums, without a NUFFT or centered-grid phase.
Every lifted column is checked for particle positions 0, 64 and 127. The
largest relative pairing discrepancy is 3.4818e-13 and largest absolute
discrepancy is 7.5157e-16. The scalar within-cell formulas are shared, but their
separate composite-quadrature tests are independent. This is a selected-column
numerical diagnostic, not a global rounding-error pad. Evidence is preserved
in `results/uncertainty/development/audit-regressions/pilot-pairing-direct-checks.json`.
