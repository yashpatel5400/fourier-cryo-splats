# Bounded CTF, gain and relative envelope sensitivity

This development extension addresses part of reviewer R8. It neither estimates
the nuisance radii nor changes the completed frozen experiments. The default
study uses declared defocus bounds of 100 and 500 Angstrom, plus a mixed setting
of 100 Angstrom defocus, 2 degrees astigmatism-axis angle, 1 degree phase,
5 percent relative gain and 10 square Angstrom relative B-factor. These are
sensitivity levels, not claimed empirical error distributions.

With amplitude contrast a fixed, the implemented transfer is

    C = sqrt(1-a²) sin(gamma) - a cos(gamma)
      = sin(gamma-arcsin(a)).

At spatial frequency s, an error of at most d in each defocus axis and at most
b radians in its astigmatism angle gives a conservative defocus error

    Delta_df <= d + |df_U-df_V| sin(min(b,pi/2)).

Indeed, change both axes while holding the perturbed angle fixed (a convex
combination with error at most d), then change the nominal angle. The difference
of the two cosine terms is bounded by twice sin(min(b,pi/2)). Thus the phase
error is at most pi*lambda*s²*Delta_df plus the supplied phase-shift radius.
Voltage, spherical aberration and amplitude contrast remain fixed.

The exact range of sine over that interval is obtained from its endpoints and
any included stationary maxima/minima. Multiplying its range by the positive
gain/envelope interval

    [(1-g) exp(-b_B*s²/4), (1+g) exp(b_B*s²/4)]

and taking the largest displacement from the nominal transfer gives e_iq.
No division by the nominal CTF is used, including at a zero. Nonzero e values
arising solely from floating-point CTF evaluation are numerical error, not part
of a validated-arithmetic guarantee.

For any allowed pose, a complex Fourier wave has unit modulus on the unit-volume
density cube. Triangle inequality and Cauchy--Schwarz therefore imply

    |<rho, [A(C,u)-A(C0,u)]* w>|
      <= ||rho||_2 sum_iq |w_iq,complex| e_iq / noise.

Using ||rho||<=B+||pilot|| yields an additive bias term that can be combined with
the existing nominal-CTF pose audit, including simultaneous CTF and pose errors.
This is a conservative uniform bound, not a new general robust-inference
principle. It loses cross-frequency cancellation and shared-field structure.

The completed higher-band experiment applies it to all six saved radius-12
fixed-pose fine-feature estimators, retaining both favorable and unfavorable
outcomes. Coherent nonlinear nuisance perturbations of the independent map
generator provide implementation diagnostics. They are not exhaustive extrema
or experimental calibration. A joint, less conservative acquisition-field
optimization remains a possible methodological improvement; the present
triangle bound must not be described as such an optimization.
