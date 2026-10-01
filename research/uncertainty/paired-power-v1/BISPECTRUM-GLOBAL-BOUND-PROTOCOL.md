# Analytic global-orientation envelope: feasibility calculation

1 October 2026 UTC. Post-outcome development, declared before computing the
new molecular curvature constants. Positive sampled margins alone do not
justify uncertainty claims. This gate measures the resolution required by
one explicit deterministic certificate, without launching a large grid or
new statistical trials. It is classical derivative/covering analysis, not a
novel recovery or uncertainty theorem.

## Uniform curvature and a finite certificate

For fixed physical density rho, write M_p = integral |rho(x)| |x|^p dx.
Its transferred Fourier coordinate is m_j(R)=t_j integral rho(x)
exp(-2 pi i q_j R x) dx. Along any unit-speed rotational geodesic, put
w_j=2 pi |q_j|. Direct differentiation under the bounded-support integral gives

    |m_j| <= |t_j| M_0 = D0_j,
    |m'_j| <= |t_j| w_j M_1 = D1_j,
    |m''_j| <= |t_j| (w_j² M_2 + w_j M_1) = D2_j.

These include physical cell extent. For orthonormal unit-cube cell
coefficients c_i at centers x_i and box size b, compute M0 exactly as
sum |c_i|/b^(3/2). Compute M2 exactly using |x_i|²+1/(4b²), and upper-bound
M1 by the square root of that quantity separately in every cell (Jensen).
Using only centers would omit the cell contribution.

For the mean contrast f(R,a)=a² sum p_j |m_j|² +
a³ Re sum b_abc m_a m_b conjugate(m_c), p=w_power/2 and
b=(w_real-i w_imag)/2. Product differentiation and the triangle inequality give

    H2 = 2 sum |p_j| (D1_j² + D0_j D2_j),
    H3 = sum |b_abc| [sum_cyclic D2_a D0_b D0_c
                      + 2 sum_cyclic D1_a D1_b D0_c],
    |d²f/ds²| <= H = A² H2 + A³ H3, for |a|<=A.

Let a finite rotation set cover SO(3) with geodesic radius r and let U be
the largest amplitude-profiled value on this set. Then the continuous maximum
is at most U+H r²/2. Proof: compactness supplies a maximizer (R*,a*).
At fixed a*, the rotational derivative vanishes at R* because SO(3) has no
boundary. Integrate the second derivative along a shortest geodesic to a
nearby grid point. Its value is at least f(R*,a*)-H r²/2, and profiling there
can only increase it. Nondifferentiability of the profiled envelope is harmless
because this argument fixes a*. If computed grid scores have certified error
eta, add eta to the upper bound. Ordinary floating-point comparisons alone
are not outward-rounded certificates.

A ZYZ Euler grid with n cyclic alpha/gamma nodes and m beta nodes including
0 and pi covers with r<=min(pi,2pi/n+pi/(2(m-1))). This follows by changing
the three factors separately and using the bi-invariant geodesic triangle
inequality. It deliberately overcounts near the Euler poles. For an allowed
remainder epsilon, choose n=ceil(3pi/sqrt(2epsilon/H)) and
m=ceil(3pi/(2sqrt(2epsilon/H)))+1 (with elementary lower size limits).
This constructs a sufficient covering. Its size is not an optimal covering
number or a lower bound on certification cost.

## Declared computation and stop rule

For all eight nonzero saved directions, compute M0/M1/M2, H2/H3 and H on
the identical region-removed cell maps and stated amplitude ranges. Record
all twelve cases including four zero directions. Preserve the existing
continuous-search margin, even if negative. For positive margins delta,
report a grid whose curvature remainder is at most delta/2. This would only
certify positive separation if the newly evaluated grid maximum plus its
score-error bound still lies less than delta/2 above the current sampled
maximum. That condition is not assumed and no grid scores are evaluated here.

If this construction needs more than 10 million rotations, do not run it.
Retain the count and require a quantitatively better envelope before more
enumeration. This rule limits one implementation; it establishes neither
intrinsic computational hardness nor statistical impossibility. Fixed-map
curvature says nothing about experimental noise, density error, amplitude
calibration or the continuous null maximum until the score grid is evaluated.
