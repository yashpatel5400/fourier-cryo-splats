# Archived design hypothesis, now implemented in development

The plan below was written before implementation. See THEORY.md Section 15,
`src/fourier_splats/uq_continuous_pose.py`, and the independent numerical tests.
The current implementation uses norm-bounded pilot contributions and explicit
polynomial-Fourier quadrature errors; full development experiments are running.
It is not a frozen confirmatory result.

# Proposed continuous pose extension — not implemented or established

This is a concrete next development hypothesis. Do not describe it as an
available guarantee until derivation, implementation and independent numerical
checks are complete. Current continuous code is fixed-pose only; current
nonlinear pose code is finite-voxel only.

## Conservative starting point

For unit pose balls, write the scalar-estimator adjoint change as

 g(u,x)=sum_i,a D_ia^* w_i(x) u_ia
       + 1/2 sum_i,a,b E_iab^* w_i(x) u_ia u_ib + remainder.

With ||delta||<=B, ||rho0||=P, residual h=ell-A0*w and rho=rho0+delta,

 |bias| <= B||h|| + (B+P)(L1 + L2/2) + remainder,

where L1 bounds ||sum_i T_i u_i||, and L2 bounds
||sum_i U_i(u_i tensor u_i)|| in continuous L2 of the cube. This deliberately
uses the pilot norm rather than its exact derivative contractions; it can be
conservative, but avoids claiming that point-sampled pilot derivatives equal
continuous ones. Weighted block spectral inequalities can bound L1 and L2,
retaining the density/pose field shared across particles. Initial weights can
be the continuous fixed-pose optimized weights, which already respect the full
cube support and avoid the bad sphere-to-cube transfer seen in the probe.

For the unit cube, R=sqrt(3)/2. With rotation radius a and translation radius
s in physical field units, each Fourier column has uniform phase bounds

 L=2*pi*(a*||k||*R+s*||q||), H=2*pi*a^2*||k||*R,
 T=2*pi*a^3*||k||*R.

The continuous Hilbert-Schmidt third-derivative bound per particle is
K3_i=sqrt(sum_q |C_iq/noise|^2*(L^3+3LH+T)^2), since the cube has volume one.
A valid scalar remainder is (B+P)/6 sum_i K3_i ||w_i||. This must be checked
carefully against independent nonlinear projections. A shift of 0.01 working
pixels on a 24 grid is s=0.01/24 in the physical field; do not change units
when evaluating quadrature nodes.

## Polynomial Fourier quadrature error

D^*w is a Fourier field with polynomial degree at most one; E^*w has degree at
most two. Their Gram integrands therefore have degree at most four times
exp(i nu x). For one coordinate and beta<=4, Leibniz gives derivative bound

 sum_{l=0}^{min(beta,2r)} binom(2r,l) beta!/(beta-l)!
      *(1/2)^(beta-l) * Lfreq^(2r-l),

with Lfreq=4*pi*max|k_axis|. Multiply by the standard C_r and sum coordinate
bounds (other monomial amplitudes <=1 on the cube). Taking the maximum over
beta=0,...,4 yields a uniform kernel error E_poly for the required monomials.
If column a has absolute Fourier-polynomial coefficient sum M_a, then its Gram
entry error is bounded by E_poly*M_a*M_b, hence the Gram operator error by
E_poly*||M||^2. Add this to the largest quadrature Gram eigenvalue before taking
a spectral norm. Block scaling must also scale M. This is a proposed bound;
verify constants, realification and tensor-product telescoping explicitly.

For first derivative columns, the polynomial-coefficient absolute sums are
bounded by sum_j |c_j| S_ja, where rotation S entries are
2*pi*a*(|k_b|+|k_c|), and translation entries are 2*pi*s*|q_axis|.
Second derivative sums can use sum_j |c_j| (S_ja S_jb + curvature_ja,b), with
curvature from 2*pi*a^2*((k_a x_b+k_b x_a)/2 - delta_ab k'x).
Keeping duplicate ordered Hessian pairs matches ||u tensor u||=||u||^2.

## Numerical plan

1. Validate polynomial Fourier quadrature remainder on small direct integral
   examples and compare with fine Gauss integration/exact derivative integrals.
2. Reuse the twenty Fourier-moment formulas at nonuniform quadrature nodes in
   physical units, using direct sums or type-3 transforms (the existing voxel
   helper's NUFFT branch assumes an integer grid and cannot be silently reused).
3. Store first/second columns multiplied by sqrt(quadrature weights), with
   positive error padding of their spectral norms. At order40, 128 particles,
   5+25 columns per particle, this is about 2 GB plus temporary matrices.
4. Evaluate continuous fixed-pose weights at 0.1 and 0.5 degrees before a larger
   sweep. Include the no-data option and report if bounds become uninformative.
5. Validate actual nonlinear changes using independent constant-cell Fourier
   generators, with rotated k, exact cell sinc factors, half-cell phase, fixed
   CTF and explicit translations. No finite-voxel coverage labels should be
   substituted for this continuum model.
6. Only after useful verified results, incorporate this into a frozen protocol.
