# Second-order log-likelihood enclosure: prospective refinement

30 September 2026, after the first full-orientation case's large gap and its
anchor diagnostic. Existing runs remain unchanged. The aim is to retain
data-dependent curvature rather than charge all mean curvature as error.
This uses standard Taylor and concave-majorization arguments, not a new
statistical inference principle.

For a single image with whitened mean m(a), residual r=y-m(a), and Euler
increment d in a box |d_j|<=h_j, let H0=sum(h). The log kernel is
`ell(a)=-.5||y-m(a)||²`. Its gradient is `g_j=<r,m_j>`, and its Hessian is
`H_jk=<r,m_jk>-<m_j,m_k>`, using real inner products of real/imaginary
concatenations. Exact Gaussian Fourier gradients/Hessians and first/second
Euler matrix derivatives give these quantities without finite differences.

Using the enclosing frequency ball from the first derivation, bound the
Gaussian third differential by
`(r³/s⁶+3r/s⁴) exp(-r²/(2s²))`. Its maximum over a radial interval is attained
at `r=clip(3^(1/4) s,lo,hi)`. The derivative of `(x³+3x) exp(-x²/2)` has
sign `3-x^4`, which establishes this maximum. Sum absolute coefficient bounds
over both conjugate centers to obtain B1, B2 and B3.

For each Fourier coordinate after absolute transfer/whitening, define
`M1=|C| K B1`, `M2=|C| (K² B2+K B1)`,
`M3=|C| (K³ B3+3K² B2+K B1)`.
Differentiating a product of Euler rotations a third time bounds frequency
jerk by K H0³. Along the segment the complex residual magnitude is at most
`|r(0)|+H0 M1`. Since
`ell'''=<r,m'''>-3<m',m''>`, a valid scalar log-kernel remainder is

`T3=H0³/6 sum_q [(|r_q(0)|+H0 M1_q) M3_q+3 M1_q M2_q]`.

Thus `ell(a+d) <= ell(a)+g^T d+.5 d^T H d+T3`.
This remains conservative but is third order in box diameter. It may be
larger than the original bounds, so retain their pointwise minimum too.

To upper-bound the quadratic on the box, choose
`lambda>=max(lambda_max(H),0)` and `Q=lambda I-H`, which is PSD. For any v,
the tangent to the concave function `g^T d-.5 d^T Q d` gives

`max_box [g^T d+.5 d^T H d] <=
 .5 v^T Q v + sum_j h_j |g_j-(Qv)_j| + .5 lambda ||h||²`.

Coordinate ascent in the concave surrogate selects v. Retain v=0 and the
best tangent bound, so its optimization need not converge. Eigenvalue
computation is ordinary floating point with a small positive numerical pad;
it is not a validated eigenvalue enclosure. The proofs concern real
arithmetic. Independent tests must check derivatives, the quadratic upper
bound and nonlinear cell envelopes before any biological-geometry probe.

This is a local bound-quality proposal. It supplies neither a global pose
cover by itself nor a learned predictor, uncertain-CTF model, noise calibration
or density-coverage guarantee. A local improvement would still require a
separate complete global computation.
