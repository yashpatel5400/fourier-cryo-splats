# Continuous Gaussian-orbit likelihood bounds: development derivation

30 September 2026. This is an implementation proposal, not an empirical result.
The normalized split-likelihood testing principle is established universal
inference. Mixture likelihood duality, finite-support maximizers and the scaled
dual gap below are classical (Lindsay, 1983, sections 2--5). We do not claim
these as new statistical theorems. The proposed work is a tractable continuous
orientation envelope for a specified Fourier Gaussian cryo-EM model; its
computational usefulness and acquisition-model adequacy remain to be tested.

## Model and continuous dual

Fix a candidate Hermitian Fourier Gaussian map, known CTFs and independent
isotropic Gaussian image noise. For now translations are exactly zero and the
noise scale is supplied. Let p_i(R) be the likelihood of image i at rotation R.
The unknown shared viewing measure G gives log likelihood
`L(G)=sum_i log integral p_i(R) dG(R)`.
For any strictly positive anchors z_i, and any M satisfying
`M >= sup_R sum_i p_i(R)/z_i`, Jensen gives

`sup_G L(G) <= sum_i log z_i + n log(M/n)`.

Indeed, apply arithmetic--geometric mean to the n positive integrals divided
by their anchors, and then interchange a finite sum and integration. This is
the optimally scaled tangent bound. For feasible anchors and M=n+delta it is
Lindsay's `n log(1+delta/n)` gap, tighter than the unscaled delta bound.
Neither the anchors nor a finite orientation grid alone certify M. A finite
mixture supplies a feasible lower likelihood, not a continuous upper bound.
The compact-orbit likelihood curve has a maximizing mixture supported at at
most n rotations (Lindsay's result includes independent nonidentical kernels).
This does not identify those support locations or make their optimization easy.

## A box covering every rotation

Use `R(a,b,c)=Rz(a) Ry(b) Rz(c)` on
`[0,2pi] x [0,pi] x [0,2pi]`. Every SO(3) rotation has a representative.
The nonunique boundaries do not invalidate an upper bound; a measurable
representative can assign every rotation to one box. Bisecting closed Euler
boxes preserves coverage. A box center and half-widths h have path length at
most H=sum(h), so frequencies of norm K move by at most
`2 K sin(min(H,pi)/2)` from their center values. Along a straight segment in
Euler coordinates, frequency speed and acceleration are bounded by KH and
KH², respectively. The latter follows by differentiating the three rotation
factors: all diagonal and mixed derivative norms sum to H².

For a real isotropic Fourier Gaussian `g(k)=exp(-||k-mu||²/(2s²))`, over the
enclosing frequency ball let distances to mu range over [l,u]. Then
`B1=max_[l,u] r/s² exp(-r²/(2s²))` bounds the gradient norm; the maximum is
at r=clip(s,l,u). The conservative Hessian bound
`B2=max_[l,u] (1+r²/s²)/s² exp(-r²/(2s²))` has the same maximizing r.
Sum these bounds with absolute complex coefficients over both conjugate
Gaussian centers. Known CTF/whitening factors multiply absolute bounds.

The image mean in the Euler box therefore satisfies
`m(theta+delta)=m(theta)+J delta+e`, with
`||e_i|| <= b_i = .5 H² || |CTF_i| (K² B2+K B1) ||_2`.
Complex norms equal the Euclidean norms of real/imaginary concatenation.
The Jacobian J is evaluated exactly from the Gaussian and Euler derivatives.
No truncation of Gaussian tails is permitted in this calculation.

## Linear residual lower bound with an explicit dual

For image residual r and Euler box |delta_j|<=h_j, the convex linearized
least-squares problem has dual lower value
`d(u)=u^T r - .5||u||² - sum_j h_j |(J^T u)_j|`.
This follows by minimizing the quadratic Fenchel representation over the box.
Every u is dual feasible. In particular, for u=r-Jv,
`d(u)=.5||r||²-.5 v^T(J^T J)v - sum_j h_j |J^T r-(J^T J)v|_j`.
Coordinate descent may choose v, but its convergence is not needed for a
lower bound; retain the best dual value including u=r and u=0. Thus

`min_box ||r-J delta-e|| >= max(sqrt(2 max_u d(u))-b_i,0)`.

The resulting negative half squared residual upper-bounds the Gaussian log
kernel, excluding its rotation-independent normalization. Also retain the
first-derivative mean-ball bound and take the larger residual lower bound.
All statements are in real arithmetic; floating-point validation is separate.

An additional coordinate-wise enclosure keeps uninformative frequency
coordinates from being absorbed freely into a large image-wide error ball.
At fixed K=||k|| and M=||mu||, a spherical cap of angular radius rho gives
`dot(k,mu)` in `[KM cos(min(alpha+rho,pi)), KM cos(max(alpha-rho,0))]`, where
alpha is the angle at the center. The conjugate pair's real basis is
`2 exp(-(K²+M²)/(2s²)) cosh(dot/s²)`, even and increasing in absolute dot;
its imaginary basis is the analogous sinh, increasing in dot. Their interval
extrema are therefore explicit. Sum intervals with signed real/imaginary
coefficients and CTFs. The Euclidean distance from each observation to the
resulting rectangle is another residual lower bound. Take the maximum with
the Taylor bound. This does not assume that frequencies or Gaussian terms can
actually attain their interval extrema simultaneously. That independence is
a conservative relaxation. Stable evaluation uses two bounded exponentials
instead of a potentially overflowing cosh/sinh factorization.

## Adaptive computation and its stopping interpretation

Each leaf Euler box carries upper log kernels for all images and an exact
center kernel. A finite mixture over evaluated centers supplies anchors.
For these anchors compute the continuous M upper bound by maximizing
`logsumexp(log_upper_cell - log_z)` over leaves. Split the largest such leaf
along its widest Euler coordinate, retain both children, and periodically
refit the finite mixture on all evaluated centers. Keep the best valid global
upper bound obtained, because any past anchor remains admissible.
The best feasible likelihood and best continuous upper bound bracket the
global likelihood at every completed cover. A wall/cell limit returns that
bracket with an unresolved status, never the finite-grid lower value as an
upper likelihood. Preserve the final covering partition and every anchor
needed to replay the best bound.

Before any biological-geometry calculation, compare analytic derivatives with
independent finite differences, box envelopes with sampled rotations, the
linear-residual dual with independent constrained least squares, and the
adaptive result with a dense orientation grid on a small phantom. These tests
can detect implementation errors but do not prove floating-point enclosure.

Unknown translations, common noise-scale optimization, colored covariance,
uncertain CTFs and heterogeneous ensembles are additional tasks. Existing
discrete-scale success does not automatically transfer to this continuous
algorithm. A useful test also needs an independently learned normalized
predictor; a good denominator alone does not supply that predictor.
