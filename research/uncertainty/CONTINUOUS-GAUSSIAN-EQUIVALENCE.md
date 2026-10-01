# Continuous Gaussian comparison: derivation and scope

1 October 2026 UTC. Response to full review 2, findings B/R7/R10. This note develops the reviewer's observation; it is not a new statistical principle or an externally reviewed theorem. No empirical result is asserted here.

Let H=L2 of the unit cube, A:H→R^m fixed, ell∈H, and whitened Y=A(rho0+h)+epsilon with epsilon~N(0,I) and ||h||≤B. Write G=AA*, a=A ell and L=||ell||². For any fixed w, the centered affine estimator has variance s²=||w||² and worst absolute bias b=B||ell−A*w||. The generalized white prior with covariance tau²I gives weights w_tau=(G+I/tau²)^{-1}a. Its functional posterior variance is

v_tau=tau²(L−aᵀw_tau)=||w_tau||²+tau²||ell−A*w_tau||².

White covariance on infinite-dimensional H is not trace class; it defines an isonormal process on test functions, not a finite-energy random density in H. The finite-dimensional joint distributions used here are nevertheless valid. At a nonsolution w, the variance expression on the right equals the optimal posterior variance plus a nonnegative quadratic optimization error. Its centre is not the exact posterior mean. Numerical residuals must therefore be reported.

## A normal-envelope inequality

For z≥sqrt(3), t≥0, let q=z sqrt(1+t²). Then P(|N(t,1)|≤q)≥2 Phi(z)−1. Thus z sqrt(s²+b²) is a uniformly valid fixed-bias half-width whenever alpha≤2 Phi(−sqrt(3)) (about .0833). In particular this includes .05 and .05/12. This does not require w to be the optimal prior weight.

Proof. Set F(t)=Phi(q−t)+Phi(q+t)−1. For t>0 put a=z t/sqrt(1+t²). Then

F'(t)=phi(q−t)[(a−1)+(a+1) exp(−2qt)].

If a≥1 the derivative is nonnegative. If a<1 this is equivalent to atanh(a)≥qt=a/(1−a²/z²). Expand both sides as absolutely convergent power series. The coefficient inequality 1/(2k+1)≥1/z^(2k) holds for all k≥0, because z²≥3 and 3^k≥2k+1 (induction). Hence F is nondecreasing and F(0)=2Phi(z)−1. Pure noise and pure bias follow by continuity or direct evaluation. QED.

At tau=B the continuous prior interval therefore has the same uniform fixed-pose class coverage guarantee. At tau=B/2 it need not. A folded-normal critical value can shorten it at the same weights, and a different weight objective can change efficiency. Neither observation makes the underlying inferential construction novel. This implication corrects the previous blanket contrast between credible intervals and uniform coverage.

With a fixed Gaussian nuisance covariance Sigma=I+JJ*, replace ||w||² by wᵀSigma w in the prior identity and normal bound. This is valid for the stated linearized Gaussian nuisance model. It is not a guarantee for nonlinear unknown-density poses, biased alignment, or weights estimated from the same noise. Even the fixed-bias normal-envelope inequality need not apply if the bias adapts to the sign of the realized noise; the safe b+z s union/triangle bound handles bounded adaptive bias only when the Gaussian noise projection itself remains fixed.

Existing primary foundations are Donoho (1994), Armstrong–Kolesár (2018), and probabilistic Fourier-slice reconstruction (Ullrich et al., 2020). Backus–Gilbert resolving kernels and Gaussian-process equivalences will be incorporated with verified primary readings; reviewer's citations alone are not treated as source verification.
