# A residual bound for an unfinished continuous Gaussian solve

1 October 2026 UTC. Numerical diagnosis only; retain the original declared CG convergence flags. This standard quadratic-error identity is not a new inference principle.

Let G=A A*, a=A ell, Sigma>=I, H=G+Sigma/tau^2 and v(w)=w*Sigma w+tau^2||ell-A*w||^2. If w*=H^{-1}a, then

    v(w)-v(w*) = tau^2 (Hw-a)* H^{-1} (Hw-a).

In exact arithmetic, a truncated pivoted-Cholesky factor F of the positive-semidefinite kernel G satisfies FF*<=G. Consequently P=FF*+I/tau^2<=H, and H^{-1}<=P^{-1}. This gives the upper bound tau^2 r*P^{-1}r, where r=Hw-a. It can be much smaller than tau^4||r||^2. The diagonal-corrected preconditioner does NOT necessarily satisfy P<=H and cannot be used for this inequality; use only the original truncated-Cholesky inverse.

If a quadrature residual r_q has a known Euclidean error at most e, the square root of r*P^{-1}r is at most sqrt(r_q*P^{-1}r_q)+tau e. Apply the existing continuous Gram-action integration remainder. Floating-point Cholesky/NUFFT errors are not enclosed by this identity; report the result as a numerical optimization-gap diagnostic. A small variance gap does not on its own bound the conditional posterior-mean error for an arbitrary observed image vector. Nonconverged intervals must not be relabelled exact posterior intervals.

The first probe uses only the first completed fixed-pose 10028 Gaussian fit, its existing design/target and weights. It evaluates no new image pixels or coverage outcomes, and does not alter the 48-fit protocol.
