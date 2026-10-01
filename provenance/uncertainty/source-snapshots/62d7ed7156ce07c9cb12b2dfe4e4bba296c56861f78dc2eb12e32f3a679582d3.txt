# What independent alignment images do and do not establish

1 October 2026 UTC. Elementary arguments for a diagnostic control; no claim of a new statistical principle.

## Conditioning obstruction

If a true pose theta is fixed and an estimated pose theta_hat is measurable with respect to a design sigma-field D, their error U=theta-theta_hat (in a fixed local chart) is D-measurable. Consequently E[U | D]=U. Conditional centering can hold only where U=0. Unconditional symmetry of an estimator's errors does not justify a conditional centered-error model after freezing its estimated design. A genuinely random latent-pose population can have a different conditional law, which must be derived rather than assumed. This observation does not rule out an unconditional analysis of dependent weights/errors.

## Phase-only counterexample

Let s>0 and, after rotating to the true frame, A=s+sigma_A(Z_1+i Z_2) and B=s+sigma_B(W_1+i W_2), with independent standard-normal coordinates. Set theta_hat=arg(A). The same-image aligned real observation is |A|, with E|A|>s by strict convexity away from a single ray. The independent-image aligned observation is Re(exp(-i theta_hat) B), whose expectation is s E cos(theta_hat)<s whenever sigma_A>0. Its variance is finite. Each bias persists as the number of independently aligned particles grows, even though the angular errors are symmetric.

By the law of large numbers, the corresponding sample means converge to the biased expectations. Any intervals around these means whose half-widths converge to zero have coverage of s tending to zero. In particular, estimating their actual sampling variance, or resampling only the aligned observations, cannot fix the bias. Independent images remove a noise-selection effect but leave attenuation. This is not an impossibility theorem for latent-pose reconstruction: invariant or marginalized methods use different estimators.

## An elementary invariant control with a finite-sample bound

For equal known coordinate noise sigma in A and B, let X=Re(conj(A) B). Unknown true phases cancel. Write m=s^2. Independence gives E X=m. For |t|<1/sigma^2, the product-normal Gaussian integral gives

    log E exp(t(X-m)) = -log(1-sigma^4 t^2) + m sigma^2 t^2/(1-sigma^2 t).

One way to verify this identity is to put U=(A+B)/(sqrt(2) sigma) and V=(A-B)/(sqrt(2) sigma). These independent real two-vectors have squared norms distributed as noncentral chi-square with two degrees and noncentrality 2m/sigma^2, and central chi-square with two degrees, respectively; X=(sigma^2/2)(||U||^2-||V||^2). This also gives Var(X)=2(sigma^4+m sigma^2).

For either sign of t, -log(1-u^2)<=u^2/(1-|u|), so the centered log-MGF is at most (sigma^4+m sigma^2)t^2/(1-sigma^2|t|). Independence across particles and the Chernoff argument for a sub-gamma variable imply, with x=log(2/alpha),

    P( |mean(X)-m| > 2 sqrt((sigma^4+m sigma^2)x/n) + sigma^2 x/n ) <= alpha.

Invert this inequality over m>=0, then take square roots of its endpoints. The resulting possibly empty interval covers s with probability at least 1-alpha for every fixed sequence of phases, under this known equal-variance Gaussian model. Its width need not be efficient. The confidence set is not converted to [0,0] when empty; that failure is saved. The script's numerical roots are ordinary floating point, not validated enclosures.

This is a standard concentration-bound specialization and an invariant-moment control. It does not solve unknown noise, CTF, heterogeneous density, motion/dose changes or the full SO(3) inverse problem. Classical alignment-envelope and cross-validation papers belong in the background, and the cryo-EM moment literature precludes presenting phase cancellation as a new idea.
