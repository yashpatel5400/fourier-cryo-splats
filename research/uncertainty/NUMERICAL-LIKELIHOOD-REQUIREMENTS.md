# Numerical likelihood error and the eventual uncertainty target

1 October 2026 UTC. These elementary facts guide the current numerical gate.
They are not claimed as new statistical theorems or as a completed cryo-EM
uncertainty method. The current adaptive calculation has not yet finished.

## A data-dependent proposal does not bias the integral itself

Fix an observed image y and a map parameter theta. Let Haar measure H be
normalized to one, f_theta(R;y) the nonnegative conditional likelihood, and
q_y a probability density relative to H, positive wherever f_theta is positive.
The proposal may depend on y, the maps and earlier computations. Conditional
on this frozen proposal, let R_j be fresh independent draws from q_y H. Then

    I_hat_theta = M^(-1) sum_j f_theta(R_j;y) / q_y(R_j)
    E[I_hat_theta | y,q_y] = integral f_theta(R;y) dH(R).

The equality follows by cancelling q_y inside each expectation. A finite
variance requires integral f_theta^2/q_y dH < infinity. A defensive Haar
component supplies support; it does not make this variance small. Reusing
proposal-fitting samples as if they were fresh draws would require a separate
adaptive-sampling argument. The implemented gate uses new draws.

Neither log(I_hat_theta) nor I_hat_1/I_hat_0 inherits that unbiasedness. For
positive estimates, Jensen gives E log(I_hat_theta) <= log(I_theta).
Using common orientations for two models can reduce ratio variance but does
not establish a finite-sample error bound. These are general importance-
sampling facts, not a special advantage of Gaussian splats.

## Agreement and large ESS can miss the same region

For a general bounded integral, let q=1 and choose a set A of Haar mass p.
At one fixed observation define

    f_0(R) = p/(1+p),
    f_1(R) = (p + 1_A(R))/(1+p).

Both functions lie in (0,1], and their integral ratio is exactly 2. If neither
of two independent M-sample banks visits A, both estimated ratios equal 1 and
both empirical importance ESS values equal M. This event has probability
(1-p)^(2M), at least 1-2Mp. Thus any finite two-bank agreement criterion can
coexist with an arbitrarily high probability of a factor-two ratio error.
The proof is direct integration and the union bound. Smooth approximations
to the indicator preserve the issue, but this example does not establish
that the specific cryo-EM maps exhibit it, or that these functions themselves
are normalized image likelihoods. It is a limitation of the integration
diagnostic, not a cryo-EM identifiability theorem.

## What certified numerical intervals would buy

Suppose for every tested parameter theta a computation produces positive
bounds L_i(theta) <= p_theta(y_i) <= U_i(theta), jointly with probability
at least 1-delta over all numerical randomness, conditional on the images.
For any fixed normalized alternative density q_i(y_i),

    sum_i log q_i(y_i) - sum_i log U_i(theta)
       <= log product_i q_i(y_i)/p_theta(y_i)
       <= sum_i log q_i(y_i) - sum_i log L_i(theta).

For a composite null, replacing its likelihood supremum with a certified
upper bound also gives a conservative lower bound on the likelihood ratio.
Under the null and an independently fitted normalized alternative, the ideal
split likelihood ratio has expectation at most one. Markov's inequality and
a union bound then control rejection using the numerical lower bound by
alpha+delta. This is the established universal-inference argument combined
with a numerical failure budget, not a new coverage theorem. A local optimizer,
an ESS diagnostic or agreement of two banks does not provide the required
joint bounds. A pointwise-in-theta numerical statement also cannot silently
be promoted to a uniform confidence-set guarantee.

If absolute log-density errors are deterministically bounded by epsilon_i,
their contribution to a product likelihood is bounded by sum_i epsilon_i.
Consequently a .01 per-image agreement tolerance is only a feasibility gate;
it is not an adequate joint error budget for a stack of 100,000 particles.
Stochastic cancellation requires its own verified argument. Statistical
sampling uncertainty, numerical uncertainty and acquisition-model error
remain different quantities and all matter for the eventual study.

The [prior-art note](LIKELIHOOD-VALIDATION-READING.md) documents the primary
universal-inference reading and the normalized-density requirement. The
[current gate](ADAPTIVE-POSE-INTEGRATION-PROTOCOL.md) tests computational
stability only. A proposed uncertainty method must state a single scientific
estimand and address its model error before these facts can support a paper.
