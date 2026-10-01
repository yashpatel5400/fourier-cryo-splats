# Direct folded-width optimization after review 3

This is a development comparison addressing R7. It changes the optimization
criterion of the fixed-pose audit; it does not resolve estimated-pose bias,
experimental calibration of the density class, or the three rejection verdicts.

## Classical reduction

Let G=AA*, a=A ell, s=||w|| and b=B||ell-A*w||. Let q(s,b) be the
folded-normal half-width at level 1-alpha, with alpha<1/2. Implicit
differentiation of the two-tail equation gives

    q_b = tanh(q b/s^2),       q_s = (q-b q_b)/s.

Both are nonnegative (q_s is positive when s>0), because q>b for alpha<1/2.
Any interior stationary optimum therefore solves

    (G+lambda I)w=a,          lambda=q_s b/(q_b B^2 s).

The endpoints include w=0 and the least-residual minimum-norm limit. More
generally the bias/variance Pareto frontier lies on this ridge path: minimizing
||w||^2 subject to a residual-norm bound gives the same quadratic Lagrangian,
except at its endpoints. Searching the path is a classical fixed-length
confidence-interval calculation, not a new uncertainty principle. The old
sum-width criterion uses q_s=z_(1-alpha/2), q_b=1 instead.

[Armstrong and Kolesar (2018), section 3.4](https://www.princeton.edu/~mkolesar/papers/optimal.pdf)
explicitly minimize the two-sided fixed-length half-width over the bias/variance
tradeoff. Their general result credits Donoho (1994). The related
[Armstrong, Kolesar and Kwon draft, section 2.2](https://tbarmstr.github.io/files/2023/08/regularized_regression_draft.pdf)
derives generalized ridge estimators under an L2 restriction. These targeted
sections were checked; this is not a new full-paper replication or a novelty
claim based on renaming ridge regression.

## Search and numerical lower diagnostics

For the exact ridge solution, s decreases and b increases with lambda.
For lambda in [L,R], q(s(R),b(L)) is a lower bound on the attainable width
on that interval. Retain this bound when adaptively bisecting the interval
with the smallest lower value. At an approximate solve with residual r,

    ||w_hat-w_lambda|| <= ||r||/lambda,
    ||A*(w_hat-w_lambda)|| <= ||r||/(2 sqrt(lambda)).

The latter inequality maximizes sqrt(g)/(g+lambda) over g>=0. Add the
continuous-Gram quadrature action bound to the computed CG residual, and
use the existing field-norm quadrature bound for residual-norm endpoints.
Rounding and NUFFT error have diagnostic pads only, not validated enclosures.
The numerical lower values must carry this qualification.

The initial tests demonstrated that simple outer-ray lower bounds can remain
loose when the target has an unobserved component. The first test output is
retained. Consequently the declared range is [lambda_initial/64,
64 lambda_initial]. Stop at relative gap .005 within this range, 80 solves,
or 180 seconds per target. The no-data endpoint is always a candidate. Keep
the outer-ray bounds and global relative gap as separate diagnostics; do not
call restricted convergence a global certificate. No inference pixels are
used for selecting lambda.

## Six molecular comparisons

Reuse the true-pose design of the saved local-alignment generators: 128
particles, frequency radius giving 40 nonredundant complex coefficients,
each of EMPIAR-10028/10049/10076. Use the center and contrast targets at width
.07 of the field, B=2 and alpha=.05. Use the same continuous quadrature order
40, preconditioner rank 1024, CG tolerance 1e-10 and 2,000-iteration limit.
The initial ridge is the last sum-objective ridge from each saved true-pose
fit. Report every case, all solves/statuses, true-pose sum-objective width,
both matched Gaussian prior-scale widths, the actual selected folded width,
bias/SD, restricted/global gaps and timing. Save the selected weights.

This is not another noise-trial coverage study. The known-pose class coverage
comes from the existing conditional normal calculation; whether an actual
experimental density lies in B=2 remains unestablished. Frozen trial results
and earlier releases are not overwritten.
