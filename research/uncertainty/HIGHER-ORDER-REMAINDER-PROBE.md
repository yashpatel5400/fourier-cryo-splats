# Prospective higher-order remainder feasibility check

30 September 2026, after the first pilot-selected one-degree audit. That case's
expanded-cube cubic remainder is 3.5496 of a total bias bound 5.1554; its sign
power remains negligible. This motivates a bounded diagnostic before investing
in a higher-order pose-field implementation.

Use exactly two existing estimators: the 10049 first pilot-selected region at
20 A Gaussian sigma (128 particles), and the previously completed 10049 central
10 A sigma estimator (1,024 particles). Both use radius 12. Preserve their
weights, one-degree/0.5 A joint pose class and B=2, P=1 density bounds. Compute
remainders for Taylor degrees 2, 3, 4 and 5 on the same expanded cube. Retain
every degree and both cases. No fitting, reference-driven degree choice or new
calibration data is part of this diagnostic.

For the known real Fourier adjoint field f, let S_j be the Frobenius norm of its
j-th derivative tensor in L2 over the enclosing domain. For a pose path
y(t)=exp(t Omega)x+t b, set M_1=hypot(a sqrt(3)/2, ||b||_bound) under the joint
ball, and M_j=a^j sqrt(3)/2 for j>=2. The rigid change of variables preserves
volume. Repeated chain differentiation and Cauchy--Schwarz therefore give

    ||(f o y)^(k)||_L2(cube) <= sum_{j=1}^k Bell_{k,j}(M_1,...,M_k) S_j,

where the partial exponential Bell polynomials sum products of M values over
set partitions. Dividing by k! bounds the remainder after degree k-1. This is
the standard higher-order chain rule, not a new general statistical theorem.
The exact Fourier seminorm identity is already implemented for arbitrary j.

Approximate detector embeddings require a residual phase pad. If its absolute
size is at most e and Q_r bounds the r-th derivative of the main exponential,
Leibniz differentiation of (exp(i t delta)-1) times that exponential gives

    e Q_k + sum_{j=1}^k binom(k,j) e^j Q_{k-j},    Q_0=1.

The Q_r are complete Bell polynomials of the nonnegative phase-derivative
bounds 2 pi ||k_frequency|| M_j. This pad is summed with coefficient magnitudes
and divided by k!. Ordinary floating-point guards remain heuristic, as in the
existing enclosing-domain probe.

Check degree two against the completed implementation, verify Bell coefficients
through order four, and compare independent direct nonlinear fields with Taylor
polynomials through degree five on small unrelated examples. Then run the two
declared estimator diagnostics. Higher derivative fields are NOT implemented
in the confidence interval yet: replacing the old cubic term with a smaller
number while omitting those fields is invalid. Do not report a new interval,
coverage probability or sign power from this diagnostic alone.
