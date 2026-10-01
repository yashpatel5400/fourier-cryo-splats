# Full-frequency paired covariance: post-audit feasibility gate

This follows the 18-cell compressed screen and the authentic focused Fable
audit `paired-statistics-audit-01`. It does not hide those results or treat
small lower values as impossibility. In particular, the rank-32 bases retain
only .642–.791 of total true signal energy and may omit the removal's signal.

## Two-sided expected-growth diagnostic

Use all 440 real coordinates from the 220 nonredundant complex frequencies.
There is no projection or selected basis. For null means m_l and alternative
second moment S, choose any coefficients lambda_l>=0 and put
A=sum lambda_l m_l m_l'. Set Delta=S-A. For every null-feasible T,

    growth(T) <= tr(T Delta/v) + (1/2)log det(I-T^2)
              <= sum_j g(eigenvalue_j(Delta)/v),

where g(d)=r-(1/2)log(1+r), r=(sqrt(1+4d^2)-1)/2. The first step uses
tr(TA)<=tr(UA)<=0. For the second, diagonalize Delta and use spectral
alignment followed by the scalar maximization over -1<t<1. The maximizer is
t=2d/(1+sqrt(1+4d^2)). This is a standard conjugate/spectral argument,
suggested as the missing diagnostic by the independent audit, not new general
theory. It holds for every feasible mixture, including an unconverged iterate.

Approximate the nearest cone point in Frobenius norm with nonnegative least
squares. Trace-normalized outer products have Gram kernel (u_l'u_k)^2, so
4,160 views require a 4,160-square kernel rather than an explicit
97,020-by-4,160 feature matrix. The symmetric dimension at 440 is
440*441/2=97,020. Use L-BFGS-B with 3,000 iterations or 120 seconds per
candidate, retaining the best nonnegative iterate and every stopping reason.
The residual supplies a separator only at the exact projection. At any
iterate repair the largest normalized finite-view violation by subtracting
that amount times identity plus a 1e-12 numerical margin. Search this repaired
direction as in the initial matrix protocol. The feasible lower and spectral
upper quantify the remaining optimization gap; neither solver success nor
an expected value is measured rejection power.

## Six cases and fresh views

The stacks are 10028,10049,10076. Use true-map control and full local removal,
all 4,160 saved views for both null and alternative, and the same first
CTF/noise profile. Assert no zero, duplicate or opposite frequency in the
nonredundant set; no Nyquist coordinate is present at radius 12 in a 64-pixel
box. Identity noise covariance remains a stipulated model, not calibrated
experimental whitening. The true-map mixture can be checked directly using
equal view probabilities; optimization is only needed for full removal.

For each removal, generate 10,000 new Haar orientations with seed
261005+dataset ID. Use the same continuous-cell Fourier evaluator and density
construction. Check the first 64 saved templates before any fresh views.
Evaluate the selected direction on all fresh views, recording the fraction
with normalized violation above 1e-10, maximum violation and a bound on the
null exponent. Repair against the union and recompute the directional growth.
This is an off-catalog falsification check, never a continuous certificate.
With unlimited amplitude any genuinely positive violation invalidates the
proposed continuous null, regardless of how small it is.

At v=1,2,4 report growth per particle, a 128-particle comparison to the old
table, and normal-approximation rejection probabilities at n=1,000,10,000,
100,000. Compute the exact variance of a single log statistic for independent
unit Gaussian exposure noise and a uniform catalog viewing law:

    Var(X'KZ)=tr(K^2)+2 tr(K^2 S)+Var_R(m_R'K m_R), K=T/v.

The normal approximation is a development diagnostic. It is not a simulation,
a finite-sample power guarantee, or evidence of calibrated real data.
Include the known-pose, known-unit-amplitude two-exposure likelihood-ratio
oracle's mean KL, mean_R ||m_true-m_removed||^2, as an information comparator.
It has stronger nuisance knowledge than the proposed method.

Preserve sources, inputs, means, mixture coefficients, directions, selected
weights, fresh Fourier values, all results and failures. Do not launch large
coverage trials before these information and fresh-view checks. Translation
distributions, dose/gain mismatch and paired empty-ice dependence remain
separate unresolved gates. No ICML acceptance claim follows from this screen.

## Numerical correction from the audit

The old scalar diagonal dual evaluator reports f(x_hat) at an approximate
stationary root; its root error can make that slightly smaller than the true
supremum. A new evaluator adds the maximum of the concave tangent at both
open-domain endpoints and a diagnostic rounding pad. Frozen outputs are
retained. Replaying those bounds with the corrected evaluator is a separate
correction, not an alteration to old trial data. The covariance spectral
formula also reports its raw value and rounding pad separately.
