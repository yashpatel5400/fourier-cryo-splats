# Algebra audit for the proposed population screen

1 October 2026. Elementary derivations, not a novelty claim. All image vectors
below are realified complex Fourier coordinates. The inner product is
Re sum(conj(u) v); each real coordinate has noise variance one.

## Exact known-pose misspecification identity

Let state A have probability p, recorded/known orientation R have conditional
laws nu_A and nu_B, and conditional image densities f_A(y|R), f_B(y|R) be
correct. A pooled conditional likelihood uses f_t=t f_A+(1-t) f_B while
ignoring dependence between state and R. Define

    J_t(R) = integral (f_A-f_B)^2 / f_t dy.

For 0<t<1, integration of f_A-f_B=0 gives

    E_A[(f_A-f_B)/f_t | R] = (1-t) J_t(R),
    E_B[(f_A-f_B)/f_t | R] = -t J_t(R).

Consequently the population score is

    S(t) = p(1-t) E_nu_A[J_t] - (1-p)t E_nu_B[J_t].

At an interior pseudo-true optimum t*, S(t*)=0 and

    t*/(1-t*) = p/(1-p) E_nu_A[J_t*] / E_nu_B[J_t*].

This is an **implicit equation**, not a closed-form constant odds correction.
Under strict concavity and identifiability, t*=p iff the two expectations of
J_p agree. Different viewing laws alone do not imply bias. Conditional score
curvature is minus the expected squared score, so the population objective is
concave; it is strictly concave when the components differ with positive
probability. The identity is for observed poses and correct conditional
densities. It is not a derivation for a latent-pose cryo-EM likelihood.

For Gaussian means m_A,m_B and covariance I, let g=m_A-m_B and s=||g||^2.
Under B/A the log density ratio is sqrt(s)Z minus/plus s/2, Z~N(0,1).
Writing h_t(l)=logistic(logit(t)+l),

    J_t(s) = E[h_t(sqrt(s)Z+s/2)-h_t(sqrt(s)Z-s/2)]/[t(1-t)].

One-dimensional Gauss-Hermite quadrature therefore evaluates the exact
known-pose score without rendering new images. The weak-separation expansion
J_t(s)=s+O(s^2) yields the approximate fitted fraction
p E_A[s]/(p E_A[s]+(1-p)E_B[s]). Its validity is checked against the exact root;
it is not silently assumed for the measured s values.

## What amplitude projection does and does not establish

At the nominal full map, a real multiplicative amplitude has tangent m_A.
Projecting the state mean difference off this tangent gives

    s_perp = ||g||^2 - (Re <g,m_A>)^2 / ||m_A||^2.

This is the Schur complement for a *Gaussian mean* perturbation with one
real amplitude nuisance at that point. Complex amplitudes would give a
different projection and are not physical here. Neither s nor s_perp is
automatically the efficient information for a binary population mixture.
The secondary projection at m_B is also recorded; neither repairs arbitrary
view- and state-dependent amplitudes. A local nonzero tangent can coexist
with nonlinear ambiguity or severe instability.

If amplitude is fixed at one and poses are known, the score

    psi(Y,R) = Re <Y-m_B(R),g(R)> / ||g(R)||^2

has conditional mean 1{state=A}. For any joint state–view law and zero-mean
unit-covariance noise independent conditional on that law,

    Var(mean psi) = [p(1-p) + E_nu(1/s)]/n,
    nu = p nu_A + (1-p) nu_B.

The law of total variance proves this formula: the conditional mean is the
Bernoulli label and the conditional variance is 1/s. It requires s>0 and a
finite harmonic moment. Replacing s by s_perp does **not** prove the same
formula under unknown amplitude. The screen uses that replacement only as
an optimistic local diagnostic, with exact nominal threshold
10000*.03^2 - .75*.25 = 8.8125.

Indeed a fixed affine function of Y cannot have state-specific expectation
uniformly over a nontrivial interval of multiplicative amplitudes: its
expectation is a <v,m_z>+b, so being constant in a forces <v,m_z>=0 for
both states and leaves the same b. A nonlinear feature response method must
explicitly allow/calibrate approximation bias. Numerical identifiability is
not a substitute for that step.
