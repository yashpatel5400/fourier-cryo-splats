# Adding noise is not a substitute for validating independence

30 September 2026, after full-review-2 submission. This is prior-art reading
and a diagnostic identity, not an implemented method or a new novelty claim.

[Neufeld et al. (JMLR 2024)](https://www.jmlr.org/papers/v25/23-0446.html)
develop data thinning for convolution-closed distributions. Their Gaussian
construction requires its covariance parameter. Section 2.3, Proposition 10
explicitly shows that using the wrong Gaussian variance leaves the two parts
correlated; adding too much variance changes the sign of that correlation,
rather than making it vanish. The source, selected reading scope and local
PDF digest are recorded separately. The entire paper and author code were
not reviewed or executed.

In a convenient rescaled Gaussian notation, let Y = mu + e, with
e ~ N(0,Sigma), and independently Z ~ N(0,Sigma_z). For fixed tau > 0 define

    Y_A = Y + tau Z,       Y_B = Y - Z/tau.

Direct covariance expansion gives Cov(Y_A,Y_B) = Sigma - Sigma_z. Joint
Gaussianity implies independence when these covariances match. A conservative
upper bound on Sigma does not itself supply equality. This is a restatement
of the classical construction and nuisance-misspecification issue, not a new
theorem. It cannot justify treating previously picked, normalized and
pose-refined particles as independent Gaussian pseudo-halves without modeling
those operations. In particular, our directional scalar variance upper bounds
are not a known full covariance matrix for such a construction.

The raw-frame possibility and this synthetic alternative therefore address
different prerequisites. Neither removes the need to define the scientific
estimand, account for nuisance/data dependence, and retain realistic negative
controls. No fission-generated data or favorable outcome is introduced here.
