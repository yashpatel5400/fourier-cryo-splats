# Prior art for the viewing-variance refinement

1 October 2026 UTC. These methodological source cards supplement the 138-entry
cryo-EM/inverse-inference candidate ledger; they are not counted as additional
fully reviewed cryo-EM papers. The frozen simulation protocol remains unchanged.

**Maurer and Pontil (COLT 2009).** The
[primary conference PDF](https://www.learningtheory.org/colt2009/papers/012.pdf)
was read selectively at Theorem 4, Section 2's variance calculation and the
Theorem 11 proof. The sample-variance normalization agrees with the unbiased
variance used in our code. Their inequality is applied to independent view
groups after rescaling the signed cross-product. We have not reproduced the
remaining uniform-learning results, experiments or full proof chain. The
concentration inequality is credited rather than presented as our theorem.

**Duchi and Namkoong (JMLR 2019).** The
[journal listing](https://jmlr.org/beta/papers/v20/17-750.html) gives volume
20(68), pages 1–55, despite the downloaded PDF's old volume-19/2018 header.
The [primary paper](https://www.jmlr.org/papers/volume20/17-750/17-750.pdf)
was read selectively in the introduction, Section 2.1 and Theorem 1. Their
chi-squared distributional-robustness construction gives a mean-plus-standard-
deviation envelope, with conditions for equality. Our population density-ratio
bound uses the same elementary Cauchy--Schwarz mechanism, not a newly discovered
robustness principle. The paper's empirical-distribution neighborhood and
convex-learning analysis are different from our specified Haar simulator and
noise-replicated conditional event means. Its remaining proofs, implementation
and classification experiments have not been replicated.

**Janon et al. (ESAIM: Probability and Statistics, 2014).**
[Published metadata](https://numdam.org/articles/10.1051/ps/2013040/) confirms
volume 18, pages 342–364 and DOI 10.1051/ps/2013040. The
[author manuscript v1](https://arxiv.org/pdf/1303.6451v1), introduction and
Sections 1.1–1.2 including Lemma 1.2, was read selectively. Holding one input
fixed while independently resampling the others yields a covariance equal
to the variance of the conditional mean. This directly precedes our paired
conditional-noise product. We use the unnormalized conditional variance with
a finite concentration bound, not their normalized Sobol-index estimator or
its asymptotic confidence intervals. The metamodel-error results and numerical
examples were not audited. Conditional-noise replication itself is prior art.

**Additional classical test context.** Dufour (2006) and Berger–Boos (1994)
are described, with their distinct reading scopes, in the preceding
[Monte Carlo protocol](MONTE-CARLO-VIEW-LAW-PROTOCOL.md). They delimit novelty
in simulation testing and nuisance maximization as well. A favorable new
simulation result must therefore establish a useful problem-specific method
and comparison, not rename these established ingredients. None of these
sources establishes our experimental noise or viewing-distribution assumptions.
