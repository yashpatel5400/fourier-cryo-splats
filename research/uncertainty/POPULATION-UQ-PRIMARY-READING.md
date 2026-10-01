# Population uncertainty: direct prior art and remaining questions

1 October 2026 UTC. This reading was performed while full review 4 was running
and is outside its immutable packet. No population method or dataset result is
claimed. The ledger gains one previously missing primary paper (139 candidates;
not 139 full readings).

## Cryo-BIFE (2021)

[Giraldo-Barreto et al., Scientific Reports](https://doi.org/10.1038/s41598-021-92621-1)
already infer population/free-energy uncertainty from individual particle
likelihoods. The primary JATS text was retrieved from Europe PMC. We read the
introduction, theory, experimental TMEM16F example, discussion, BioEM and MCMC
methods, experimental-data methods, and data availability. Synthetic-results
sections, supplementary material, images and implementation were not audited.

The method fixes a conformational path, marginalizes nuisance parameters using
BioEM and samples the resulting free-energy posterior. Orientations have a
uniform prior; numerical evaluation uses a coarse search and local refinement.
The paper reports 5th–95th posterior quantiles, synthetic tests, and an
EMPIAR-10278 application. It discusses pose accuracy, noise, path choice and
selection bias. These are direct predecessors to any proposed population UQ
method, not merely reconstruction baselines. The reported intervals are Bayesian
credible intervals under that model; the passages read do not establish uniform
frequentist coverage over arbitrary viewing distributions. The paper points to
public BioEM code but asks readers to contact the author for its MCMC code.

## Counting-particles study (2026)

We re-read the [primary published main text](https://doi.org/10.1038/s42003-026-09859-6)
and figure caption. This was already in the ledger. Evans et al. distinguish
hard/soft classification counts from solving the whole-stack population inverse
problem, with ensemble reweighting and RECOVAR deconvolution comparisons.
Their experimental 80/20 mixture is constructed from extensively classified
subsets and explicitly is not known biological truth. Supplementary derivations,
code and raw data were not inspected in this reading.

## Implication for method development (our interpretation)

Posterior population sampling or a two-state mixture likelihood is not a new
contribution. A worthwhile proposed advance needs a precise difference, such as
validated sensitivity to unknown viewing distributions or template errors, and
comparisons against the existing likelihood/reweighting methods. It must not
replace weak particle classification with a straw baseline while omitting those
methods. CAHRA's separate population/pose target is relevant, but does not itself
supply a new algorithm or establish coverage. We have not concluded that the
desired robust-population contribution is absent from the literature.
