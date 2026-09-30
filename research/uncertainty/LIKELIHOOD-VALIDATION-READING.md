# Likelihood validation and inference: targeted primary readings

30 September 2026. This deepens the existing Ortiz entry and adds three
likelihood-validation/statistical candidates to the earlier ledger.
Retrieval and targeted reading are distinguished below; none is a reproduced
baseline or a full-paper-plus-supplement review.

**BioEM (Cossio and Hummer, 2013).** The primary indexed methods describe
per-image evidence for structural ensembles, integrating nuisance parameters
and allowing ensemble weights. Their model includes orientation, displacement,
intensity, offset and noise. This is direct prior art for validating structures
against individual particles, rather than only comparing reconstructed maps.
Read: indexed sections 2.1 and 4. PMC live access returned a browser challenge;
the Europe PMC XML attempt returned 500. The metadata DOI is
[10.1016/j.jsb.2013.10.006](https://doi.org/10.1016/j.jsb.2013.10.006).

**Independent-particle validation (Ortiz et al., 2020).**
[Primary article](https://doi.org/10.1016/j.yjsbx.2020.100032), XML locally hashed.
Read: methods, the validation protocol and conclusions. The paper compares
held-out map evidence and a normalized Jensen–Shannon statistic across frequency
cutoffs, including EMPIAR-10049. BioEM first searches 36,864 orientations, then
uses 1,250 candidates near ten selected orientations. Its noise/alignment
controls and caveats about picking, classification and heterogeneity matter
for any proposed image-based validation. These diagnostics do not constitute
a stated frequentist confidence set for an arbitrary local density feature.
Figures and supplement were not independently inspected in this pass.

**CryoLike (Tang et al., 2025).** The
[published abstract and figure 1](https://doi.org/10.1107/S2059798325009350)
and [author repository](https://github.com/flatironinstitute/CryoLike) were checked.
It accelerates image-to-structure likelihood calculation with Fourier--Bessel
representations. Its workflow includes CTFs, templates, correlations, optimized
image parameters and likelihood outputs. This makes it a relevant computational
baseline for particle-based validation. The main publisher PDF returned 403;
the main methods have not yet been read. Do not infer calibrated
feature intervals or a specific numerical-integration guarantee from the abstract.
The pinned author code's mathematical-framework and likelihood-computation
documentation, plus its integrated-likelihood kernel, were also inspected.
The documentation itself labels parts incomplete; normalizing constants and
integration versus profiling would need auditing before reuse in a calibrated
test. No author code was executed.

**CryoLike supplement follow-up, 30 September 2026.** The web text of the
[supplement](https://journals.iucr.org/d/issues/2025/12/00/bar5002/bar5002sup1.pdf)
became accessible. Targeted reading covers S1.12--S1.15, not the entire supplement.
The image model has independent Gaussian pixel noise, intensity and offset.
Equation 52 integrates intensity/offset over the real line; equation 55
integrates the noise standard deviation over the positive line. Equation 58
uses exponent (N_pixels-3)/2. Equation 60 approximates orientation integration
with quadrature for a specified, resolved viewing density. S1.15 relates the
normalized-image alignment objective to cross-correlation. Direct PDF retrieval
still returned 403, and two web screenshot attempts timed out: no local PDF
hash or visual inspection is claimed.

Our inference: unbounded Lebesgue nuisance integrals do not by themselves give
a normalized density on raw image space. Tonelli's theorem integrates each
conditional density to one, leaving infinite nuisance measure. Reusing such a
score as a universal-inference numerator therefore needs proper priors or an
explicit alternative sample-space/normalization argument. This does not impugn
its use for the paper's likelihood comparisons. The pinned code's kernel uses
exponent N_pixels/2-2; the convention behind its difference from equation 58
remains unresolved, not an established implementation defect. The code commit
is a413ffd265e2815c9f81c492a3d8bfedaf35d737; the earlier manifest hashes its kernel.

**Universal inference (Wasserman, Ramdas and Balakrishnan, 2020).**
[Published paper](https://doi.org/10.1073/pnas.1922664117);
[arXiv v4](https://arxiv.org/html/1912.11436v4) sections 2 and 6 read.
The split likelihood-ratio construction, nuisance profiling and replacing a
null likelihood maximum by an upper bound are established methods. They can
provide finite-sample testing under a correctly specified likelihood, but do
not make the noise law correct or a nonconvex local fit a global upper bound.
Using these ingredients in cryo-EM would require a substantive computational
and empirical contribution; their general coverage argument is not new.

**Open question suggested by these readings.** Can continuous pose searches be
bounded tightly enough to turn held-out image scores into useful tests of
specified structural hypotheses, while retaining a calibrated noise model?
This might avoid externally supplied small pose balls, but it changes the
inference target and still faces model error, dependence and computational
limits. It is a research direction, not an implemented result in this paper.


## 30 September: original-code CPU experiment completed

The subsequent protocol executes the unchanged author comparator and
integrated-score kernels on all 24 prescribed real-particle/map/grid cases.
See CRYOLIKE-BASELINE-PROTOCOL.md and CRYOLIKE-BASELINE-RESULTS.md, the pinned
commit, adapter tests and saved per-image arrays. Earlier statements that
original code had not yet been executed describe the reading-stage status;
this dated follow-up supersedes that status without changing the interpretation
of the score's normalization or claiming a matched e-value baseline.
