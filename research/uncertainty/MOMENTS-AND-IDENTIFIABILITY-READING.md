# Moments, latent orientations and validation: targeted readings

30 September 2026. Four additional primary works, not four full replications.
PDFs are locally hashed; only bibliographic records and these notes are published.
This reading addresses an important limitation of our known-map information
diagnostic: difficulty estimating one particle's orientation does not establish
that the molecule cannot be recovered collectively from many particles.

**Sharon et al. (2020), nonuniform method of moments.**
[Published article](https://doi.org/10.1088/1361-6420/ab6139);
[arXiv v2](https://arxiv.org/abs/1907.05377v2), sections 3 and 5 read.
Theorem 1 gives unique recovery with a known, totally nonuniform orientation
law under explicit rank, distinct-eigenvalue and bandlimit conditions.
The unknown-law cases report numerical Jacobian-rank evidence for generic
finite recovery modulo rotations over specified expansion sizes. These are
not uniform uniqueness guarantees for arbitrary molecules or viewing laws;
the paper explicitly distinguishes numerical rank evidence from exact proofs.
It also reports conditioning. Its discussion retains missing CTF/centering
effects and synthetic-data scope. Nonuniformity can improve particular
low-order identifiability questions without making preferred orientations
universally beneficial for experimental resolution.

**Hoskins et al. (2026), SubspaceMoM.**
[Published article](https://doi.org/10.1137/24M1699644), online 27 February 2026;
[arXiv v5](https://arxiv.org/abs/2410.06889v5), November 2025.
Read: introduction, section 2.4, the reduced formulation in section 3, and
discussion. Low-rank compression and quadrature reduce the cost of first,
second and third moments. The method jointly fits a volume and viewing law,
without assigning a pose to each particle. Its principal reconstruction
experiments are synthetic, with CTFs and specified white noise. Centering,
selection, heterogeneity, detector response and high-resolution memory use
remain material limitations. Compression need not preserve the full moment
objective when the retained range is incomplete. This is relevant prior art
for both scalable reconstruction and any proposed confidence calculation
after statistic compression. We have not run its author implementation.

**Fan et al. (2024), high-noise orbit likelihood.**
[Published article](https://doi.org/10.1214/23-AOS2292);
[arXiv v2](https://arxiv.org/abs/2107.01305v2).
Read: model, Remark 2.1, Theorem 2.7 and surrounding interpretation, not all
appendix proofs. The model uses Haar-uniform independent rotations and Gaussian
noise. It relates high-noise Fisher eigenvalue scales to degrees of invariant
polynomials at generic signals. The projected result has a finite-orbit
identifiability condition; mirror ambiguity and other projection ambiguities
are distinguished. Its reported protein information simulations omit
tomographic projection. This theory concerns a marginalized signal likelihood,
not our per-particle known-map pose information. Neither calculation supplies
finite-sample experimental pose-confidence radii.

**Zhang et al. (2024), moment-based molecular metrics.**
[Published article](https://doi.org/10.1017/S2633903X24000023);
[author preprint](https://arxiv.org/abs/2401.15183v1).
Read: section 3.2, experimental protocol 4.4, and limitations 5.
The image metric compares moment estimates with a proposed structure while
fitting the unknown viewing-law expansion by linear least squares; positivity
of that expansion is relaxed. The experimental example separates five
EMPIAR-10076 classes, uses deposited shifts and rescales moments. It is not
an independent validation of our pooled homogeneous interpretation of 10076.
The paper retains limited discrimination among similar states and sensitivity
to unmodeled acquisition effects. No finite-sample local-density coverage
claim should be inferred from its structural ranking results.

**Implication for this project.** Our worst-case local pose audit and its
negative power results are properties of that procedure and declared class.
They are not an impossibility result for latent-pose reconstruction. A possible
alternative is calibrated structural validation using marginalized likelihoods
or moments while profiling the viewing law. Making that useful would require
controlling acquisition/model error, compression and optimization error,
especially on experimental particles. It is not an implemented contribution
in the present paper, and the general statistical constructions already exist.
