# Gaussian representation and validation scope

Targeted reading, 30 September 2026, of [GaussianEM version 2](https://arxiv.org/html/2512.21599v2),
sections 2.1–2.2, 2.8 and 3–5. The January revision supersedes the December title.

GaussianEM learns density changes and displacements of real-space Gaussians.
The fitted consensus map supplies initialization and poses. Evaluation includes
conformational subspace overlap and cluster-reconstruction FSC. Its discussion
acknowledges missing initial regions and possible pose inconsistencies.

Our assessment: these checks concern heterogeneity recovery, not confidence
coverage for a density functional. The initialization limitation supports
checking errors outside a fitted representation, without establishing the
usefulness of our proposed bounds. A matched heterogeneity comparison would
require an appropriate new target and experiment.

The 36-page PDF is retained locally with its hash. Figures, supplementary
analyses and original code have not been independently checked.
