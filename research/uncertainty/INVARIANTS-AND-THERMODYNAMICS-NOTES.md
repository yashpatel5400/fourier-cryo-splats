# Additional targeted primary readings

These notes deepen the existing orbit-recovery entry and add one recent
thermodynamics paper. They are not full-paper-plus-supplement replications.

## Moment-based posterior sampling

[Janson and Andén](https://arxiv.org/html/2510.12651v1) condition a diffusion prior
on the power spectrum of randomly shifted one-dimensional signals. The Gaussian
observation model gives a noncentral chi-square likelihood for these moments.
Power spectra omit phase information; the authors explicitly distinguish prior
concentration from consistency as data increase. Their sampler also makes
approximations to intermediate diffusion likelihoods. This is relevant to
representation/identifiability audits, but is not an experimental cryo-EM
reconstruction baseline. For our project, a future comparison could separate
uncertainty lost by statistic compression from uncertainty suppressed by a
learned prior. Read: sections 1–4; the existing five-page PDF is archived
and hashed in the background manifest.

## Thermodynamic inference and the 10076 stack

[Yamazaki et al.](https://doi.org/10.1038/s42004-026-02077-5) propose cryoTWIN/PaStEL,
using an isometric latent representation to estimate conformational free-energy
landscapes and pathways. Their assumptions include adequate denoising, known
poses, unbiased sampling and equilibrium. They evaluate heterogeneous states
on EMPIAR-10076 and discuss finite-data, pose and particle-selection biases.
The main text points to supplementary uncertainty analyses, which we have not
read. This sharpens an open question: calibrating population or energy uncertainty
requires more than reconstructing one consensus density. Our current shared-
density audit does not answer it. Read: theory assumptions, 10076 experiment
and discussion; public XML hash recorded separately.
