# Alignment, collective recovery and reference bias are different questions

1 October 2026 UTC. This post-review addendum checks the references suggested
by full review 3 and follows two recent primary citations. It does not replace
the frozen reviewer packet or imply a fourth favorable review. Source/access
details and local download hashes are recorded separately.

**Perry et al. (2019).** [Published metadata](https://doi.org/10.1137/18M1214317)
and [author manuscript](https://arxiv.org/pdf/1707.00943), sections 1–4 and
Theorems 1/3 read selectively. In cyclic MRA, generic signals with suitable
nonvanishing Fourier coefficients can be recovered collectively at low SNR,
with sample complexity scaling as SNR^-3. Their model distinguishes estimating
the signal from assigning a latent shift to each observation; poor individual
alignment does not prove signal nonidentifiability. Power spectra alone lose
phase, whereas third-order invariants can identify generic cyclic signals.
These results concern the stated finite cyclic model, not an arbitrary
continuous 3-D viewing law, experimental covariance, or finite-sample coverage
for local cryo-EM density. Our inference is that the failed local optimizer and
quadratic moment tests cannot support a general impossibility claim. Bispectrum
features are established prior art, so adding them alone is not novelty.

**Henderson (2013).** [Primary author abstract and figures](https://pubmed.ncbi.nlm.nih.gov/24106306/)
and the earlier accessible [PMC text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3831464/)
discuss selection of noise resembling a template, misleading reference-driven
reconstructions, and failures of apparent validation under shared processing.
This directly supports care with particle selection and half-map independence.
It is not a proof of the particular conditional-centering failure or attenuation
law in our fixed-data experiment. Current PMC requests triggered verification;
the Europe PMC full-text request returned HTTP 500, which is preserved as a
failed JSON response rather than described as a retrieved article.

**Neyman and Scott (1948).** The [official issue contents](https://www.jstor.org/stable/i332816)
verify the title *Consistent Estimates Based on Partially Consistent
Observations*, authors, volume 16(1), pages 1–32, and DOI 10.2307/1914288.
Full-text access was unsuccessful. This is a historical bibliographic lead;
we have not checked its original proof or proved that our cryo-EM optimizer
satisfies a Neyman–Scott asymptotic model. The review's suggested connection
must not collapse incidental parameters, alignment attenuation, and template
selection into a single established theorem.

**Janson and Andén (ICASSP 2026; preprint v1, 2025).**
[Primary manuscript](https://arxiv.org/pdf/2510.12651v1), sections 1–6 read;
[conference listing](https://cmsworkshops.com/ICASSP2026/view_paper.php?PaperNum=14410)
and publisher DOI 10.1109/ICASSP55912.2026.11463714 verify publication context.
Their moment-based posterior sampler combines a sample-power likelihood with
a learned diffusion prior. The likelihood is noncentral chi-squared under the
cyclic Gaussian model; the diffusion conditioning uses explicit approximations.
The authors acknowledge those approximations and the possibility of nonvanishing
data-limit error with power-only conditioning. Experiments compare aligned
point-estimation error on two synthetic 41-coordinate signal families against
EM and bispectrum inversion. Those experiments are not an experimental cryo-EM
density-coverage benchmark. This is direct prior art for posterior uncertainty
from invariant summaries and learned priors; our inference is that visual
posterior spread alone would not settle calibration under prior misspecification.
We have not run their implementation.

**Al-Ghattas et al. (2026).** [Published metadata](https://doi.org/10.1137/25M1765602)
and [v2, 18 May 2026](https://arxiv.org/pdf/2506.12201v2), introduction and
sections 2.1–2.4, including the identification theorem/proof, read selectively.
Their functional MRA model uses compactly supported signals and centered,
compactly supported noncyclic shifts with a specified Gaussian-process
covariance. A derivative of the full second-order Fourier moment identifies
the signal through a Kotlarski-type formula under nonvanishing conditions.
This is not an invariant diagonal-power statistic: retaining off-diagonal
derivative information and centering assumptions changes the model. It does
not contradict our characterization of bilinear statistics invariant to every
common shift. The paper's broader recovery/zero-handling proofs and MATLAB
implementation have not been audited. No direct transfer to unknown 3-D
projection orientations or empirical covariance follows from the passages read.

Baldwin–Penczek (2005) and Park et al. (2011), previously recorded in the
[alignment-moment note](ALIGNMENT-MOMENT-FOLLOWUP.md), are now also entered in
the candidate ledger with their original reading-depth qualifications. The
ledger increases from 131 to 138 candidates; this is not 138 full readings.

## Consequence for the next experiment

The declared [bispectrum gate](paired-power-v1/BISPECTRUM-PROTOCOL.md) tests a
specific restriction left by the quadratic construction. It preserves common
translation invariance but still needs continuous rotation control, specified
noise and amplitude assumptions. Exact conditional noise variance helps assess
whether extra discrimination survives its statistical cost. Sampled separation,
even with correct variance algebra, is not a calibrated continuous-pose test.
Experimental processing, CTF variation, heterogeneity and valid source-group
splits remain separate requirements. A new paper contribution must establish
usefulness beyond the existing moment and posterior methods, not relabel them.
