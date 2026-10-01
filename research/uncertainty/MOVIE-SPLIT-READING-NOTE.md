# Movie-frame splitting: prior art and a remaining independence question

30 September 2026, after full-review-2 packet submission. This is a targeted
reading and mathematical planning note, not an implemented experiment or a
claim that movie splitting is new. The 119-candidate submitted survey remains
unchanged pending the next revision.

[Bepler et al., Topaz-Denoise (2020)](https://www.nature.com/articles/s41467-020-18952-1)
already use odd/even movie-frame pairs for Noise2Noise training. Figure 1
describes separate processing and summation of the halves. The training-data
methods include direct sums and MotionCor2-aligned pairs, and explicitly note
that motion correction can change the noise distribution. The source also
reports that its selected denoised-particle ab-initio trials were less reliable
than the corresponding raw-particle trials. These are denoising and
reconstruction observations, not density confidence guarantees.

Reading scope: primary abstract, introduction, Figure 1 caption, selected Results
passages, and the indexed primary training-data preparation methods. A later
direct Nature fetch failed through its identity redirect; the same methods
were available from the indexed [PMC version](https://pmc.ncbi.nlm.nih.gov/articles/PMC7567117/).
No Topaz code was executed, supplementary material was not inspected, and no
new raw movies were downloaded.

## Our conditional argument, not a claim from that paper

Let frame groups A and B have independent noise conditional on fixed physical
signal/trajectory parameters. Let every selection, pose, interpolation and
weight-design decision D be measurable using A and independent training data.
If held-out raw noise e_B is Gaussian with covariance Sigma_B, processing B
only by a D-fixed linear operator T_D gives

    T_D e_B | A ~ N(0, T_D Sigma_B T_D^T).

This is an elementary conditioning/linear-transformation identity. It explains
which design dependence a split could remove. Merely producing two separately
motion-corrected sums does not supply this exact conditional model if processing
of B itself fits transformations to B's noise. Likewise, preprocessing using
the full movie before splitting can reintroduce cross-half dependence.

This does not calibrate pose radii, a pilot-centered density class, Gaussianity,
dose-dependent signal changes or covariance stability. Fixed detector artifacts
and frame correlations also need evaluation. Picking on A avoids direct
selection on B noise under the stated conditional-independence premise, but
false picks and heterogeneous particles can still violate a shared-density
model. A future protocol would need a fully specified processing graph, raw
frame access, treatment of motion/dose/CTF uncertainty and realistic negative
controls. It must not replace the present unresolved input conditions with
an unsupported independence assertion. No next experiment is selected here.

## Cached source metadata feasibility

The existing 10028 archive metadata lists 1,081 unaligned 16-frame, 4096-square
float32 movies (about 1 GiB of pixel values per movie, excluding headers). The
listed image sets on 10049 and 10076 contain one frame each, even where their
category describes a multiframe origin. Thus an odd/even raw-frame procedure
cannot simply be retrofitted onto the three downloaded particle stacks. The
10028 entry is a possible source to inspect further, not a selected experiment.
The exact cached metadata hashes and arithmetic are recorded in
`provenance/uncertainty/movie-split-metadata-feasibility.json`. No current remote
file enumeration, source-file availability check or movie download has occurred.
