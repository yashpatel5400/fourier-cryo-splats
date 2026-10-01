# cryoTWIN uncertainty and comparison protocol: supplementary reading

1 October 2026 UTC. Targeted reading of the [scientific supplement](https://static-content.springer.com/esm/art%3A10.1038%2Fs42004-026-02077-5/MediaObjects/42004_2026_2077_MOESM2_ESM.pdf)
(SI-3, SI-4 and selected GMM descriptions) deepens the existing ledger entry.
The first attachment was the published peer-review file; its distinct provenance
is retained. Neither source is our external Fable review. No full proof audit,
code execution or inspection of all figure images is claimed.

SI-4 varies training subsets and seeds across three fits, reporting variability
of correlations evaluated on training images. Its path experiment varies search
seeds on a fixed learned landscape and reports dispersion and energy standard
deviations. These are stability measurements, not repeated-sampling coverage
of population or density confidence intervals. SI-3 compares free-energy
correlations, using different density constructions across methods and a
lower-resolution RECOVAR condition because of memory limits. We must specify
matched resolutions and estimands when constructing our own comparisons,
rather than importing an ordering from this study.

For this project, uncertainty from optimization, finite sampling, viewing-law
error and template error need distinct evaluation. Adding an ensemble alone
would not answer that requirement. This reading closes the previously recorded
supplement-access gap for the specified sections; the overall source remains
only partially audited. Hashes and exact reading scope are in
`cryotwin-supplement-source.json` and `cryotwin-peer-review-source.json`.
