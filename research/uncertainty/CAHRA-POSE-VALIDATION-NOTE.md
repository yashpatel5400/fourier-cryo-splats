# CAHRA pose-convention data and access check

1 October 2026 UTC. The live [pose-entanglement challenge page](https://heterogeneity.notion.site/pose-entanglement)
was inspected in the browser, including the lower-page attachments not exposed
in the initial web text extraction. Its public supplement was downloaded and
all ZIP member CRC checks passed. It contains two 216-image, 256×256 float32
stacks (open/closed), matching STAR metadata and angle-labelled PNGs. Both
stacks have 1.2 Angstrom pixels and zero supplied shifts. Each Euler angle
uses the six values 0,30,60,90,120,150 degrees.

The accompanying convention note specifies cryoJAX-generated poses exported
to RELION conventions, a z-aligned C11 symmetry axis, and a symmetry-minimized
rotation-angle error. The supplement is described as high-signal, without CTF;
it is a coordinate-validation resource, not the main noisy challenge or a
population coverage test. No templates or poses have been fitted here. A subsequent descriptive in-plane
coordinate check is recorded below.

The main collection link redirects to Globus login in this browser. No main
8 GB stack or template maps were downloaded, and no account was created.
The live public pose/state leaderboard was readable, but contains submission
identifiers and ordinal ranks, not enough information to identify algorithms
or infer uncertainty calibration. Its existence supplements the earlier
preprint reading, which describes dataset construction without full outcomes.

Download provenance, archive and member sizes, extracted-file hashes, MRC
headers and STAR column inventories are recorded in
`cahra-pose-validation-download.json` and `cahra-pose-validation-inventory.json`.
The external ZIP remains local; source URLs are linked rather than silently
redistributing a benchmark under an unverified license.

## Subsequent in-plane coordinate check

For all 108 available +90-degree Psi pairs in each state, compared both
quarter-turn signs using the DFT origin at pixel 128, retaining every pair.
No scale, shift or filtering was fitted. The positive NumPy quarter turn,
with its explicit one-pixel center correction, has median centered correlation
.98777 / .98776 and relative L2 error .15575 / .15584 (closed/open).
The opposite sign has correlations .73713 / .74491 and errors .72227 / .71132.
The residual is not zero: the supplied images are described as high-signal,
not noiseless. This identifies the in-plane array convention but does not
validate a complete 3D projector, recover unknown poses, estimate a physical
noise law or establish population uncertainty. Both signs and all 432 rows are
saved in `results/uncertainty/development/cahra-inplane-convention-v1/`.
