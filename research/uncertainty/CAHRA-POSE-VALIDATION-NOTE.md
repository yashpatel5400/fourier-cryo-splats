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
population coverage test. No images have been fitted or scored here.

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
