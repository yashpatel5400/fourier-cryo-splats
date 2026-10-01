# Published orientation-prior discussion: targeted follow-up

1 October 2026 UTC. Read Appendix F of the published
[Bayesian orientation paper](https://doi.org/10.1107/S2059798326001415),
an existing ledger entry. This is not a new paper count or a full reading.
PMC HTML presented a browser check; the public Europe PMC full-text XML
was retrieved successfully (277,417 bytes, SHA256
`ac4f61f24ce9ee05eaa64f676bca0547d8bcf1bc0ac1a650bf4a91529f3ce2d0`).
An initial local parser import failed; the standard-library XML parser then
extracted Appendix F. No figure, supplement or author code was examined.

Appendix F.2 distinguishes joint structure/viewing-law estimation in
likelihood refinement from moment-based joint recovery. It discusses
smoothed priors from an initial refinement and sensitivity to prior choice,
while leaving principled prior estimation and robustness theory open.
This supports investigating viewing-law sensitivity, but it does not establish
a gap for simply adding a learned prior: that direction is already explicit.
An empirical pose histogram is still an estimate, and using it as a fixed
known prior would omit an uncertainty source central to this project.
