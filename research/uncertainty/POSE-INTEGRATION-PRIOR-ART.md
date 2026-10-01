# Pose integration: targeted primary reading

1 October 2026 UTC. [Brubaker, Punjani and Fleet (CVPR 2015)](https://www.cs.toronto.edu/~mbrubake/projects/CVPR15.pdf)
is newly added to the candidate ledger. Read sections 3.1, 3.3 and the
importance-sampling experiment in section 4; other sections were retrieved
but are not counted as a full reading. Figures and supplement were not
visually inspected. Author code is public but not executed here.

The paper marginalizes pose and shift, then importance-samples discrete
quadrature. Proposals use earlier likelihood evaluations, smoothing and a
uniform defensive component. Its footnote explicitly allows continuous
importance sampling. The numerical experiment compares log-likelihood
approximation with the complete discrete sum, which does not by itself
certify the underlying continuous integral. These are direct precedents for
our numerical gate; adaptive pose integration and defensive sampling are
not new contributions.

The local PDF is 7,790,960 bytes, SHA256
`dc5df4804c70c8d598cd0b3613cf8b3210891b9b203b859e0b27e87680abccb8`.
The author website supplied it successfully. We do not redistribute it.

Tyler (1987), DOI [10.1093/biomet/74.3.579](https://doi.org/10.1093/biomet/74.3.579),
is a bibliographic candidate for the angular central Gaussian distribution.
Publisher retrieval failed. Our protocol derives its sampling density by
radial integration and tests it independently; it does not claim a full
reading of Tyler or a novel distribution.
