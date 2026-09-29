# Research status

This is an experimental research project, not a validated reconstruction package.

Requested work: method development, ICML-format manuscript, literature collection,
real-data fitting on at least three cryo-EM datasets, matched baseline/FSC analysis,
and publication to a Git remote.

Hardware: Apple M4 Pro, 24 GiB unified memory, MPS available, initially 67 GiB free.
GitHub CLI is installed but was unauthenticated when the project began.

Initial study: real, deposited extracted particle images from EMPIAR-10028,
EMPIAR-10049 and EMPIAR-10076, with supplied cryoDRGN consensus poses and CTFs.
These are particle stacks, not unprocessed detector movies. Known-pose experiments
do not establish ab initio performance. Consensus-pose half-map FSC is conditional
on those poses and must not be called fully gold-standard FSC.

No numerical results were available when this file was first created.
