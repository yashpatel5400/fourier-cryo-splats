# Post-review-4 stack nuisance inventory

1 October 2026 UTC. Descriptive analysis of existing metadata; no new score,
test calibration or method-selection experiment. Field names, counts and ten
pose-convention comparisons per stack were inspected before this declaration.
No histogram or full-field summary was inspected before freezing this protocol.

Use all records in the three archived cryoSPARC files and the author-provided
`poses.pkl`: 105,247 / 108,544 / 131,899 rows for 10028 / 10049 / 10076.
Also report the published filtered subsets where supplied. Retain the original
half labels. These are metadata populations, not the downloaded 8,192-particle
subsets and not new biological replicates. Validate rotation orthogonality and
the transpose of the cryoSPARC rotation-vector convention against the stored
matrices across all records before interpreting angles.

For full SO(3), use equal-Haar-volume ZYZ bins: alpha and gamma uniform on
[0,2 pi), cos(beta) on [-1,1]. Report all b=4,8,12 grids, b^3 bins. Also report
unoriented detector-plane normals (row 3 of the stored rotation) in the upper
hemisphere, with 2b azimuth bins and b uniform-z bins. This removes antipodal
equivalence but does not assert or quotient an unverified molecular symmetry.

For each grid/subset, retain every count by recorded half, maximum empirical
density ratio to uniform, plug-in chi-square divergence, the usual iid
multinomial finite-count correction, cross-half inner-product divergence,
half-to-half total variation and empty-bin fraction. The iid corrections are
descriptive references: shared processing and pose fitting can invalidate their
sampling interpretation. No continuous-density upper bound, pose-error bound
or experimental kappa certificate follows. Orientation estimates, symmetry
representatives and coordinate frames can affect these diagnostics.

Summarize all numeric CTF/3-D-alignment scalar fields using finite/nonzero
counts, extrema and 1/5/25/50/75/95/99 percentiles. Do not equate constant stored
CTF scales with known constant physical amplitudes, or residual powers with
pure-noise variances. For alpha, weight and CTF scale, report weighted between-
bin fraction of variation using the b=8 normal bins only where the field varies;
this is association, not a calibrated physical nuisance model. Report defocus
and astigmatism in Angstrom, and all unique scalar instrument settings. Source
noise spectra are not present in these metadata; this inventory cannot supply
them. Later image-based noise diagnostics must have their own stated limits.

Record input hashes, code revision, all failures and numerical discrepancies.
The result can falsify the plausibility of a nearly uniform *recorded* pose
distribution. It cannot establish the true law needed by the earlier theorem.
