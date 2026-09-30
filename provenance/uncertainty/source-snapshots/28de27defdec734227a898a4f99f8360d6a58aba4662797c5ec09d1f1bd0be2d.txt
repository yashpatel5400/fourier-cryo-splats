# Original CryoLike CPU comparison on three experimental stacks

30 September 2026. Declared after the initial continuous-denominator failures
and after verifying the author's CPU cross-correlation test, before any scores
from this comparison. This is exploratory reuse of an existing prediction
cohort, not fresh confirmation or uncertainty calibration.

Use the unchanged author repository at commit
`a413ffd265e2815c9f81c492a3d8bfedaf35d737` from
https://github.com/flatironinstitute/CryoLike. The convenience file runner hardcodes
CUDA; call its original CPU-capable `template_first_comparator` and
`compute_optimal_pose` instead. Preserve the external source commit and hashes.
The only local adapter changes array axes, passes explicit physical units and
known CTF values, prepares inputs, and enforces a wall budget. It does not
replace the author's Fourier--Bessel comparison or integrated-score kernel.

For each of EMPIAR 10028, 10049 and 10076, choose the first row for each source
exposure in the previously frozen prediction-v1 selection CSV, sort by source
exposure string, and take up to the first 128 distinct exposures (all if fewer are available).
Metadata inspection before scoring finds 115, 71 and 185 available exposures,
so the three selected counts are 115, 71 and 128. This deterministic
selection is fixed before the new scores. Read the corresponding stored 64x64
raw-particle arrays and CTFs. Do not use deposited orientations or translations
in scoring. These images have already been examined in the earlier prediction
study, and deposited CTFs/preprocessing remain external assumptions.

Compare four volumes per geometry: the mean of the two locked Fourier Gaussian
reconstructions, the mean of the two locked voxel reconstructions, the mean of
the two locked stock cryoDRGN reconstructions at their previously frozen epochs,
and the deposited EMDB map resampled to the same field and 64-cube. The deposited
map is an approximate comparison structure, not truth and not independent of
these public particles. All learned reconstructions used supplied training
poses; this comparison does not make them ab initio.

Use physical [x,y] images and [x,y,z] volumes as required by CryoLike, explicitly
transposing the project's [y,x] and [z,y,x] arrays. Independently test both
transforms against asymmetric direct Fourier sums before any biological score.
Supply physical field widths on returned template objects: the low-level
volume constructor otherwise leaves a normalized default box, which would
misinterpret angstrom displacements. Retain author centering and normalization:
zero-mean/max-norm physical images, then unit L2 polar images, and unit L2
polar templates. Template generation retains the author's DC-mask subtraction.
These normalized scores are not raw-pixel normalized predictive densities.

The maximum box Fourier radius is 12 (CryoLike polar radius 6), radial spacing
0.25 in its units, and 128 in-plane angles. Use both viewing separations 0.4
and 0.2 radians (173 and 470 viewing directions in the pinned version), without
selecting a preferred setting afterward. Search a 5x5 displacement grid from
-2 to +2 downsampled pixels per axis, retaining zero. Use double precision,
NUFFT tolerance 1e-12, CPU torch threads 2, template batches 8 and image batches
16. The CTF has the same project's deposited convention, sampled at physical
frequencies 2*(polar_x,polar_y)/field_A. The radius corresponds to physical band
endpoints about 40.2, 19.68 and 34.93 angstrom, not atomic-resolution validation.

Run all 24 map/grid cases. Budget each score computation at 900 seconds, checked
between original API batches; setup and final saving are separately timed.
Save a failed or limited case and continue the declared remaining cases. Record
template counts, all per-image maximum correlations and integrated log scores,
optimal searched poses/shifts, view/quadrature coordinates, source/input hashes,
all statuses, elapsed times and peak resident memory. Do not interpret any
local search as a global maximum on continuous SO(3).

For each grid and geometry, report paired mean differences from the neural map
for Gaussian, voxel and deposited reference, separately for maximum correlation
and integrated log score. Every selected exposure contributes (115, 71 and 128 by geometry). Use 20,000
paired exposure bootstrap resamples with seed 940001+dataset+1000*grid_index,
preserving all method scores. Report marginal 95% percentile intervals and
Bonferroni intervals within each metric's eighteen-contrast family. These
bootstrap intervals are descriptive approximations, not finite-sample density
or model-validity guarantees. Retain both grids and unfavorable comparisons.

No likelihood-ratio rejection threshold or e-value is assigned to CryoLike's
integrated score: its nuisance marginalization/normalization conventions remain
distinct from the candidate universal-inference construction. Differences in
noise, amplitude and pose treatment prevent claiming this is a matched-power
comparison. This study adds original-code structural scoring and computational
evidence while the proposed continuous test remains unresolved.
