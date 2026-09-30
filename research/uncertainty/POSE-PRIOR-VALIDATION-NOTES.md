# Targeted primary reading: pose correction and structural priors

The two primary XMLs are retained locally with hashes in the background manifest.
Reading covered the listed methods, validation and discussion sections; this is
not a claim to have reviewed all supplementary results. The curated ledger now
contained 100 candidates at this reading checkpoint, not 100 fully read papers.

CryoPROS (2025) uses generated auxiliary images during alignment and reconstructs
from experimental particles afterward. Its preferred-orientation experiments
include untilted/tilted HA-trimer accessions EMPIAR-10096/10097. Validation varies
reference models and uses noise-only reconstruction controls, local refinement
and map/model agreement. These are relevant model-bias tests, rather than a
per-particle confidence-region calibration. The results also separate poor
orientation coverage from erroneous pose assignment: neither should be treated
as a synonym for the other. A future uncertainty benchmark could hold out the
tilting intervention and vary initialization while checking whether proposed
uncertainty tracks the resulting failures. That is our proposed experiment,
not a guarantee made by CryoPROS.
Primary: https://doi.org/10.1038/s41467-025-59797-w

CoCoFold (2026) fine-tunes an AlphaFold structural model against particle images
through a Gaussian-mixture forward model. Its intended setting includes limited
particles and missing views; it uses supplied poses and CTFs. Reported static
Gaussian pose-perturbation controls show stability near a one-degree perturbation
scale and degradation at larger scales. Those controls make pose sensitivity
relevant to the present study, but they do not turn an assumed angle radius into
a simultaneous confidence set for experimental poses. The method targets atomic
structure, so scoring it as if it promised density-functional coverage would be
misleading. A matched uncertainty comparison would need the sequence/prior inputs,
coordinate target, upstream refinement and possible alternative conformations.
Primary: https://doi.org/10.1038/s42004-026-01899-7

These additions strengthen two open questions in the survey: how to validate
prior-assisted pose corrections against acquisition interventions, and how to
separate uncertainty from a structural prior from evidence supplied by the
particles. The current conditional audit does not resolve either question.

## Archived pose-summary availability

A separate source-file inventory checks all 105,247 / 108,544 / 131,899 metadata
rows in EMPIAR-10028/10049/10076. Both `alignments3D/pose_ess` and
`alignments3D/shift_ess` are identically zero in every file. These stored zeros
must not be interpreted as zero alignment uncertainty. The nonzero alignment
error fields do not have a verified conversion here to an angular confidence
radius. Thus the metadata inventory does not supply the independently calibrated
pose bounds requested by the reviewer. Exact metadata hashes and field summaries
are in `results/uncertainty/development/pose-metadata-inventory.json`; reproduce
with `scripts/audit_uq_pose_metadata.py` in a fresh output checkout.
