# Deposited-map frame audit and independent-pilot registration

30 September 2026. Post-outcome development, motivated by the first completed
RELION continuation evaluation. On 10028, native half agreement and agreement
with all three existing reconstruction methods are much stronger than every
method's unregistered deposited-map agreement. This motivates checking the
unverified assumption that consensus-pose and deposited-map frames coincide.
It does not prove a registration problem or justify selecting a favorable
reference transform from inference maps. Preserve all original comparisons.

Phase A: for each of the same three geometries, reconstruct the existing old
pilot from its saved Gaussian coefficients exactly as the RELION evaluator
does: 64-cube, support radius 0.35 of the field, physical pixel spacing.
Register the corresponding deposited map **to this pilot only** using the
already tested `align_map_to_pilot` implementation. Keep its 40 A lowpass,
32-cube working grid, all 24 proper principal-axis starts in each of two hands,
300 evaluations per start, and linear interpolation of the final full map.
Select by pilot correlation only. Save all starts, the selected rigid transform
and hand, both input/output hashes and the resampled reference. This is a
local alignment heuristic with a 48-start search, not a global-optimum claim.
No learned inference-half map, experimental interval, or FSC selects the
transform. Run all three geometries. The fixed work budget is 48 starts per
stack, no post-outcome restart/extra seed. Record failures.

Phase B: after all three transforms are committed, apply each frozen reference
to the already locked Gaussian, voxel and neural averages, and to completed
RELION runs as they become available. Use the same transformation for every
method within a stack. Retain original unregistered reference curves as well
as registered curves, the sampling-limit censoring and the pilot dependence.
Never label reference FSC as independent reconstruction truth or density
coverage. A poor 48-start pilot fit remains poor; no reference-selected fitting,
checkpoint, pose, reconstruction or target is allowed.

Historical interval comparisons and simulation generators are not overwritten
or reinterpreted as using registered references. A later target comparison, if
needed, requires its own transparent development record. This frame audit does
not calibrate experimental noise/pose classes or the homogeneous model on 10076.
