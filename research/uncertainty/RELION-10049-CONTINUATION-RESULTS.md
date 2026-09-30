# Completed 10049 unknown-pose RELION continuation

30 September 2026. The declared initializer resumes iteration 40 to 100 in
4,547.06 seconds. Auto-refinement finishes in 921.03 seconds; its final
optimizer reports convergence and the exposure half labels are verified.
No settings, seed or checkpoint are selected by its FSC. This result is
weaker than the first stack and remains in the record.

Native half FSC crosses 0.5 at 22.73 A and 0.143 at 18.95 A, with mean shell
FSC .35685. RELION joins halves below 1/40 A, so agreement there is coupled.
After the fixed pilot-only alignment, RELION's mean shell agreement with
Gaussian, voxel and neural consensus maps is .20903, .20237 and .20126.
Registered-reference agreement is .20158 (original .07411), with 0.5/0.143
crossings 45.43/19.85 A. These reference/cross-method curves include interpolation
and frame-selection dependence; none is a density confidence coverage label.

The final local model diagnostics estimate rotational accuracy 8.87 degrees
and translational accuracy 4.579 A. They are not calibrated uniform radii.
The strong supplied-pose Gaussian/voxel/neural reference agreement on this
stack is therefore not reproduced by this declared small-cohort ab-initio
baseline. A convergence flag does not establish correct structure. Limited
particles, initial-model quality, local alignment and model/data mismatch are
possible explanations; this experiment does not distinguish them.

The existing evaluation runners are unchanged. Inputs, maps, alignment starts,
original/registered curves, local model diagnostics and exact source snapshots
are retained in `relion-evaluation-v2/10049` and
`reference-registered-comparison-v1/10049-relion.json`. No failed earlier
initializer record is replaced. The third-stack continuation is still running.
