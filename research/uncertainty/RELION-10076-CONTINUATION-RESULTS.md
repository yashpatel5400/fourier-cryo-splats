# Completed 10076 unknown-pose RELION continuation

30 September 2026. The declared initializer resumes iteration 30 to 100 in
5,425.20 seconds. Auto-refinement finishes in 7,144.11 seconds and reports
convergence; the exposure half labels are verified. The historical timed-out
initializer remains, and no settings or checkpoint are selected by FSC.

Native half FSC crosses .5 at 24.73 A and .143 at 15.92 A, with mean .53759.
RELION joins halves below 1/40 A, where their agreement is coupled by fitting.
Registered-reference mean FSC is only .12169 (original .06678), with .5/.143
crossings at 101.87/47.64 A. Cross-method mean agreement with Gaussian, voxel
and neural consensus maps is .33154/.32522/.33091. All curves are retained.

EMD-8434 is a Class A assembly intermediate, not ground truth for this
heterogeneous population. Native half agreement does not establish the correct
state or remove reconstruction, registration, preprocessing or selection error.
The local model reports estimated rotational accuracy 2.715 degrees and
translation accuracy 2.9082 A; these are not calibrated uniform pose bounds.

All three declared unknown-pose CPU continuations and original/registered
comparisons are now complete. The Gaussian/voxel/neural baselines instead use
supplied poses, so these are different reconstruction tasks. The all-case FSC
figure preserves all twelve method/stack pairs and both frames. Neither FSC
nor a convergence flag is a density confidence coverage measure.

Authoritative outputs are `relion-evaluation-v2/10076/metrics.json` and
`reference-registered-comparison-v1/10076-relion.json`; the fit record's
historical `evaluation_complete: false` predates these separate evaluations.
