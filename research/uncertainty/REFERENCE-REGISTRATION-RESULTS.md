# Reference frame sensitivity and first completed unknown-pose baseline

30 September 2026. All three phase-A pilot-only registrations complete in
17.28, 16.76 and 16.23 seconds. Each retains all 48 starts. Their transforms
and hashes were committed at `becc3a6` before any phase-B FSC computation.
No inference map, interval, reference-FSC score or reconstruction checkpoint
selects these transforms. The protocol is explicitly post-outcome development.

Ten phase-B cases complete: Gaussian, voxel and neural on all three stacks,
plus the completed RELION continuation on 10028. The remaining RELION fits are
still active; no outcomes are inferred for them. Every original reference
curve remains alongside the newly registered curve.

| Stack | Method | Original mean shell FSC | Registered mean shell FSC | Registered 0.5 crossing (A) |
|---|---|---:|---:|---:|
| 10028 | Gaussian | 0.0704 | 0.7269 | 19.01 |
| 10028 | Voxel | 0.0706 | 0.7066 | 20.40 |
| 10028 | Neural | 0.0716 | 0.7364 | 18.73 |
| 10028 | RELION | 0.0748 | 0.7024 | 21.32 |
| 10049 | Gaussian | 0.0524 | 0.7566 | 8.58 |
| 10049 | Voxel | 0.0499 | 0.6918 | 9.89 |
| 10049 | Neural | 0.0511 | 0.7994 | censored at 7.872 |
| 10076 | Gaussian | 0.0858 | 0.1268 | 95.73 |
| 10076 | Voxel | 0.0861 | 0.1271 | 96.45 |
| 10076 | Neural | 0.0869 | 0.1273 | 94.27 |

Every registered 0.143 comparison on 10028/10049 is censored at the sampled
limit, 16.08/7.872 A. Those are not measured resolution values. On 10076 the
0.143 crossings remain 51.52--52.28 A. A mean over shells is descriptive and
does not weight shells by particle information or voxel count. These results
confirm an important frame effect on the first two stacks; they do not establish
a general registration optimum or independent reconstruction accuracy.

The native 10028 RELION half-map 0.143 curve is censored at 16.08 A; its 0.5
crossing is 21.21 A. Half labels are verified at the exposure level. RELION joins
halves below 1/40 A, and the source pool/preprocessing remain shared with other
methods. Its continuation takes 5,492.27 seconds for the remaining initializer
and 4,012.07 seconds for auto-refinement; the final optimizer reports convergence.
That flag is not inferred from FSC. Pilot-selected alignment and interpolation
also affect cross-method/reference curves. Native half FSC is retained separately.

The primary EMDB entry identifies EMD-8434 as the Class A bL17-depleted large
ribosomal-subunit assembly intermediate, not the complete heterogeneous
population: https://www.ebi.ac.uk/emdb/EMD-8434 (checked 30 September 2026).
The author's EMPIAR-10076 tutorial describes assembly states and impurity
filtering: https://ez-lab.gitbook.io/cryodrgn/cryodrgn-empiar-10076-tutorial.
The present homogeneous reconstruction is a consensus approximation. Its weak
agreement may combine state mismatch, frame-search limits and reconstruction
error; the data here do not separate those explanations.

The generic alignment helper's nested `scope` string says "no reference
access" because it was originally used to align RELION to a pilot. In these
phase-A records its *source input is explicitly the deposited reference map*;
its target and selection criterion are only the old pilot. The parent record,
input hashes and protocol make that role explicit. The old inherited string
is not a claim that this registration did not read the reference. Historical
records and source bytes are retained unchanged.

Historical reference-based uncertainty/power calculations use the original
unregistered generators and are not silently corrected. Those are conditional
simulation/implementation outcomes. They cannot support biological localization
in the consensus frame without a separately reported registered-target study.
