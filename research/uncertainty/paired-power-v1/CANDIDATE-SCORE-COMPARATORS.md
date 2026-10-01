# Additional candidate-score comparisons declared before calibration completion

The three candidate-score processes are running from 932ea23. Only training
means/directions and progress have been inspected; calibration is not yet
complete and no held-out rejection result has been computed. To carry the
focused Fable audit's comparator requirement into this study, apply the same
four independently valid classical risk procedures (`cvar_dkw`, `cvar_split`,
`mean_only`, `unpaired_variance`) to each candidate's saved grouped arrays,
using center .5, the same five kappa values and delta=.001. The implemented
formulas and separate confidence budgets are those frozen at 82f8403. No
method is selected by taking an uncorrected minimum after inspecting outcomes.

Preserve the original three-method runner unchanged. Its 540 rows remain the
primary archived output. Add 720 comparator rows (12 scores, five ratios, four
procedures, three sample sizes), each with all 15 deletion/amplitude outcomes.
Together these give 1,260 rows and 18,900 scalar rejection projections, with
both raw and scale-orthogonal scores on every stack. Each method/score is a
separate test, not an aggregate multiple-testing decision. Also retain the
actual repeated-group outcomes for the four added procedures, holding the
same calibration fixed.

**Interpretation boundary:** a rejection is a test against the complete fixed
candidate under the supplied image simulator. A statistic designed for a
regional deletion does not prove that the region is uniquely responsible for
a mismatch; other structural or imaging errors can affect the same moments.
Neither a non-rejection nor a directionally designed test certifies the density
or occupancy of that region. The fitted Gaussian maps are experimental inputs,
but all new image observations are simulated and their noise/CTF normalization
is inherited from a specified test fixture. This is not application to held-out
experimental particles or an experimental noise-calibration claim.
