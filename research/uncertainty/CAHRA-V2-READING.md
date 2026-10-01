# CAHRA version 2: a distinct uncertainty target

1 October 2026 UTC, during full review 4. This new reading is not in that
review's frozen evidence. CAHRA was already in the 138-candidate ledger;
this deepens that entry rather than increasing the candidate count.

The [primary manuscript](https://doi.org/10.64898/2026.09.15.751515), version 2
dated 24 September 2026, was retrieved through the publisher's JATS XML.
Sections 1–4 and figure captions were read. Figure images, supplementary
methods, data and author code have not been inspected. The source record
distinguishes successful API/XML access from failed browser/PDF/article access.

The three challenges address different estimands. Challenge 1 mixes separately
acquired experimental species, with source labels and radiation-damage variants;
these are not exact density truth. Challenge 2 uses molecular dynamics and
micrograph simulation, then consensus refinement perturbs the supplied poses.
Challenge 3 uses two known CALHM2 structures, a fixed population proportion,
strong preferred views, and variable Gaussian noise. Its 30,000 particles were
simulated with CryoJAX Gaussian-mixture projection. State and viewing direction
are sampled independently, with a 3% uniform-view component. The requested
outputs include population proportions and particle state/pose assignments.
The manuscript describes dataset construction, not full competition outcomes.

## Consequences for this project (our interpretation)

Population uncertainty for two nominated structures is a concrete possible
research question. It differs from both a confidence interval for homogeneous
local density and our current whole-candidate deletion test. A useful method
would have to report how population uncertainty changes with unknown viewing
distribution, amplitude/noise variation, imperfect templates, and violations of
state–view independence. The published proportions are now known, so any new
evaluation is post-release research, not a blind challenge submission. No
CAHRA result, dataset download or method implementation is claimed here.

The previous density-ratio cap is especially questionable for this setting:
a small uniform component supplies a lower density floor but does not bound
the peak of the remaining preferred-view distribution. A known 3% uniform
component therefore does not justify our earlier cap of 1.1.

Challenge 1's first/last 30-frame constructions are not automatically
independent exposure halves. The methods describe dropping the first or last
10 frames. Overlapping retained frames would share shot noise, and joint motion
correction introduces further dependence. The exact overlap and preprocessing
must be checked before treating such images as independent measurements.

The supplementary methods and released file inventory are still needed to
define a reproducible CAHRA experiment. This reading supplies a candidate
target and falsifiable assumptions, not a novelty argument or an acceptance
claim. Any pivot must first address the fourth review's scientific objections.
