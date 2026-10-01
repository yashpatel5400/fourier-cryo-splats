# Equal-noise-budget allocation diagnostic, version 1

Declared after inspecting the fresh Fisher-score-v1 results. That study finds
limited sensitivity to smaller errors even with a covariance-aware contrast.
This follow-up tests an allocation choice, not a newly claimed estimator.

Freeze all four trained Fisher-study scores and thresholds on each of the
three candidate maps. Retain both matched/Fisher and power/combined families.
Compare the already frozen M=32,768 views, two L=32 noise groups calibration
with a new M=8,192 views, two L=128 groups calibration. Both use 2,097,152
conditional noise draws. They have different numbers of Fourier evaluations
and different costs; do not call this equal wall time or an optimized allocation.
Retain all seven probability-bound procedures, all five density caps and all
three sample sizes separately. Never select the smaller bound as an
unadjusted combined procedure.

Increasing L can reduce upward bias from maximizing noisy amplitude-cell
averages, but reducing M increases outer Monte Carlo uncertainty. Neither
allocation is uniformly tighter by theory. The convexity identity
max_j (x_j+y_j)/2 <= (max_j x_j+max_j y_j)/2 explains why the nested mean
envelope decreases in expectation under doubled replication at fixed view;
it does not promise a smaller finite-calibration upper confidence bound.
This is classical nested-simulation bias, not a new theorem.

The new calibration retains the same 16 amplitude cells [.9,1.1], b=.5,
delta=.001. Fresh held-out testing uses 131,072 independent Haar rotations,
all five deletions 0/.1/.25/.5/1 and amplitudes .9/1/1.1. Both allocations
use these same new scores, hence neither reuses viewed test outcomes.
Actual conditional testing partitions the new draws into 128 groups of 1,024.
Report 37,800 scalar projections and 12,600 repeated-group cells across the
three stacks, and every correct-null control. Do not pool dependent cells.

Seeds 261018+int(dataset)+stage*1,000,000 for stages 1 calibration and
2 testing; batches 8 and 512. Preserve source/input hashes, saved Gaussian
maps, region, directions and thresholds, all RNG states, grouped/individual
calibration envelopes, and all held-out scores. Independent direct physical
sums check the first view of both stages. Calibration repetitions and actual
experimental-noise validation remain outstanding.

Relevant primary literature on this numerical issue includes Giles and Goda,
Decision-making under uncertainty: using MLMC for efficient estimation of
EVPPI (arXiv:1708.05531), and Giles and Haji-Ali, Multilevel nested simulation
for efficient risk estimation (arXiv:1802.05016). Their abstracts and the
author's nested-expectation overview were consulted at protocol declaration;
no MLMC implementation or complexity claim is made here.
