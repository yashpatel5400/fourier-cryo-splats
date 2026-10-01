# Frozen arithmetic gate for a population-uncertainty direction

1 October 2026 UTC. Written after the focused Fable consultation and metadata
field inspection, before computing new regional information, odds, roots or
harmonic moments. This is a bounded development screen, not preregistered
experimental validation. No new method is selected by writing this protocol.

## Fixed inputs and target

Use all 65,536 stage-0 Haar rotations and their saved physical means from
`candidate-fisher-score-v1` on each of 10028, 10049 and 10076. A is the already
fitted full normalized Gaussian map; B deletes its previously fixed 20-Angstrom
peak region completely. Do not relocate, resize, rescale or select a new region
after outcomes. Use the existing 220 frequencies, one CTF and independent
unit-variance real/imaginary white noise. Population p=.75, nominal n=10,000.
This is a two-template known-pose calculation, not a reconstruction experiment.

Replay all rotations from seed 261017+dataset, 128 batches of 512, consuming
the original real and imaginary noise draws between batches. Assert the
archived first rotation and independently evaluate direct cell sums at indices
0, 511, 512, 32767, 32768 and 65535 against both saved mean arrays. Freeze
source/protocol in Git before execution; record hashes and retain failures.

## Recorded law proxies fixed in advance

Use the same ZYZ equal-Haar-volume bins as the recorded metadata inventory.
The primary grid has 8 bins along each of alpha/(2pi), (cos beta+1)/2 and
gamma/(2pi), hence 512 cells. Grids 4 and 12 are sensitivity analyses. Within
each cell define the law to be Haar, with cell mass equal to its metadata
frequency. Evaluate that law with all saved training rotations in the cell.
No nearest-neighbour interpolation, smoothing, clipping or omitted empty cells.
An occupied metadata cell with no training rotation makes that grid unavailable.
This explicitly specified piecewise-Haar law is only a proxy for recorded
orientations, which themselves are estimates rather than latent truth.

Prelisted pairs on every applicable stack:

1. The published source half 0 versus source half 1, using all metadata.
2. Published filter retained versus its complement on 10028 and 10076.

Report both assignments of each pair to A/B. Do not select view bins, alpha
thresholds, acquisition clusters or new subgroups to maximize bias. Record
the full metadata law and Haar as additional harmonic summaries. Half labels
and filter membership are **not** conformational labels; the source `.cs`
class fields are all constant 0 with posterior 1. A measured biological
state-law pair is unavailable. These proxies can motivate a stress test, not
demonstrate state–view confounding in an actual population.

## Quantities and gate

Compute s=||g||^2, the full-map amplitude-tangent residual s_perp and the
secondary residual at m_B. Use real inner products. Save every per-rotation
value. Report quantiles, arithmetic and harmonic moments for Haar and every
law, and the exact known-amplitude nominal SE where its assumptions hold.
Report local projected nominal SE as a **surrogate**, not a confidence limit
or proof of identifiability with unknown amplitude.

For each pair and assignment, report the weak-signal odds approximation and
the exact known-pose Gaussian mixture population-score root derived in
`POPULATION-SCREEN-ALGEBRA.md`. Use 64-node Gauss-Hermite quadrature for root
finding on [1e-7,1-1e-7], root tolerance 1e-10. Independently check each root's
score with 128 nodes. If any normalized score check exceeds 1e-6, report
numerical uncertainty and do not declare a pass from that root. Do not change
the quadrature after inspecting an outcome. Report arithmetic weak-signal
results as approximations even when they agree numerically.

The primary arithmetic gate passes a stack only when **the same prelisted
pair and assignment on the 8-bin grid** has:

- exact population bias |t*-.75| >= .02; and
- the mixture-weighted harmonic surrogate
  .75 E_nu_A[1/s_perp]+.25 E_nu_B[1/s_perp] <= 8.8125.

At least two of three stacks must pass. Sensitivity grids cannot rescue a
primary failure. No screening significance test, biological inference or
novelty conclusion follows from crossing these descriptive thresholds. Report
within-cell Monte Carlo SEs of linear summaries conditional on the fixed
metadata histograms; they do not include metadata error, pose error or model
misspecification. No resampling or extra rotations are used to force a pass.

## Decision and stopping rule

If fewer than two stacks pass, do not build the proposed response-calibrated
population method on these compact regions. Preserve the negative result and
its scope. A failure does not prove pooled likelihood adequate for other
states or global impossibility of occupancy inference. A pass authorizes
only a separately frozen, nuisance-aware feasibility design; it does not
establish that biological state laws differ or that G2/G3 have been run.
There is no repair, new region, modified threshold or favorable subset inside
this screen. No external GPU or monetary spend is needed.
