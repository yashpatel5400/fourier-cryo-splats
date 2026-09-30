# Prospective utility check at pilot-selected density regions

Status: proposal only; target coordinates are not yet selected and no reference
outcomes have been computed for this design. This does not alter any prior
experiment, class or result. The existing center and symmetric z-contrast targets
are coordinate-defined mathematical functionals, not regions established to
contain a biologically meaningful positive feature. Low sign power can therefore
reflect target choice as well as loose or intrinsically wide intervals. That
possibility motivates a new test, not a retroactive explanation asserted as fact.

Select targets from the independent, already fitted Gaussian pilot on each of
10028, 10049 and 10076. Use its same unit-L2 24-cell representation and exact
Gaussian cell integrals to find the strongest smoothed feature within 0.30 field
of the center. Use one physical Gaussian standard deviation of 20 Å for the
initial matched-scale probe. Select the three largest positive pilot features
with pairwise center separation at least 40 Å; if fewer exist, record that fact
rather than reduce the separation after viewing outcomes. Choose from a fixed
24-cell center lattice and use lexicographic tie-breaking. Freeze coordinates,
pilot/checkpoint hashes, selection source and parameters before evaluating any
deposited map. These are pilot-selected mathematical features, not annotated
atomic motifs or prevalidated biological conclusions.

Proposed first stage: same 128-particle/radius-12 development geometry across
all three stacks, known simulated noise inherited from the original radius-five
setup, B=2/P=1 continuous unit-cube class, full Gaussian cell target, order-80
quadrature, and a fixed 0/1/2-degree joint pose sensitivity grid with 0.5 Å shift.
Use existing fixed-pose continuous fitting followed by the sharp, cross-term,
known-pilot and enclosing-cube audit without silently narrowing the density class.
Compare every feature with the old central feature at the SAME physical width,
sampling geometry and noise. All selected features and comparison outcomes must
be retained. Since all are chosen from the pilot, not inference noise or a
reference, a conditional target statement remains meaningful under the original
pilot-independence assumption. If presenting simultaneous conclusions across
features, allocate alpha across them explicitly; do not relabel pointwise
intervals as a confidence map.

A radius-12 band has different physical endpoints on the three stacks. Record
them and the Gaussian target's missing-band sensitivity. A 20 Å standard deviation
is not a 20 Å FSC resolution claim. A later 10 Å study should first ensure the
physical band is comparable (roughly radius 24/12/21 for these fields) and account
for its higher compute cost. No favorable density-coverage claim can follow from
reference-map checks alone, and this design does not calibrate experimental pose
radii, the density ball or the noise law.

Before launch: implement and numerically check the pilot selection rule; freeze
all nine coordinates and the matched-center controls; implement a generic-target
runner with source-class guards; bound runtime for a single development case.
Only then execute the grid. Escalating particle count or GPU use should depend
on the resulting bottleneck, not on an assumption that more data fixes a biased
pose or representation model.
