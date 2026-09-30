# Reviewer-motivated disk-envelope diagnostic

30 September 2026. Declared after focused audit mixture-audit-01 and its
three independent implementation tests, before the following calculations.
Reuse all three completed continuous-mixture-v1 simulations. No new particle,
signal, or noise sample and no validation rejection rate is computed.

1. Replay all 54 oracle-centered local boxes from the curvature diagnostic,
   with their existing feasible local-search values, using the disk envelope.
   Report every old/new upper, feasible value, and within-disk primal/dual gap.
   These oracle-centered boxes are local diagnostics only.
2. For each original final full cover, compute the column dual scores at its
   saved best anchor. Sort stably by score. Select 80 evenly spaced sorted
   ranks including both endpoints, together with the 16 largest scores.
   This selects at most 95 unique cells. Re-evaluate all 128 images in every
   selected cell. Record per-image old, recomputed, and disk bounds and all
   disk gaps. Replacing only those cells by the minimum with their new
   envelopes still bounds the full cover; report its unchanged-anchor upper
   and the original, with no assertion that either is tight.
3. Draw 256 independent Haar rotations with seed 950001+dataset for numerical
   evaluation of the fixed map, not for generating new observations. Record
   every log kernel and compare each image's best previously evaluated center
   with the 50th, 90th, and 99th percentiles over these rotations. Compare with
   the existing cell-envelope slack. These finite contrast diagnostics do not
   prove failure at every computational budget or statistical impossibility.

Keep known noise, known CTFs, zero translation, common viewing-law and
homogeneous-map assumptions explicit. Validate each half-plane has no DC,
duplicates or conjugate duplicates. Save arrays, source/input hashes, all
prescribed records, and failures. The total calculation budget is 1,800 seconds,
checked between completed cells/evaluations; do not extend it after outcomes.
The old source, ongoing refinement and original results are not modified.
