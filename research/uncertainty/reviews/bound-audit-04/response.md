# Response to the pilot-pairing audit

The authentic `claude-fable-5-1` report confirms the known-pilot inequality,
constant-cell moment derivation and pose contraction. It closes B9 and B10
and re-confirms the positive-scale design argument. This is a focused
mathematical audit, not full review round 2. The original full-paper rejection
and unresolved experimental-calibration/usefulness objections remain in force.

**B13, stale inventory:** regenerated the summary after the pilot batch. It
now explicitly identifies `pilot_pairing` for each variant and provides a
separate pilot-refinement table, without changing the scope of the original
joint-only comparison. At this response snapshot, there are 35 completed
audit variants, including eleven pilot-refined cases. The source-fit inventory
also retains unattempted and unfinished cases. These are inventory counts,
not calibration evidence.

**B14, pilot numerical evaluation:** added direct exponential sums over actual
physical cell centers, without a NUFFT or centered-grid phase. For each of the
first ten archived pilot audits, all twenty lifted columns for particles 0,
64 and 127 agree with the fast evaluation. The maximum relative discrepancy
is 3.4818e-13; maximum absolute discrepancy is 7.5157e-16. See
`results/uncertainty/development/audit-regressions/pilot-pairing-direct-checks.json`.
Original audit records and their source snapshots remain unchanged. Future
pilot audits, including the larger-band run, include the check in their own
record and preserve a failed outcome before raising. The direct sums share
scalar within-cell formulas, separately tested against composite quadrature.
This check is expressly selected-column numerical evidence, not an all-column
rounding-error bound. No interval-arithmetic guarantee is claimed.

**B15, source residual assumption:** the higher-band post-audit now rejects
sources containing pose-polynomial or cubic bias terms, and requires the
source to identify its fixed-pose/CTF scope. Thus `fit['bias']/B` is used only
for the intended fixed-pose residual. The larger-band audit started after
this guard was installed; it is still running at this snapshot.

The numerical changes and bounded-optimizer option pass the full suite of
93 tests in 7.15 seconds. The new pilot formulas and table have been compiled
in the 36-page manuscript; pages 20 and 21 were rendered and visually checked.
The bounds remain conditional on supplied density, pose and noise assumptions.
