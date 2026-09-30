# Response to the enclosing-domain audit

This authentic focused audit confirms the real-arithmetic identity, enclosures,
chain rule and residual-phase correction. It is not full ICML review 2. All
original records and the unmodified report are preserved.

- E1: the actual source record explicitly uses the joint scaled five-dimensional
  ball. The missing guard was real even though this case did not have a class
  mismatch. `validate_source_class` now rejects product or unknown pose sets,
  passes `joint_ball=True` explicitly, and records the class and flag. Separate
  product-ball tests exercise the sum-speed bound and show the flag changes it.
- E2: the source fit explicitly states radius 2 around a unit-norm pilot in its
  class string. The legacy loader requires that exact string, or explicit
  numeric fields, and rejects missing/unknown classes. New source audits record
  both numbers. The probe recomputes both original joint/triangle branches from
  saved scales, residual/cross norms and pilot pairing, checks equality, then
  recomputes both with the replacement remainder. It no longer subtracts an
  assumed additive term without that check. The v2 records include the replay.
- E3: targeted tests now cover a near-tight pure translation (bound/exact norm
  between 1 and 1.04), a closed-form antipodal cosine field that isolates the
  sum-frequency sign, a detector embedding scaled by 0.3, product directions,
  and angles at and above pi. A deterministic grid covers three rotation axes
  and eight translation directions. This is falsification testing, not a proof
  of a supremum computed by sampling.
- E4: an independent central finite-difference test checks first and second
  derivative fields. The product-grid test also constructs its own phase
  derivatives from the skew matrix instead of calling the project helper.
- E5: independent mpmath 60-decimal evaluations of third-derivative squared
  norms on particles 0, 512 and 1023 check both ball and cube kernels. All six
  discrepancies lie inside the heuristic pads; the maximum relative discrepancy
  is 4.501e-16. The original float inputs and domain radii are represented at high
  precision. This is a spot check, not validated rounding for every particle,
  embedding SVD or downstream operation. See the immutable precision record.
- E6: `cell_forward` delegates to `VoxelObservationOperator.forward`, whose
  FINUFFT call uses `isign=-1` in `uq_data.py:39`; its cell-center phase is also
  negative. Thus the positive adjoint and negative translation phase agree.
  Those sources were omitted from the focused packet, so the review correctly
  identified a verification limit rather than an established sign error.
- E7: the v2 records save domain, pose flag, B/P metadata and per-particle path
  speeds. The original ball-only record is not rewritten; the reproduction
  guide identifies its older schema. NumPy integer derivative orders now work.
  The proposal preserves its initial derivation and dates the subsequent joint
  speed refinement, rather than rewriting the chronology.

All 14 focused tests pass in 0.80 seconds; the complete suite passes 116 tests
in 6.57 seconds. Both guarded v2 probes complete and reproduce the original
bounds: ball-selected cubic 32.6233 with no data, cube cubic 21.1103 with width
0.75139 and zero reference feature power. No scientific-usefulness claim is
inferred from these fixes. The reviewer did not execute any tests; these are
project checks, not execution attributed to Fable.
