# Response to the focused mixture audit

30 September 2026. The unmodified review and full provider stream are retained.
The invocation names and returns `claude-fable-5-1`; it read supplied materials
without execution. This is a mathematical development audit, not a second full
ICML review. The full round-1 rejection remains unchanged.

- **L1, common viewing law.** Accepted. The current simulated null assumes a
  shared G independent of known imaging operators. This cannot be silently
  applied to experimental defocus/orientation dependence or dependent particles.
  A future experimental null must allow separate viewing laws in independently
  declared operator groups and model within-exposure dependence. No experimental
  likelihood-test result is currently claimed. Grouping alone does not prove
  independence within groups or between particles; these remain open inputs.
- **L2, independent Fourier coordinates.** Implemented an explicit no-DC,
  no-duplicate, no-conjugate-duplicate guard in `uq_mixture_disks.py`, with
  falsification tests. The declared disk diagnostic applies it to all three
  existing simulation arrays. It cannot establish independence of whitening,
  noise, CTF estimation or experimental observations. Those remain assumptions.
- **L3, upper-bound names.** The historical `independent_image` value is exactly
  `normalization + sum_i max_c log U_ic`. The historical `row_max` value is
  instead the **scaled Lindsay bound using row maxima as anchors**:
  `sum_i a_i + n log(max_c sum_i exp(log U_ic-a_i)/n) + normalization`,
  with `a_i=max_c log U_ic`. Its second term is nonpositive, so it can be below
  `independent_image`. The values are consistent; the short labels were unclear.
  Preserve the old records and document these formulas. The review's coupling
  comparison should use independent-image minus fitted-envelope bound:
  289.525, 253.101 and 110.083 log units, rather than 101.655, 56.330 and 32.656.
  This remains small compared with spatial slack and does not change its main
  diagnosis.
- **L4, trivial upper cap.** Accepted as a free computational refinement.
  Existing anchor diagnostics already report the Gaussian peak explicitly.
  The original solver and currently running refinement are retained unchanged;
  a subsequent solver revision should include the cap in every reported upper.
  Their uncapped values are loose, not invalid, and are not silently overwritten.
- **L5, parent inheritance.** Already implemented in the independently published
  refinement; original results preserve their weaker implementation.
- **L6, final feasible refit.** Accepted. The currently running prescribed fit
  keeps its stopping behavior. A subsequent revision must refit at wall exit
  and report that extra work separately, without changing the saved old result.
- **L7, feasible support.** Accepted. Both evaluated-center and oracle-support
  values remain reported, each only a feasible lower bound. The oracle mixture
  is not the true Haar likelihood or a known optimum. No gap is called a known
  error of the upper bound.
- **L8, floating point.** Accepted and disclosed. There is no interval-arithmetic
  guarantee; priority clipping cannot supply one.
- **L9, nested stopping.** Agreed. Outer and inner gaps retain separate fields.
  Nonconverged inner fits still yield upper bounds in real arithmetic.

For Change 1, the new disk derivation preserves each complex frequency's
two-dimensional error block. An independent conic optimum check, exact example,
dual replay, and nonlinear checks all pass (three tests). An exact product-disk
minimum dominates the old whole-image-ball minimum, but finite dual iterates
need not; retaining the minimum of old and new likelihood uppers prevents a
regression. No empirical gain has been established at the time of this response.
The protocol specifies all 54 earlier local boxes and at most 95 selected
full-cover cells per geometry before running this diagnostic.

For Change 2, the same protocol computes fixed-map likelihood contrasts at 256
seeded Haar orientations, retaining all values. Such a finite diagnostic can
identify a bottleneck at specified resolutions or budgets; it cannot prove
that adaptive refinement fails at every possible budget or SNR. For a fixed
smooth finite Gaussian map, compact SO(3), fixed finite observations and positive
noise scale, arbitrarily fine covers remain theoretically available. No
computational impossibility theorem follows from a quantile comparison.

Remaining statistical development includes a normalized independent predictive
numerator, translation and noise nuisance treatment, heterogeneous maps and
experimental acquisition validation. Neither this audit nor the new numerical
tests resolve the full review's usefulness and novelty objections.


## Completed diagnostic follow-up

All 54 local boxes and 285 sampled coarse cells complete. Two-degree median
local gaps improve slightly on 10028/10049; every coarse-cell envelope remains
unchanged across all 36,480 image/cell coordinates. All three frequency guards
pass. The declared Haar contrast diagnostic and all source/array hashes are
retained in continuous-mixture-disks-v1. The larger ongoing refinement has
also completed with gaps 12,898--15,892. These results confirm that the tested
changes do not resolve the coarse global enclosure; no acceptance or useful
continuous-test claim follows.
