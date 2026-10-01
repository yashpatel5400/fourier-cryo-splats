# Response to full independent review 4

1 October 2026 UTC. The [unaltered review](round-04/review.md) recommends
**reject**, confidence 4/5, and says the work is **not a strong ICML contender**.
The provider reports the exact requested `claude-fable-5-1` model. All 33 pages
were visible, all 2,684 evidence copies remained unchanged, and the call returned
success. The reviewer read source and records, but ran no code, opened no NPZ
arrays and consulted no third-party papers. All four full reviews remain public.

## Scientific decision

Freeze the moment-test branch. No new score, shrinkage, calibration bound,
replica allocation, viewing cap or sample-size sweep is justified by the
current results. The next work is a diagnosis of existing constructions and
measurement assumptions, not another attempt to improve their headline power.
The overall research task remains unfinished; ending an unproductive branch
does not mean completing it.

The review identifies a substantive problem: for important design alternatives,
the estimated alternative event probability is below the conservative null
bound. Adding particles then reduces rejection probability. A covariance-aware
score or a smaller Monte Carlo radius cannot by itself remove the physical
nuisance ambiguity. We accept this as the priority over further simulations.

## Finite next work

1. **R17/R22: matched information and nuisance ledger.** Use the actual three
   candidate maps and identical nominated regions for known-pose Gaussian
   discrimination, numerical Haar-marginal likelihood comparisons, and the
   saved moment scores. Report the full decomposition of amplitude/view/noise
   effects, finite-calibration slack, power-block weight sums, and the
   bispectrum block's retained separation. Do not compare different maps or
   masks as if they supplied an information ratio. No newly optimized score.
2. **R18/R19: inspect all three experimental stacks.** Quantify the recorded
   viewing directions with bin-resolution and split-half checks, defocus
   variation, available amplitude fields, and existing noise-spectrum evidence.
   Consensus poses are estimates: their histogram is not a certified bound on
   the true continuous viewing density. Constant metadata fields do not certify
   constant physical amplitudes. Report what cannot be measured from the
   archived inputs instead of converting placeholders into uncertainty bounds.
3. **Decide the estimand once these diagnostics finish.** The present paper
   cannot combine density intervals and a whole-candidate point-null test as
   one method (R16). A latent-pose likelihood direction must first have a
   precisely stated nuisance model and a written prediction of useful power.
   Regional occupancy requires nuisance hypotheses that actually isolate the
   region; otherwise use an omnibus compatibility claim (R20). Population UQ
   is a distinct possible target, with direct cryo-BIFE/reweighting prior art,
   not a presumed new contribution.
4. **R21 only after a useful candidate survives.** Freeze independent
   calibration-plus-test repetitions, boundary and adverse viewing controls,
   the primary size/power/specificity criteria, and experimental held-out
   controls before outcomes. Retain every outcome and do not repeatedly rerun
   calibration until a favorable size estimate appears. Further repetitions of
   an already insensitive test are not the immediate priority.
5. **R11/R14 if alignment inference is retained.** The .72 finding is specific
   to a degrading local aligner. It needs a working-aligner gate and realistic
   pose radii before any broader claim. No more repeats of the same failed
   refinement. The earlier exact envelopes quantify a particular realized
   error, not a general information lower bound.
6. **R16/M31/M37: one paper after the decision.** Move development history and
   review commentary to the repository. Rewrite the abstract, method and
   experimental section around the surviving estimand. Preserve earlier PDFs
   and all failed branches as immutable artifacts. Another full acceptance
   review is warranted only after a substantive scientific result.

Adopt the review's branch stopping rule: every proposed change must target the
largest measured obstacle and state a falsifiable prediction before execution;
two consecutive missed predictions end that branch. A failed preregistered
gate after one substantive revision also ends it. This is a research decision,
not permission to select a favorable subset or claim the overall goal achieved.

## Qualifications and factual checks

- Numerical mixture likelihoods are not certified information ceilings. Their
  integration and discrimination estimates need convergence/error diagnostics.
  Failure of an approximate likelihood implementation is not proof that every
  method fails. A claim of impossibility needs a valid lower-bound argument.
- A finite-simulation estimate of a population envelope is not its exact value.
  Retain uncertainty when describing the inferred asymptotic detection floor.
  The conditional failure to cross the implemented bound is directly checkable.
- **M32 is a transcription error.** The allocation CSV must supply the corrected
  10076 number in the next manuscript revision. The v0.7.6 PDF is preserved.
- **M36 is factually mistaken.** The round-3 response is present at
  `round-04/evidence/research/uncertainty/reviews/response-to-round-03-development.md`
  and indexed in the review packet; a Read call in the raw event stream accesses
  it. This correction does not affect the scientific rejection.
- The 23-condition author-code population replay was frozen and started during
  this review, before receiving the verdict. It is a separate retrospective
  comparator check on published likelihoods, not another candidate-test variant
  or a result on CAHRA. Its primary-reading notes were outside this review's
  immutable evidence.

R1/R2/R7/R9/R11/R14/R15 and R16–R22 remain scientifically open. Correct algebra,
reproducible code and a larger experiment inventory do not resolve them.
