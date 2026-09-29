# Experiment protocol draft — not yet frozen

The current results are development results. Method, dimensions and tuning may
change. Freeze a versioned confirmatory protocol before accessing its outcomes.
Every final table must point to exact config, code revision and data checksums.

## Estimands and labels

- Fixed linear density averages/contrasts, selected from an independent pilot;
  coordinate frame and smoothing scale specified in physical units.
- Image prediction intervals and likelihood, reported separately from density
  coverage. The real unknown density cannot supply empirical coverage labels.
- Homogeneous density uncertainty is distinct from physical conformational
  variation. Treat heterogeneous observations as a deliberate assumption stress
  unless a heterogeneous target is explicitly defined and implemented.
- For map-based simulations, a deposited map is the exact simulation generator,
  not the experimental ground truth. Gaussian-model, voxel-map and atomic/other
  independent generators must be clearly distinguished.

## Controlled coverage experiments

1. Exact Gaussian linear model, fixed known poses: sample several independent
   dictionaries/geometries and signal draws, with both uniform and anisotropic
   views. Include well-determined and weakly identified targets.
2. Both Bayesian prior-predictive and frequentist fixed-signal coverage. Preserve
   the correct prior for the former. Use random interior/boundary signals and
   each estimator's own adversarial directions for the latter. Do not generalize
   a deliberately adversarial result into a typical-case claim.
3. Physical nonlinear rotations and shifts: zero, random centered, coherent and
   adversarial errors. Evaluate correct bounds, overestimated bounds and deliberate
   violations. Compare retaining/dropping the density-pose interaction and the
   second-order remainder.
4. CTF/scale errors, correlated noise, heavy tails, heterogeneity and representation
   mismatch are separate stress axes. Current proofs cover fixed CTFs and known
   Gaussian noise. Either extend the model or label these as outside-theorem
   tests. Never claim a theorem handles an unimplemented nuisance parameter.
5. Sweep density-energy and pose sensitivity budgets. State whether budgets were
   chosen before outcomes, estimated using a pilot, supplied by an external
   assumption, or derived from simulation truth. Oracle budgets are diagnostic
   controls and cannot support practical calibration claims.
6. Report actual coverage with binomial intervals for independent replications,
   widths, normalized widths, sign-certification power, prediction error,
   reconstruction/FSC, runtime and memory. Correlated voxels are not independent
   Monte Carlo replications. Where exact Gaussian coverage is available, label
   it as analytic coverage and also cross-check simulation numerics.
7. Correct multiplicity for simultaneous maps or selection over sensitivity
   budgets. Data-dependent target selection on the inference images requires
   a valid simultaneous procedure or another independent split.

## Baselines

- Full Gaussian posterior in the same representation, with correctly specified
  prior predictive controls and several prior strengths.
- Diagonal variational covariance with the same posterior mean, documenting the
  relation to Ullrich et al. rather than claiming to have run their full code.
- Joint Gaussian/Schur-complement pose-volume covariance, including a correctly
  specified stochastic-nuisance control; cite joint Hessian prior work.
- Unregularized sampling covariance on identifiable subspaces; show extremely
  wide intervals rather than quietly dropping unfavorable targets.
- Ridge sampling covariance; particle and source-group bootstrap; where practical,
  full refitting/alignment versus conditional fixed-pose resampling.
- Classical bias-aware fixed-pose intervals, nuisance projection, coarse
  remainder certificate, structured interaction certificate and no-data bound.
- Representative actual reconstruction systems, including neural cryoDRGN where
  feasible. Clearly distinguish a trained neural model from its backprojector.
- Map confidence/enhancement tools only on matched tasks/targets. A LocScale
  baseline-map score is not a particle-derived density confidence interval.
- Conformal or empirical calibration only with disjoint calibration structures
  and honest exchangeability/domain-shift tests. Do not calibrate on the final
  truth or report noisy-image coverage as structural coverage.

## Data and splitting

The development pool contains 8,192 extracted particles from each of EMPIAR-10028,
10049 and 10076. All have now been mapped to repeated source exposure identifiers.
Pilot, tuning, two inference halves and test labels are assigned by a stable hash
of the group. These labels apply to development experiments; earlier work and
the representation audit have already inspected outcomes from this pool.

For final confirmation, use fresh simulated draws and select additional original
source groups absent from this entire development pool if real-data predictive
confirmation is needed. Record the original dataset's published particle filter
and any missing/corrupt source exclusions. Fixed published consensus poses and
CTFs remain a conditional-input limitation even after source-group splitting.
An independently estimated-pose experiment must actually rerun the alignment
pipeline on separate data; a different random split does not establish it.

## Compute and stopping

Use the Mac for development and measure actual scaling before requesting a GPU.
The provisional first cloud batch has a $100 ceiling as detailed in COMPUTE.md.
Stopping a compute batch is not scientific completion. A favorable result on a
small synthetic case is not grounds to stop baseline or failure-control work.

Before external review, the candidate must include complete results, proofs,
limitations, reproducibility instructions and an honest novelty comparison.
The requested independent reviewer is the actual Claude Fable 5.1; its availability
has been checked, but no scientific review has yet occurred.
