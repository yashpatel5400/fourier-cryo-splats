# Three-stack RELION reconstruction comparison

Development protocol, 30 September 2026; declare before fitting. This extends
the successful native runtime probe to an actual unknown-pose baseline. It
does not turn reconstruction agreement or estimated angular accuracy into
density coverage. Existing data and earlier reference results have been seen.

## Inputs and fixed choices

Run stock RELION 5.0.1 on EMPIAR-10028, 10049 and 10076, in that order. Use all
old `pilot` particles for a de novo initializer, then all old `inference_half0`
and `inference_half1` particles for refinement. Preserve their whole-exposure
half labels in `rlnRandomSubset` and check that RELION retains them. These are
the same refinement-image pools as the fixed-pose neural comparison. They are
64-pixel extracted images, not a new movie-processing or full-stack benchmark.
Use the original CTF metadata but supply no consensus orientations, translations,
reference density, learned regularizer, or deposited atomic model to fitting.

Center and scale each raw particle using background pixels outside 0.43 of the
field, and use minus its recorded cryoDRGN data-sign multiplier. Write one
optics group after verifying constant voltage, Cs and amplitude contrast.
Hash all input files and the fixed preparation/runner source. The target mask
diameter is 0.8 of the field on every stack. All fits use C1 and one class;
10076 is heterogeneous, so its result is a consensus approximation.

For each dataset, seed = 670100 + accession number. Fit 100 VDAM mini-batches,
initial/final batch sizes 128/min(1024, available pilot particles), two CPU
threads, initial resolution 60 A, final limit max(15 A, 3 pixels). Other search
and masking choices follow the runtime probe. Keep a model every ten batches.
The initializer has a 30-minute wall limit. If it fails or times out, retain
that outcome and skip that stack's dependent refinement; continue other stacks.

Refinement uses the single predetermined initializer, low-passed to 50 A,
with RELION auto-refine and separate halves. Use three MPI ranks, two threads
per worker, oversampling one, initial HEALPix order two, local searches from
order four, shift range five pixels/step two, normalization and scale correction,
and padding two. Join halves only below 1/40 A as in the stock job. Cap at
50 iterations and two wall-clock hours per stack. These bounds are workload
limits, not convergence criteria. No retry or parameter selection based on the
deposited reference is permitted within this protocol.

## Evaluation and reporting

Keep the return code, wall time, all logs/checkpoints and RELION convergence
fields. A timeout, iteration ceiling, nonzero return, or missing convergence
flag must remain explicitly unconverged. Report all three outcomes. Parse the
final/last model STAR records, including angular/shift accuracy; those estimates
are local model diagnostics, not calibrated uniform pose radii.

Compute unmasked half FSC on the native half maps, displaying the shared
low-frequency region and the sampling limit. If no threshold crossing occurs,
report it as censored. Do not relabel half agreement as reconstruction accuracy.
For cross-method/reference comparisons only, resolve the global frame using
the old independent Gaussian pilot, low-passed to 40 A. Align the summed RELION
map once and apply the same transform to both halves. Enumerate both hands and
the 24 proper signed permutations of principal axes, followed by local rigid
optimization; retain every start and choose only by pilot correlation. The
deposited map does not select a hand, transform, initializer or checkpoint.
Report alignment failures and the limits of this local search. Reference FSC
is approximate-reference agreement, with interpolation and shared-preprocessing
caveats. Compare against the existing Gaussian and neural maps on the same
inference pools. A successful initializer alone is not a converged baseline.

This protocol supplies an additional conventional reconstruction comparison,
not a full probabilistic RELION posterior or an experimental validation of our
uncertainty assumptions. Multiple independent starts and higher-resolution/full-
stack refinements remain outside this bounded development batch.

Primary command/convergence documentation:
https://relion.readthedocs.io/en/release-5.0/SPA_tutorial/InitialModel.html
https://relion.readthedocs.io/en/release-5.0/SPA_tutorial/Refine3D.html
https://github.com/3dem/relion/blob/5.0.1/src/pipeline_jobs.cpp
