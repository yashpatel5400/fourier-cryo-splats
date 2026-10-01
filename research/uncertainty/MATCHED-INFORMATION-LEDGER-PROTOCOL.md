# Matched information ledger after review 4

1 October 2026 UTC. Diagnose existing frozen candidates and scores. Do not fit
new directions or vary the declared test, cap, amplitude range or particle
counts. This is a post hoc information comparison, not a new power experiment.

Stage A replays the first 8,192 held-out views/noises from each Fisher study
(stage-2 seed 261017 + dataset + 2,000,000, batches 512). Verify the replay
against every saved amplitude-one score in that prefix. Use the identical
candidate, region, frequencies and transfer profile for every comparison.
Retain deletions 0/.1/.25/.5/1 and the existing four score directions. For each
combined direction also inspect its power and bispectrum components without
refitting or re-normalizing them. Component diagnostics are not new tests.

For each score calculate exact conditional Gaussian means and variances using
the already verified finite Hermite expansion; average over the replayed
views. Report noise and view contributions to null variance and squared mean
separation `(E_alt T - E_null T)^2 / Var_null T`. These view integrals remain
Monte Carlo estimates. For known poses and unit real/imaginary noise variance,
twice the conditional Gaussian KL is `deletion^2 E ||region_mean(R)||^2`.
This identity uses the same map/region, unlike the earlier cross-map comparison.
It is also the squared separation of a centered known-pose linear score; it is
not the raw unconditionally varying log-likelihood ratio's variance formula.

Record the power-block weight sum. A constant real/imaginary noise-variance
error epsilon shifts its expectation by `epsilon sum(w_power)`; report the
absolute epsilon matching the 25%-deletion mean gap. This assesses a specified
white-noise variance error, not arbitrary colored-noise misspecification.
Also report the null *mean-score* envelope under one global amplitude and under
view-dependent amplitude in [.9,1.1], using the analytic quadratic/cubic mean.
These are different from the event-probability envelopes in the original test.
Do not subtract them from event-probability quantities or identify finite-L
bias from the two previous unmatched calibration realizations.

Stage B (separate runner, declared here before stage-A outcomes) compares
numerical Haar-marginal likelihood ratios on the same 8,192 paired null/altered
images, amplitude one and the four positive deletions. Use two disjoint banks
from the existing 65,536 training rotations, each with nested prefixes
4,096/8,192/16,384/32,768. These are quadrature convergence checks, not changes
to particle-count or statistical-test choices. Save per-image log ratios and
effective sample size / maximum quadrature weight. Compare bank/prefix
stability and signed mean/variance separation against the matching moment and
known-pose quantities. This approximate integration is not a certified
information ceiling. An unstable approximation cannot establish impossibility
or trigger a claimed tenfold information advantage.

The original simulator still has one CTF and known independent white noise.
Recorded nuisance distributions from the separate inventory are not magically
incorporated in this best-case comparison. No new experimental validity or
regional-specificity claim follows. Later conditional event-envelope diagnosis,
noise-spectrum work and a statistically calibrated latent-nuisance method
remain separate unfinished requirements. Preserve every result and failure.
