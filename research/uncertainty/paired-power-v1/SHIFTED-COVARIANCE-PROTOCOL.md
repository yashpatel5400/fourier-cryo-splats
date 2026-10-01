# Common-translation gate after the full-frequency probe

This post-outcome development step responds to P3/P4 of the focused audit.
The uncompressed first probe finds positive finite-catalog growth on all
three stacks, but directions violate 26.3–44.0% of 10,000 new orientations.
Repair leaves positive achieved growth only on 10049. That neither defeats
the full family on the other two nor establishes continuous validity on 10049.

Retain all 14,160 original-plus-fresh null orientations for a new cone fit.
For each stack use shift SDs 0,.5,1,2 pixels of the downsampled 64-pixel box:
twelve cells. Use the same shared standard-normal shift draws across scales,
seed 261006+dataset ID. These are simulated common translations of both
exposures, not relative frame motion. The null permits arbitrary continuous
translations (equivalently a periodic 64-pixel torus at these integer Fourier
frequencies); the training catalog is a relaxation, not a bounded-shift proof.

For the alternative use the exact second moment averaged over the original
4,160 true-map orientations and an independent isotropic Gaussian shift of
the specified SD. With C_jk=E[F_j conj(F_k)] and P_jk=E[F_j F_k], multiply
C by exp(-2 pi^2 sigma^2 ||q_j-q_k||^2/64^2) and P by the analogous sum-frequency
factor, then realify both moments. This is exact for the stated shift law.
It avoids substituting noisy Monte Carlo moments for the known alternative.
Validate the implementation against independent Gauss-Hermite quadrature.

Fit a full 440-coordinate cone using the same kernel NNLS settings, with
120 seconds and 3,000 iterations per cell. The previous spectral upper is
valid for any nonnegative mixture even if the solve is incomplete. Retain
all statuses and gaps. Include a direct true-map control: its exact
shift-averaged moment is itself a continuous null mixture, giving zero
residual without relying on a sampled shifted catalog.

Generate a third set of 10,000 Haar rotations per stack, seed 261007+dataset,
and independent common shifts using seed 261008+dataset. Reuse these draws
across SDs. Save all Fourier values and rotations. Test each direction on
this set, retain the fraction/max violation, repair against the union and
recompute its achieved growth. A positive untested continuous-orbit supremum
still invalidates an amplitude-free test. The sampled union cannot replace
a continuous certificate.

No new noise observations, experimental coverage or rejection-power claim
occurs. The output is an information and discretization diagnostic. Both
phase/gain agreement between exposures and a defensible Loewner noise bound
remain required physical inputs. Do not launch a large confirmation study
merely because a finite-catalog mean score is positive.
