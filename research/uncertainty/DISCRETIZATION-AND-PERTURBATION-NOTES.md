# Discretization and perturbation uncertainty: targeted reading

Checked 30 September 2026. These are additions to the reading ledger, not three
new reproduced baselines. Exact downloaded-file hashes are in the background
manifest; third-party full texts remain local and are not republished.

## Population discretization

[Mordant et al., arXiv:2609.01688v1](https://arxiv.org/abs/2609.01688v1),
1 September 2026. Read the introduction, model, Theorem 1, Sections 5.2–5.4,
Theorem 3, and the simulation setup; not a complete proof audit.

The analysis fixes candidate conformations and their likelihoods. Its simplified
forward operator is identity with known isotropic Gaussian noise. Theorem 1 gives
a uniform high-noise KL lower bound from the first unmatched moment order.
Boundary weight estimates have a constrained Gaussian limit; ordinary interior
normal approximations do not cover every case. Hsp90 simulations use one shared
pose without CTF or translations. This is population reweighting, distinct from
our continuous spatial-density audit. It supports reporting candidate exclusion,
statistical fluctuation and incomplete optimization separately; it supplies no
experimental pose-confidence radius for our application.

## Dynamic NeRF uncertainty as a prior proposal

[Patel, *Conformally Robust Decision Making*](https://www.ambujtewari.com/theses/Yash_Patel_Thesis_2025.pdf),
University of Michigan dissertation, 2025. Read Chapter 5's introduction and
Section 5.2, printed pages 57 and 59–62 (PDF pages 69 and 71–74); the dissertation
as a whole was not reviewed.

Section 5.2 appears under **Future Directions**. It proposes extending the
Bayes' Rays spatial perturbation construction to a dynamic deformation field,
with cryo-EM as an eventual application. The underlying fitted field stays
fixed while a displacement grid receives a local Gaussian approximation. This
section supplies neither a completed cryo-EM experiment nor a demonstrated
distribution-free density-coverage theorem. Its existence precludes treating
the general dynamic-NeRF-to-cryo-EM uncertainty idea as unclaimed territory.
The earlier conformal decision-making chapters and the future proposal must
not be conflated as a proved cryo-EM guarantee.

## Bayes' Rays and the target of spatial uncertainty

[Goli et al., CVPR 2024](https://openaccess.thecvf.com/content/CVPR2024/html/Goli_Bayes_Rays_Uncertainty_Quantification_for_Neural_Radiance_Fields_CVPR_2024_paper.html).
Read Sections 3–4 and 6, plus the discretization and cleanup ablations in
Section 5.3; not a reproduction of its published experiments.

The method perturbs a frozen NeRF through an interpolated displacement field
and approximates its uncertainty using curvature. The practical approximation
discards off-diagonal curvature. The reported target is spatial epistemic
uncertainty, with depth-error ranking and artifact cleanup evaluations; the
paper explicitly excludes aleatoric inconsistency. A displacement variance is
not an additive cryo-EM density interval. A comparison would require a declared
Fourier imaging likelihood, noise scale, target propagation and treatment of
reference and pose errors. We therefore cite it as adjacent methodology, not
claim to have run a matched external reconstruction baseline.

## Remaining access limit

The CryoDiff primary abstract remains accessible through search, but a fresh
attempt to open its primary full-text page on 30 September returned HTTP 403.
Its Monte Carlo confidence metric is an abstract-level claim in this survey;
the exact calibration procedure and full experimental design remain unchecked.
The earlier failed PDF attempt and its status are retained in the manifest.
