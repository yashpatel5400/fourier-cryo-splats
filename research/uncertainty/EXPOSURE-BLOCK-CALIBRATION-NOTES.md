# Possible exposure-block calibration: theory and metadata only

30 September 2026. This is a possible subsequent research direction, **not an
implemented method, experimental protocol or favorable result**. Finish the
current joint/RELION/application batch and full review before deciding whether
this addresses the substantive objections. No additional particles were
selected or downloaded for this note.

The current grouped bound uses Cauchy--Schwarz within each inference exposure:
Var(sum_i w_i^T epsilon_i) <= m_g sum_i Var(w_i^T epsilon_i).
It allows arbitrary within-exposure dependence using calibration that contains
only one particle per exposure. On the selected cubic cohort, m_g reaches 12.
This is a genuine information limitation of that calibration design, not proof
that the actual particle noises are perfectly correlated.

Suppose a future acquisition model instead justifies one common joint Gaussian
covariance Sigma_M for M ordered raw-coordinate particle noises per exposure.
For smaller inference groups require the appropriate principal marginal of
this same covariance. Pad each group's weight vector to length M*d with zeros,
and collect the vectors as columns of K. Independent inference exposures then
give exactly

    Var(sum_g w_g^T epsilon_g) = tr(K^T Sigma_M K).

If independent calibration exposure vectors Y_h have this same covariance,
the projected rows K^T Y_h have covariance K^T Sigma_M K. Their means may be
arbitrary. The existing noncentral Gaussian trace bound therefore applies
directly, without the m_g inflation. Fixed orthogonal contrasts across
calibration exposures are also allowed under those assumptions. This is a
direct reuse of the existing Gaussian identity, not a new statistical theorem.
No improvement is guaranteed: signal variation still inflates observed energy.

The scientific cost is a **stronger and different covariance assumption**.
Arbitrary common particle-marginal covariance does not imply common ordered
block covariance. Particle spacing, crop overlap, preprocessing, selection and
exposure-dependent backgrounds can violate it. An exchangeable random-effect
model would be an additional assumption, not something supplied by the formula.
Nor does this remove supplied-pose dependence, picking bias, heterogeneity,
or the need to justify density/pose bounds. An empirical check and a defensible
block-sampling rule are required before choosing such a model.

A metadata-only inventory (no pixels accessed) finds that the already reserved
128 calibration exposures each contain at least 12 eligible particles on
10028 and 10049; 122 of 128 do on 10076. Counts use the published source filters
and known truncated-source exclusion. The full old inference-half-0 pools
contain respectively 2,030/1,677/2,283 particles but only 58/30/93 exposures.
Thus larger particle counts alone do not create independent exposures. Exact
input hashes and counts are in
[the metadata inventory](../../provenance/uncertainty/exposure-block-metadata-feasibility.json).
These are availability counts, not calibration outcomes or a new independent
confirmation cohort. The original representative pixels have already been seen.
