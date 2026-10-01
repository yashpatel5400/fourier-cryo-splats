# Reproducing the post-review-4 diagnostic checkpoint

The v0.7.7-dev release contains four archives: common evidence and one numerical
bundle per stack. Each adjacent manifest lists every member, source commit,
byte length and SHA256. Archives are independently streamed and checked after
creation. The v0.7.6-dev per-stack score/calibration archives are dependencies;
their hashes are recorded in every new manifest. Extract archives at a clean
repository root and follow the environment setup in the README. For acquisition
metadata, use the pinned cryoDRGN EMPIAR repository and existing download
manifests; no third-party PDFs or external model reasoning are redistributed.

The paper PDF remains the rejected v0.7.6-dev snapshot. This release is an
experimental evidence increment, not a new manuscript or acceptance claim.

## Read-only reconstruction of the saved evidence

The stage-A replay saves the actual 8,192 candidate/region means and Gaussian
noise arrays per stack. Stage B saves all two-bank/four-level/four-deletion
likelihood ratios and diagnostics. The adaptive gate saves all 384 observed
images, proposal modes, every integration quaternion, proposal density and
both physical kernels. The 3,072 local optimizer outcomes, two failing stack
gates and initial background-reference-check failure remain included.

`verify_matched_information.py` reconstructs selected variances using complex
Gaussian monomials and mixture likelihoods using direct Euclidean residuals.
`verify_adaptive_pose_integration.py --dataset DATASET` independently checks
every proposal/importance summary and selected direct cell sums. Verification
scripts intentionally refuse to overwrite their original records; for a new
audit, change only their output destination in a separate checkout and retain
both source hashes. Do not delete archived results to make a rerun look fresh.

All three original runners require committed source/protocol and an unused
output directory. Their commands use `PYTHONPATH=src` and the project Python:

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/analyze_matched_information.py --dataset 10028
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 .venv/bin/python scripts/analyze_matched_haar_likelihood.py --dataset 10028
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/probe_adaptive_pose_integration.py --dataset 10028
```

Use an unused output version for actual reruns; repeat for 10049/10076 without
choosing favorable seeds. The source contains all fixed seeds, indices and
thresholds. Image indices are development data, not a new independent
statistical test cohort. The gates concern the declared Haar/one-CTF/unit-noise
simulator. Background spectra and recorded metadata document why these
assumptions cannot simply be asserted for the experimental images.

The corner-spectrum NPZ files contain patch energies and powers; the complete
second moments and centered covariances are in their summary. The author-code
population replay includes outcomes and input hashes; its upstream likelihood
arrays remain in the pinned public author repository. Consult the separate
replay protocol for its isolated dependency environment and unchanged kernels.
