# Reproducing the post-review method gates

The v0.7.2-dev increment preserves completed negative and positive development
outcomes. It is not a new blinded confirmation, experimental calibration, or
favorable full review. All three full Fable reviews still reject.

Use the tagged repository and its pinned Python environment. Download the
release's `v0.7.2-dev-method-gates.tar.gz` and adjacent manifest. Before extraction:

```bash
.venv/bin/python scripts/verify_uq_release_artifacts.py output/artifacts/v0.7.2-dev-method-gates.tar.gz
```

Extract the verified archive at the repository root. It includes every new
numerical array, all three EMDB input maps, the preceding 4,160-view arrays,
the original moment-screen arrays, and the three local-refinement generator
files. The EMDB maps and original moment arrays are contextual dependencies,
not new downloads or new experimental particles. Earlier raw-particle
reconstruction reproduction still uses the older releases and their manifests.
Regenerating the earlier initial moment/orbit preparation also requires the
metadata, splits and prior target/noise inputs listed in its original record;
the new bundle supplies its exact saved outputs for the subsequent analyses.

The following checks use the saved arrays and need no new EMPIAR download:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_revision5_diagnostics.py --output /tmp/revision5-replay.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/verify_uq_adversarial_witnesses.py --output /tmp/revision5-direct-witnesses.json
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q tests/test_paired_power.py tests/test_power_cone.py tests/test_power_cone_witness.py tests/test_paired_covariance.py tests/test_covariance_cone_bounds.py tests/test_folded_ridge.py tests/test_shift_moments.py tests/test_pose_cone_probe.py
bash scripts/build_paper.sh
```

Choose unused verification output paths: the scripts refuse to overwrite an
existing scientific record. The first verifier replays covariance mixtures,
computes feature ranks/discrimination retention and independently evaluates
all six selected ridge widths using direct sinc integrals. The second evaluates
the four worst continuous-pose counterexamples by physical-cell exponential
sums, independent of the interpolation guide and NUFFT forward evaluator.
Neither is validated interval arithmetic. The focused test run has 32 passes;
it does not replace the earlier historical full-suite result.

For exact original optimization runs, use a clean checkout of the
`git_head` recorded in each result summary, with its input artifacts but
without its output directory. The runners require committed matching source
and refuse to replace prior output. The original sequence is:

| Run | Script | Original source commit |
|---|---|---|
| Compressed covariance | `scripts/probe_uq_paired_covariance.py` | `a07d4fe` |
| Folded-width ridge | `scripts/probe_uq_folded_ridge.py` | `7ff8b11` |
| Full covariance/fresh views | `scripts/probe_uq_full_covariance.py` | `fa6f939` |
| Common shifts/third view set | `scripts/probe_uq_shifted_covariance.py` | `d2f9859` |
| Adversarial continuous poses | `scripts/probe_uq_adversarial_pose.py` | `98d2654` |

The corrected-bound script is a post-audit replay and mass-sensitivity
analysis, not a replacement for its predecessor's outcomes. Likewise, the
current LP implementation adds dual checks and zeroes a degenerate direction;
the earlier compressed outputs retain the source commit that produced them.

The authentic focused review's prompt, all assistant text, raw provider event
stream, model identity, source hashes and response are included. It checks
mathematics from text and did not execute experiments. Its acceptance status
must not be inferred from favorable algebra checks.
