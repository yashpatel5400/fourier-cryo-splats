# Registered-reference target sensitivity and experimental cubic application

30 September 2026. All previously declared registered-target checks completed
on three stacks in 8.76 seconds. The follow-up experimental application also
completed, in 0.51 seconds. Both are post-outcome development, with all previous
results preserved. They do not establish experimental density coverage.

Registration changes the simulated density generator, not the estimator,
target, noise, bias bound, or existing interval. All twelve original nominal
centers replay to the declared tolerance. All 120 conditional target/pose
records and all nine cubic records are retained. Fixed-pose minimum powers
across each target's scenarios are at least .99999 on 10028 and 10076, and
.95632 on 10049. Every original one-/two-degree fixed-weight audit still has
negligible sign power. The 48 unchanged experimental intervals contain all
registered, pilot-amplitude-scaled reference targets; the existing six
fixed-pose and four shift-only zero exclusions remain unchanged. Inclusion
in these broad intervals is an agreement check, not coverage evidence.

On the already selected 10049 region 1 (target Gaussian SD 20 A, one-degree /
0.5 A joint pose class), minimum simulated sign power changes as follows:

| Completed estimator | Original reference frame | Registered frame |
| --- | ---: | ---: |
| Original cubic design | .00654 | .95725 |
| Coordinate-metric design | .000000773 | .26042 |
| Reduced-subspace design | .00663 | .95770 |

The original and reduced-subspace weights are identical to numerical precision;
their small width/power differences arise from separate numerical audits.
The frame change explains the larger power, not an estimator improvement.
The large unresolved optimization gaps remain. These supplied-noise calculations
are conditional implementation studies, not real-data calibration.

The separately declared experimental application retains all three designs,
reuses the existing 128 calibration exposure representatives, and replays the
fixed estimator's observed center, calibrated SD and half-width before computing
new outcomes. Calibration permits arbitrary jointly Gaussian within-exposure
noise dependence: the 128 inference particles occupy only 26 exposures, with
up to twelve particles in an exposure. This conservatism is part of the
declared procedure, not proof that such dependence occurs in the data.

| Estimator | Observed center | Calibrated SD upper | Half-width | Relative to no data | Excludes zero |
| --- | ---: | ---: | ---: | ---: | --- |
| Original cubic | 2.46178 | 1.01002 | 4.30107 | .35374 | No |
| Coordinate metric | 2.48844 | 1.01167 | 4.82603 | .39692 | No |
| Reduced subspace | 2.46178 | 1.01002 | 4.29998 | .35365 | No |

The original cubic SD upper is 4.55146 times its supplied simulation SD.
All three intervals contain both the original reference target (.71012) and
the registered target (2.48846); none falls back to no data. The same target's
fresh fixed-pose interval has half-width 2.96825 about 2.45941. Thus this
application does not establish useful rotationally robust sign detection.
It identifies noise magnitude and exposure concentration as obstacles to
transferring the registered simulation result to experimental images.

The application allocates Gaussian tail error (.045-1e-6)/12, calibration
error .005/12, and spectral error 1e-6/12 to each procedure. A union bound
for these three predeclared alternatives is .0125 under the fixed-design and
noise assumptions. Reuse of inspected data and unverified experimental input
conditions prevents asserting that guarantee empirically. No minimum-width
selection over historical alternatives is reported.

Records: registered-target-sensitivity-v1 (three geometry files and summary),
cubic-experimental-application-v1/10049.json, exact array hashes, archived
source bytes, and the two separately committed protocols. The v0.6 release
precedes these results and retains its historical manuscript unchanged.
