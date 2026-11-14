# CARMA IMU lever arm: off-CoG installation study

**Status: experimental research.**  The CARMA3 and TFG-CARMA3 filters and their
sketches are research code; `Kalman3D_Wave_OU_III` / `SeaStateFusionFilter_OU_III`
remain the shipping estimator, and nothing here changes it.  All results are
simulation results.

## What was added

An IMU rigidly mounted at body offset `r` from the vessel centre of gravity
(CG) measures

    a_IMU = a_CG + alpha x r + omega x (omega x r),

a tangential and a centripetal term on top of the CG specific force.  Both
are deterministic and both are correlated with attitude.

* `SeaStateFusionFilter_CARMA3::setImuLeverArmBody(r)` (one wrapper template,
  so it serves CARMA3 and TFG-CARMA3) removes the term from the measured
  specific force **at the wrapper's input**, before the vibration guard.  The
  startup tilt proxy, the magnetic reference it learns, the tuner, the
  wave-direction estimator and the inner filter all see CG specific force.
  This is what the sketches' `IMU_LEVER_ARM_X_M / _Y_M / _Z_M` enable
  (body frame: x forward, y starboard, z down; zero by default, which is the
  filter without the correction, bit for bit).
* The model (`src/kalman_common/ImuLeverArm.h`) is the deployable model of the
  OU-III study (`docs/ou-iii-imu-lever-arm.md`, its `gyro` arm): the measured
  rate through a critically damped two-pole low-pass at 15 Hz
  (`setImuLeverArmCutoffHz`), the angular acceleration as the causal
  second-order backward difference of what is left, both terms evaluated on
  the filtered rate.
* Both inner filters also have OU-III's placement,
  `Kalman3D_Wave_CARMA3::set_imu_lever_arm_body` and
  `Kalman3D_Wave_TFG_CARMA3::set_imu_lever_arm_body`: the term enters the
  accelerometer measurement model only, with its gyro-bias Jacobian (TFG-CARMA3
  maps it into the invariant coordinates as `H_beta_g = -Rhat J_L Rhat^T`, its
  one estimate-dependent block).  Their default model is the band-limited one;
  `set_imu_lever_arm_model(LeverArmModel::OuIII)` selects OU-III's own model
  (`Kalman3D_Wave_OU_III::set_imu_lever_arm_body` term for term: a one-pole
  low-pass, time constant 0.05 s, on the first difference of the rate, and
  `d alpha / d b_g = -k_alpha I` in the bias Jacobian).  Both exist for the
  comparison below and for using an inner filter on its own; do not combine
  either with the wrapper's correction.

`tests/kalman_carma/carma_lever_arm-test` checks the model (tangential and
centripetal terms within the low-pass lag, the bias Jacobian against a central
difference, OU-III's recursion for the ablation), that a zero offset is bit
for bit the filter without one (both inner filters and both wrappers), and
that a rigid body rocking with the IMU 1 m above and 0.5 m ahead of the CG is
recovered by both placements of the default model.  With a zero offset the
gated replays `carma3-sim` and `tfg_carma3-sim` reproduce their 240 000-sample
timeseries byte for byte.

## Method

The installation runs in the simulator exactly as in the OU-III study: the
rigid-body term is added to the accelerometer truth of the eight scored
records before sensor corruption (`W3D_IMU_LEVER_ARM_M`).  The filters are
CARMA3 and TFG-CARMA3 as `carma3-sim` (A41) and `tfg_carma3-sim` (T41) run
them, through `tests/kalman_carma/carma_lever_arm-sim`, a separate binary so
the gated replays keep their arithmetic.  Noise realization, configuration,
startup, adaptation and scoring (trailing 900 s of each 1200 s record) are the
same in every arm:

| Arm | Correction |
| --- | --- |
| baseline | IMU at the CG |
| unmodeled | none |
| exact | the simulator removes the installed term using truth kinematics (oracle) |
| wrapper | `setImuLeverArmBody`, the deployed correction |
| mekf | ablation: the same model in the inner filter's accelerometer update (OU-III's placement) |
| mekf_ou3 | ablation: OU-III's own model in OU-III's placement |

Directions are forward, starboard and above the CG; distances 10, 20, 30 cm
(the OU-III study's), 50 cm, 1 m and 2 m.  Every number is an RMS divided by
the same filter's CG baseline, either pooled over the eight seas or for the
worst single sea against its own baseline.  A sweep varies the
wrapper's derivative corner from 1 to 100 Hz with the IMU 1 m above the CG,
where the installed term is largest.

`tools/carma_lever_arm_study.py` runs it and writes
`reports/results/carma_lever_arm_study/`: the per-run and pooled CSVs, the
report with every cell (`lever_arm_report.md`), five figures and a manifest
with data provenance and SHA-256 hashes.  `.github/workflows/carma-lever-arm-study.yml`
re-runs it on changes to the filters and fails if any published ratio moves
by more than 0.01; `tests/validation/test_carma_lever_arm_study.py` checks the
committed evidence and the table below against each other.

## Results

<!-- BEGIN carma_lever_arm_headline -->
**CARMA3** (RMS / CG baseline: eight seas pooled, and the worst single sea)

| IMU | Uncorrected 3-D | worst sea | Corrected 3-D | worst sea | Corrected tilt, worst sea | Corrected yaw | Inner-filter tilt, worst sea | OU-III model tilt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 30 cm forward | 1.00 | 1.27 | 1.000 | 1.003 | 1.01 | 1.02 | 1.04 | 1.34 |
| 30 cm starboard | 1.02 | 1.07 | 1.000 | 1.000 | 1.02 | 0.98 | 1.02 | 1.05 |
| 30 cm above | 1.04 | 1.23 | 1.000 | 1.001 | 1.02 | 0.97 | 1.03 | 1.01 |
| 100 cm forward | 1.05 | 2.29 | 1.001 | 1.016 | 1.03 | 1.14 | 1.11 | 2.83 |
| 100 cm starboard | 1.09 | 1.75 | 1.000 | 1.002 | 1.03 | 0.99 | 1.04 | 1.90 |
| 100 cm above | 1.20 | 2.41 | 1.002 | 1.005 | 1.05 | 1.20 | 1.08 | 1.04 |
| 200 cm forward | 1.20 | 3.81 | 1.008 | 1.061 | 1.05 | 1.21 | 1.20 | 6.16 |
| 200 cm starboard | 1.28 | 3.11 | 1.005 | 1.031 | 1.00 | 1.07 | 1.01 | 7.18 |
| 200 cm above | 1.53 | 4.41 | 1.008 | 1.033 | 1.41 | 1.67 | 2.11 | 1.06 |

**TFG-CARMA3** (RMS / CG baseline: eight seas pooled, and the worst single sea)

| IMU | Uncorrected 3-D | worst sea | Corrected 3-D | worst sea | Corrected tilt, worst sea | Corrected yaw | Inner-filter tilt, worst sea | OU-III model tilt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 30 cm forward | 1.00 | 1.27 | 1.000 | 1.003 | 1.00 | 1.01 | 1.04 | 1.34 |
| 30 cm starboard | 1.02 | 1.07 | 1.000 | 1.000 | 1.01 | 0.99 | 1.01 | 1.05 |
| 30 cm above | 1.04 | 1.23 | 1.000 | 1.001 | 1.01 | 0.98 | 1.02 | 1.01 |
| 100 cm forward | 1.05 | 2.29 | 1.001 | 1.016 | 1.02 | 1.05 | 1.11 | 2.54 |
| 100 cm starboard | 1.09 | 1.75 | 0.999 | 1.002 | 1.01 | 0.98 | 1.02 | 1.62 |
| 100 cm above | 1.20 | 2.41 | 1.002 | 1.004 | 1.02 | 1.01 | 1.06 | 1.03 |
| 200 cm forward | 1.20 | 3.81 | 1.007 | 1.062 | 1.05 | 1.07 | 1.09 | 3.51 |
| 200 cm starboard | 1.28 | 3.10 | 1.004 | 1.031 | 1.00 | 0.99 | 1.00 | 4.20 |
| 200 cm above | 1.53 | 4.41 | 1.007 | 1.033 | 1.06 | 1.18 | 1.13 | 1.05 |

Derivative corner of the corrected arm, IMU 100 cm above the CG (3-D RMS / CG baseline):

| Corner [Hz] | CARMA3 | TFG-CARMA3 |
| ---: | ---: | ---: |
| 1 | 1.026 | 1.026 |
| 2 | 1.007 | 1.007 |
| 5 | 1.002 | 1.002 |
| 10 | 1.001 | 1.001 |
| 15 | 1.002 | 1.002 |
| 25 | 1.009 | 1.007 |
| 50 | 1.024 | 1.023 |
| 100 | 1.024 | 1.024 |
<!-- END carma_lever_arm_headline -->

The table is rendered from `reports/results/carma_lever_arm_study/lever_arm_summary.csv`
by `headline_markdown()`; `lever_arm_report.md` has every cell, every arm and
the removed fractions.

### What an off-CoG IMU costs

On the moving-boat records, pooled uncorrected 3D RMS ranges from 1.001–1.039
of the CG baseline at 30 cm, 1.047–1.197 at 1 m and 1.200–1.532 at 2 m,
across both filters and all three directions. The worst individual sea reaches
1.271, 2.410 and 4.412 respectively. Pooled tilt at 2 m ranges from
1.014–1.140 of baseline, with a worst-sea ratio of 1.267.

### What the correction recovers

At the wrapper input, the deployed correction keeps pooled 3D RMS within
0.2% of the CG baseline through 1 m and 0.8% at 2 m. Worst-sea 3D ratios
reach 1.016 at 1 m and 1.062 at 2 m. Worst-sea tilt ratios reach 1.048 and
1.407 respectively. The exact-kinematics oracle reproduces 3D and tilt to
within 0.1%; the deployable correction still includes gyro noise and filter lag.
Pooled yaw can reach 1.197 at 1 m and 1.670 at 2 m: displacement recovery
does not imply equally good heading recovery.

In the derivative sweep, 5–15 Hz gives pooled 3D ratios of 1.0011–1.0019
for the IMU 1 m above the CG. The ratios rise to approximately 1.026 at
1 Hz and 1.024 at 100 Hz. The deployed 15 Hz setting is retained; this is
one installation and one set of records, not a universal optimal cutoff.

### Why at the wrapper's input

Applying the same correction inside the inner filter gives worst-sea 3D
ratios of 1.064 at 1 m and 1.063 at 2 m. Its worst-sea tilt ratio,
however, reaches 2.114 at 2 m, versus 1.407 for the wrapper placement;
pooled yaw reaches 2.997 versus 1.670. The wrapper placement corrects the
specific force before startup, adaptation and the inner filter consume it.
The matrix compares complete placements and does not isolate those consumers.

### OU-III's model

The OU-III-model ablation inside these CARMA filters gives pooled tilt ratios
up to 2.825 at 1 m and 7.178 at 2 m; worst-sea tilt reaches 9.969 and
pooled yaw reaches 38.816 at 2 m. Its pooled displacement score can remain
close to baseline despite severe attitude degradation. It uses a one-pole
filtered difference and the alpha bias Jacobian described above. The matrix
does not isolate their separate effects. It runs CARMA3 and TFG-CARMA3,
not OU-III itself; the separate OU-III study provides that estimator's results.

## Limits

Simulation only: one 28 ft hull and its RAOs, one noise realization per sea,
a rigid hull, and an offset known exactly.  A mis-measured offset was not
studied; a wrong sign doubles the term instead of removing it.  No hardware,
on-target timing or sea data stands behind any of this, and the CARMA filters
remain research estimators.

The current moving-boat matrix was regenerated with g++ 13.3 (Ubuntu 24.04) at source commit `5b66d7a`, after the CARMA wrapper stopped attenuating wave variance on still-water detection; `manifest.json` records every published file hash. The simulator is built `-march=native`, so CI compares ratios to a 0.01 tolerance, and only when the committed and current input hash maps match.
