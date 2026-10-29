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
| 30 cm forward | 1.01 | 1.35 | 0.999 | 1.005 | 1.00 | 1.00 | 1.09 | 1.23 |
| 30 cm starboard | 1.02 | 1.42 | 1.001 | 1.002 | 1.01 | 1.00 | 1.03 | 1.03 |
| 30 cm above | 1.03 | 1.37 | 1.001 | 1.002 | 1.01 | 0.98 | 1.06 | 0.98 |
| 100 cm forward | 1.09 | 2.50 | 1.000 | 1.013 | 0.99 | 1.03 | 1.35 | 2.00 |
| 100 cm starboard | 1.14 | 2.81 | 1.007 | 1.008 | 1.01 | 1.02 | 1.03 | 1.25 |
| 100 cm above | 1.21 | 3.00 | 1.002 | 1.015 | 1.04 | 1.03 | 1.37 | 0.95 |
| 200 cm forward | 1.32 | 4.27 | 1.017 | 1.068 | 1.05 | 1.06 | 1.71 | 4.03 |
| 200 cm starboard | 1.42 | 4.98 | 1.033 | 1.038 | 1.02 | 1.06 | 0.99 | 4.20 |
| 200 cm above | 1.64 | 5.61 | 1.021 | 1.065 | 1.08 | 1.18 | 1.87 | 0.95 |

**TFG-CARMA3** (RMS / CG baseline: eight seas pooled, and the worst single sea)

| IMU | Uncorrected 3-D | worst sea | Corrected 3-D | worst sea | Corrected tilt, worst sea | Corrected yaw | Inner-filter tilt, worst sea | OU-III model tilt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 30 cm forward | 1.01 | 1.35 | 0.999 | 1.005 | 1.01 | 1.00 | 1.09 | 1.24 |
| 30 cm starboard | 1.02 | 1.42 | 1.001 | 1.002 | 1.00 | 1.00 | 1.03 | 1.03 |
| 30 cm above | 1.03 | 1.37 | 1.001 | 1.002 | 1.01 | 0.99 | 1.05 | 0.98 |
| 100 cm forward | 1.09 | 2.50 | 0.999 | 1.013 | 1.02 | 1.01 | 1.35 | 2.00 |
| 100 cm starboard | 1.14 | 2.81 | 1.007 | 1.008 | 1.00 | 1.01 | 1.03 | 1.22 |
| 100 cm above | 1.21 | 3.00 | 1.002 | 1.015 | 1.02 | 0.99 | 1.35 | 0.95 |
| 200 cm forward | 1.32 | 4.27 | 1.016 | 1.067 | 1.10 | 1.01 | 1.72 | 2.48 |
| 200 cm starboard | 1.42 | 4.98 | 1.032 | 1.037 | 1.05 | 1.01 | 0.98 | 2.54 |
| 200 cm above | 1.64 | 5.61 | 1.021 | 1.065 | 1.06 | 1.05 | 1.84 | 0.94 |

Derivative corner of the corrected arm, IMU 100 cm above the CG (3-D RMS / CG baseline):

| Corner [Hz] | CARMA3 | TFG-CARMA3 |
| ---: | ---: | ---: |
| 1 | 1.129 | 1.129 |
| 2 | 1.030 | 1.030 |
| 5 | 1.007 | 1.007 |
| 10 | 1.003 | 1.003 |
| 15 | 1.002 | 1.002 |
| 25 | 1.021 | 1.021 |
| 50 | 1.114 | 1.114 |
| 100 | 1.106 | 1.106 |
<!-- END carma_lever_arm_headline -->

The table is rendered from `reports/results/carma_lever_arm_study/lever_arm_summary.csv`
by `headline_markdown()`; `lever_arm_report.md` has every cell, every arm and
the removed fractions.

### What an off-CoG IMU costs

Uncorrected, the installation costs displacement in every direction and at
every distance, and most with the IMU above the CG.  Pooled over the eight
seas it is small to 30 cm (1.01-1.03x) and grows to 1.09-1.21x at 1 m and
1.32-1.64x at 2 m.  The pooled RMS is dominated by the largest seas, where the
installed term is smallest against the wave signal; in the calm seas
(`H_s` 0.27 m) the same installations cost 1.35-1.42x at 30 cm, 2.5-3.0x at
1 m and 4.3-5.6x at 2 m (`carma_lever_arm_sea_state.svg`).  The two filters
agree to the second decimal: the installed term reaches both through the same
specific force, and they are first-order equivalent
(`docs/tfg-carma3-research.md`).

Attitude is not consistently degraded: per sea the uncorrected tilt RMS moves
in both directions, and pooled it stays within 0.89-1.01x.  The installed
force combines with each record's existing residual; a lower score is not
evidence that the installation helps.  Yaw moves by tens of percent per sea in
every arm, the uncorrected one included (0.46-2.0x with the IMU 2 m above the
CG): heading depends on the startup magnetic acquisition, which any change to
the accelerometer stream perturbs.  Read the yaw column as that sensitivity,
not as a property of the lever arm.

### What the correction recovers

The deployed correction returns displacement to within 0.7 % of the CG
baseline to 1 m (worst sea 1.5 %) and within 3.3 % at 2 m (worst sea 6.8 %),
in all three directions and both filters; worst-sea tilt stays within 4 % to
1 m and 10.3 % at 2 m.  The oracle returns to the baseline to within 0.1 %, so
what is left is the model's own: differentiated gyro noise, which grows with
`r`, and the lag of the 15 Hz low-pass.  Pooled yaw stays within 3.3 % to 1 m; at
2 m it is 1.06-1.18x for CARMA3 and 1.01-1.05x for TFG-CARMA3, inside the
acquisition sensitivity above.

The derivative corner is two-sided, as in the OU-III study: 5-15 Hz is a flat
basin (1.002-1.007x with the IMU 1 m above the CG), 1-2 Hz loses to lag and
25-100 Hz to differentiated noise.  15 Hz is the OU-III study's value and is
not refitted here.

### Why at the wrapper's input

The same model in the inner filter's accelerometer update (OU-III's placement)
recovers pooled displacement about as well, but not in the calm seas: with the
IMU forward of the CG the worst sea keeps 1.15x of displacement error at 1 m and
1.39x at 2 m, and tilt reaches 1.35-1.37x at 1 m forward or above and
1.7-1.9x at 2 m, where the input placement stays within 1.04x and 1.10x.  Corrected inside the filter, the rest
of the wrapper still sees the IMU's motion: the startup proxy learns tilt and
the magnetic reference from it, and the tuner sets the CARMA priors and the
accelerometer noise weighting from it.  The calm seas are where that term is
largest against the wave signal, which is where the cost appears.  The matrix
does not separate those consumers; it shows that removing the term before all
of them recovers what removing it in one does not.

### OU-III's model

OU-III's own model, in OU-III's placement, is the one arm that is worse than
no correction at all.  Pooled tilt is 1.23x at 30 cm forward, 2.0x at 1 m and
2.5-4.2x at 2 m forward or starboard, with yaw several times the baseline at
2 m, and it removes only about half of the forward displacement excess.  It
differs from the inner-filter arm in two ways: the angular acceleration is a
one-pole low-pass of the first difference of the unfiltered rate, whose
differentiated noise does not roll off, and the bias Jacobian carries
`d alpha / d b_g = -k_alpha I` (19 per second at 200 Hz and 0.05 s), which
couples every accelerometer residual into the gyro bias at `k_alpha |r|`.  A
constant gyro-bias error has no effect on a differenced rate, so that term
models a sensitivity that does not exist.  The matrix does not separate the
two; a one-sea diagnostic during development (JONSWAP `H_s` 4 m, 1 m arms,
CARMA3) did: without the `k_alpha` term, OU-III's estimator kept tilt at the
CG level on the forward and starboard arms; with it, pitch RMS tripled on the
forward arm.

This study does not run OU-III, so OU-III remains untested and the study says
nothing directly about its numbers.  It does show this model degrading
attitude in both CARMA filters it ran (CARMA3's multiplicative attitude filter
and TFG-CARMA3's two-frame invariant one), and the OU-III sketch exposes the
same model through its `IMU_LEVER_ARM_*_M`.  An OU-III run of this matrix
belongs before anyone enters a nonzero offset there.

## Limits

Simulation only: one 28 ft hull and its RAOs, one noise realization per sea,
a rigid hull, and an offset known exactly.  A mis-measured offset was not
studied; a wrong sign doubles the term instead of removing it.  No hardware,
on-target timing or sea data stands behind any of this, and the CARMA filters
remain research estimators.
