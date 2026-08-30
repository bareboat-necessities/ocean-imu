# OU-III IMU installation lever-arm study

This study measures the error introduced when the IMU is not mounted at the vessel center of gravity (CG), and how much of that error is removed by modeling the same installation vector inside the filter.

## Physical model

For a rigid IMU fixed at body-frame vector \(r\) from the CG, its translational acceleration differs from the CG by

\[
a_{\rm IMU}=a_{\rm CG}+\dot\omega\times r+\omega\times(\omega\times r).
\]

## Where the two halves run

The installation and the filter's model of it sit on opposite sides of the sensor, so the simulator implements them there rather than by rewriting records. Both stages are off unless `W3D_IMU_LEVER_ARM_M` is set, and the default realization is bit-for-bit unchanged.

| Variable | Meaning |
| --- | --- |
| `W3D_IMU_LEVER_ARM_M` | `"x,y,z"` body z-up offset of the IMU from the CG, metres. Enables both stages. |
| `W3D_IMU_LEVER_ARM_MODEL` | `none` (default), `exact`, `gyro`, or `estimated`: the lever-arm model the filter's input stage applies. |
| `W3D_IMU_LEVER_ARM_CUTOFF_HZ` | Corner of the two-pole low-pass the `gyro` model runs ahead of its derivative. Default 15 Hz. |
| `W3D_IMU_LEVER_ARM_PRIOR_STD_M` | Per-axis prior std the `estimated` arm starts from. Default 0.5 m. |
| `W3D_IMU_LEVER_ARM_TRACE_SEC` | Period of the `LEVERARM` convergence trace on stderr. Off by default. |
| `W3D_IMU_LEVER_ARM_ALPHA_TAU_S` | Overrides the filter's angular-acceleration smoothing time constant. |

- **Installation stage** — runs on the record's noiseless truth, before sensor corruption. It applies the rigid-body term above using the record's own body rate and a causal second-order backward difference for \(\dot\omega\). This is what a sensor at \(r\) physically measures.
- **Model stage** — runs after sensor corruption, immediately before fusion. This is the stage a firmware implementation owns. Nothing downstream of it knows the IMU is off the CG.

The simulator prints `IMU_LEVER_ARM` when the stage is installed and `IMU_LEVER_ARM_RESULT` when the record finishes, reporting the RMS of the term the installation injected and of the residual the filter's model left behind.

## Matrix

The full run uses the eight standard stationary JONSWAP / PM-Stokes records and three canonical installation directions: x (athwartships), y (fore-aft), and z (vertical). Each direction is tested at 10, 20, and 30 cm from the CG. The CG displacement/velocity/acceleration truth, attitude, gyro truth, magnetometer, sensor-noise realization, OU-III tuning, startup, vibration guard, pseudo-measurements, and scoring are identical in every arm.

Five arms are reported:

- **baseline**: IMU at the CG.
- **unmodeled**: off-CG specific force reaches OU-III with no filter-side model.
- **gyro**: the deployable model. It sees only the corrupted rate and reconstructs \(\dot\omega\) with a two-pole low-pass at 15 Hz followed by a causal second-order difference.
- **exact**: the oracle. It removes the same term using the record's own angular kinematics, and bounds what any lever-arm model can recover.
- **estimated**: self-calibrating. Nothing is removed before fusion; OU-III carries \(r\) as three of its own states, seeded at zero with a 0.5 m per-axis prior, and has to find the installation from the motion.

### The self-calibrating arm

Both `exact` and `gyro` are handed \(r\). `gyro` replaces the *kinematics* with a reconstruction from the measured rate, not the *geometry*. The `estimated` arm removes that input, which is possible because the injected term is linear in \(r\):

\[
\dot\omega\times r+\omega\times(\omega\times r)=M(\omega,\dot\omega)\,r,
\qquad M=[\dot\omega]_\times+[\omega]_\times[\omega]_\times,
\]

so \(\partial/\partial r\) of the term is \(M\) itself — exact, and cheap enough to sit in the same accelerometer update that corrects attitude and the biases. `Kalman3D_Wave_OU_III` grows three states behind a fourth template parameter (`with_lever_arm`, default false), and the lever arm is read back through `get_imu_lever_arm_body()` with a covariance from `get_imu_lever_arm_covariance()`.

Two properties bound what this can do, and both are measured rather than asserted:

- **Observability.** \(M\) annihilates any \(r\) parallel to the instantaneous rotation axis, so a vessel turning about one fixed axis never reveals the component of \(r\) along it. A seaway helps because the axis wanders, but *how much* it wanders is a property of the sea state. `tests/kalman_ou_iii/lever_arm_estimation-test.cpp` pins down both halves: rich three-axis rotation recovers a 39 cm arm from a zero prior to a few millimetres, single-axis rotation recovers only the perpendicular part and keeps its prior width along the axis.
- **Over-confidence.** The estimate's error on the wave records is set by model error — residual tilt, the accelerometer bias competing for the non-zero mean of the centripetal term, an \(\dot\omega\) differentiated from a noisy rate — none of which is in the innovation covariance. The filter reports millimetres of uncertainty on an estimate that is centimetres or decimetres out. `set_lever_arm_rate_noise_modelled()` carries the white rate noise of the lever term into the innovation covariance (using the same matrix that already carries the gyro-bias coupling, since \(\omega=\mathrm{gyr}-b_g\)), which is correct but small next to the model error. **The covariance describes the conditioning of the regression, not the accuracy of the answer.**

The practical use is therefore as a check on a survey rather than a replacement for one: seed \(r\) with the measured value and a prior of a few centimetres, and the filter will move an arm that was entered wrong and hold one that was entered right. Where the lever arm can be measured, `gyro` remains the better instrument — it recovers the whole penalty without having to discover the geometry.

A separate sweep varies the `gyro` model's derivative band over 1-100 Hz on a 30 cm fore-aft arm. The band is the model's one design parameter and it is two-sided: too narrow and the low-pass phase lag misaligns a correction whose amplitude is already right, too wide and differentiated gyro noise exceeds the term being removed.

## Outputs

`tools/ou3_lever_arm_study.py` writes:

- `lever_arm_runs.csv`: one row per sea / axis / distance / modeling arm;
- `lever_arm_summary.csv`: pooled RMS over equal-duration scored windows, plus the injected and residual specific force and the fraction of the unmodeled excess each arm removes;
- `lever_arm_cutoff_runs.csv` / `lever_arm_cutoff_summary.csv`: the derivative-band sweep;
- `lever_arm_report.md`: publication-ready numerical summary and interpretation boundary;
- `ou3_lever_arm_penalty.svg`: 3-D displacement penalty per direction, and its removal by both models;
- `ou3_lever_arm_tilt.svg`: the same for maximum roll/pitch;
- `ou3_lever_arm_mechanism.svg`: the injected term against offset, and the residual each model leaves;
- `ou3_lever_arm_sea_state.svg`: absolute per-sea error at the worst direction and offset;
- `ou3_lever_arm_cutoff.svg`: the derivative-band sweep against both bounds;
- `ou3_lever_arm_calibration.svg`: the self-calibrating arm scored on its calibration — how much of the arm it found, and its error against the uncertainty it reports;
- `manifest.json`: source commit, data provenance, model statement, and SHA-256 hashes.

With `--mirror-doc` the six figures are copied byte-for-byte into `doc/kalman_ou_iii/` for the article; `tools/ou3_lever_arm_tex.py` generates the section's numeric fragment from the two summary CSVs. `tests/validation/test_ou3_lever_arm_article.py` fails if the committed article stops matching the committed evidence.

The full study is run in CI by `.github/workflows/ou3-lever-arm-study.yml` against `oceanography-waves-lib` release `v1.1.3`. The simulator stage has its own unit test in `tests/kalman_ou_iii/imu_lever_arm-test.cpp`, and the estimated arm's filter states in `tests/kalman_ou_iii/lever_arm_estimation-test.cpp`.
