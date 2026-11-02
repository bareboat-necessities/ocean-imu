# Calibration requirements audit

The source implements the requested calibration extensions. The complete
delivery is **not yet qualified**: physical operator/device validation,
ESP32-S3 runtime and stack measurements, and qualification of the changed
magnetic measurement policy remain open. Passing synthetic replays does not
close those requirements. Build/test results and reproduction commands are
in [calibration-validation.md](calibration-validation.md).

| Requirement | Implementation and reuse | Verification / remaining work |
| --- | --- | --- |
| 1. Review and reuse | Existing capture, numerical utilities, wizard, correction and Preferences store were extended. The accelerometer and magnetometer capture/fitter headers are byte-for-byte identical to reference main `d2d34116e3a1a1ca7e26cde49e0f2473dcb14a08`. | Retention independently checked against the reference files; applicable numerical and persistence tests retained. |
| 2. Robust accelerometer | `AccelCalCapture.h` and `AccelCalFit.h` retain multi-orientation quiet blocks, outlier rejection, PolarSPD fitting, conditioning, physical bounds, thermal evidence, residuals, cross-validation and stored-float validation. | Existing 30-seed campaign and release-reference checks; no simpler INSLIB replacement. Mounting is applied separately in `RuntimeCals`. |
| 3. IMU-TK gyro extension | `GyroCalCapture.h` supplies stationary bias/endpoints. `GyroCalFit.h` integrates measured intervening rates against calibrated gravity endpoints and fits the existing `GyroCalibration::S`. `GyroTurnCapture.h` controls acquisition. | Analytic sensitivity checks, exact/noisy full-matrix recovery, final-accelerometer endpoints, independent validation rejection, missing/frozen data and guided controller replays. Independent real device rotations remain required. |
| 4. Magnetometer and physical frames | Existing hard/soft-iron fitting is retained. `solveAlignment()` adds an undamped observability check and independent dip validation. `FrameAlignment.h` solves known-direction registration. `RuntimeCals` composes intrinsic, magnetometer-to-IMU, IMU-to-housing and housing-to-vessel transforms once. | Noncommuting transforms, reflection/degeneracy rejection and held-out magnetic poses are tested. MOUNT accepts two measured rotation matrices over serial, or an explicit confirmation that the three frames coincide. It does not infer unknown installation yaw from gravity. Physical fixture/installation measurements remain required. |
| 5. Existing wizard evolution | Existing UI, acquisition, quality, tasks, retry and verified save paths remain. Added rotation, mounting and separate noise stages. Initial gyro set is 12 turns, with optional blocks of 4 up to 20. | Production controller/worker/UI replays cover slow LCDs, orientation, pauses, late hold jitter and local retry. Instructions and start prompts are untimed; completed turns remain on retry. Initial/endpoint holds use separate rate profiles; hand-supported 4/6/8 Hz micro-tremor is replayed with coupled acceleration and rate. Placement gets its own full 3.5 s window; misplaced/unsettled endpoints wait on the same target instead of retrying the turn. Real operator acceptance remains open. |
| 6. New format, no old-record fallback | One 680-byte schema 6 record/key stores all intrinsics, gyro matrix, three rotations, quality/noise metadata, identity and CRC. Existing NVS verified writes and verified current-record recovery are reused. Only schema 6 is loaded. Full recalibration starts without previous slopes or gyro-reference coefficients. | Old keys/versions, truncation, CRC corruption, invalid rotations and coefficient/quality bindings are tested. A failed or cancelled new acquisition cannot publish a partial replacement. Preserving a valid current record after an aborted write is storage recovery, not old-format migration. |
| 7. Measured noise | Separate 40-minute sensor-clock acquisition, bounded-memory dyadic Allan deviation, timing/frozen/motion/temperature/magnetic checks, qualified white/plateau/walk fits and explicit unidentifiable terms. Runtime reset paths configure measured gyro process noise and accelerometer scatter. | White and mixed stochastic fixtures, short-data rejection, wraparound and the production worker's complete simulated acquisition/loss/cancel paths pass. Record-to-runtime tests verify units and identified-axis masks. A real warmed-up 40-minute trace has not been acquired. Marine OU/CARMA priors remain distinct. |
| 8. Magnetic innovation and recovery | `kalman_common/MagneticInnovationGate.h` supplies correlated 3x3 vector NIS, chi-square thresholds, covariance inflation, gross rejection and recovery hysteresis to OU-II/III, TFG, both CARMA variants and qMEKF. Existing residuals/Jacobians, Joseph updates and injections remain in each filter. Trusted-reference refinement also requires accepted healthy magnetic service. | Correlated Gaussian and frame-invariance checks; all six cores reject gross field-strength changes even with large measurement covariance; 90-degree disturbance/recovery and covariance tests for all six cores; actual compass/facade reference-learning regressions. Changed magnetic-service/covariance behavior still needs stability qualification; archived evidence fails closed. |
| 9. Cleanup and coherent API | `RuntimeCals` is the single correction path; `runImuCalWizard` is the application entry and the store accepts one schema. Gyro fitting extends the existing calibration object rather than introducing another correction layer. The accelerometer-only UI was retired because it could leave dependent gyro/alignment qualifications stale. | Driver mapping and M5 offset clearing remain upstream. Numeric solvers use existing Eigen; no second optimization dependency. Standalone library/test helpers are not additional active firmware pipelines. |
| 10. Tests, comparison and ESP32 cost | Existing tests and quantitative baseline comparisons are retained, with new numerical, persistence, controller, worker and fusion regressions. `atomS3R_calibration_benchmark` and production task logs report runtime, allocations and stack high-water marks. | Host numerical checks and representative crossbuilds are evidence, but a benchmark harness/static frame estimate is not a measured device benchmark. Full-repository validation and current firmware build status must be reported explicitly. Actual ESP32-S3 runtime, heap/PSRAM and task stack peaks remain open. |

## Mathematical changes

The bias-only gyro workflow was insufficient to estimate scale and axis
misalignment. The existing correction matrix is now determined from diverse
gravity-referenced rotations, with a full-rank nine-parameter information
test and a frozen independent validation set. The shorter wizard plan changes
sample redundancy, not conditioning, coefficient plausibility or residual
limits. Additional turns are offered when the measured geometry is insufficient.

Magnetic alignment retains the constant-dip model and its intended small
registration range. A numerical ridge cannot make unobservable geometry
observable, so rank/conditioning is checked before regularization. Physical
mounting rotations are stored separately from intrinsic SPD matrices.

The shared magnetic policy changes which measurements are applied and the
effective covariance. Existing stability certificates bound another operator;
their fingerprints cannot be updated merely to make tests pass. This is a
qualification failure (class E), not a demonstrated instability result.

## Completion evidence still needed

1. Complete the repository's required native and MCU checks for the final source.
2. Run the wizard on the actual AtomS3R, including USB-up turn 2, deliberate
   mistakes, long reading pauses, retries and a complete independently validated
   calibration. Record accuracy against independent rotations and fixture frames.
3. Acquire an uninterrupted warmed-up 40-minute stationary trace and retain the
   raw data, integrity results, Allan table and model identifiability report.
4. Run the ESP32-S3 benchmark and full wizard at the documented board settings;
   retain fit times, sample timing, heap/PSRAM and actual stack high-water marks.
5. Qualify the new magnetic health/inflation/reference state and actually applied
   service in the shipping stability/evidence path before promoting it for release.

No physical result or device timing is inferred from a synthetic test, and no
whole-system completion claim follows from the implementation table alone.
