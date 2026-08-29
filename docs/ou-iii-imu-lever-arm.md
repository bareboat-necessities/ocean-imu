# OU-III IMU installation lever-arm study

This study measures the error introduced when the IMU is not mounted at the vessel center of gravity (CG), and the improvement available when the same installation vector is modeled exactly.

## Physical model

For a rigid IMU fixed at body-frame vector \(r\) from the CG, its translational acceleration differs from the CG by

\[
a_{\rm IMU}=a_{\rm CG}+\dot\omega\times r+\omega\times(\omega\times r).
\]

The study derives \(\omega\) from the versioned record's body-rate truth and \(\dot\omega\) by centered finite differences. Only `acc_bx`, `acc_by`, and `acc_bz` are changed. The CG displacement/velocity/acceleration truth, attitude, gyro truth, magnetometer, sensor-noise realization, OU-III tuning, startup, vibration guard, pseudo-measurements, and scoring remain unchanged.

## Matrix

The full run uses the eight standard stationary JONSWAP / PM-Stokes records and three canonical installation directions: x (athwartships), y (fore-aft), and z (vertical). Each direction is tested at 10, 20, and 30 cm from the CG.

Two matched arms are reported:

- **Unmodeled**: OU-III receives the off-CG accelerometer signal as if it were measured at the CG.
- **Exact model**: the same rigid-body term is removed immediately before fusion. This is an ideal upper bound with exact angular kinematics and exact lever vector; normal stochastic sensor corruption remains enabled after the reconstructed input record is formed.

The exact arm is deliberately not a claim about a deployable noisy gyro differentiator. It answers a narrower question: how much of the observed installation penalty is deterministically recoverable if the lever-arm kinematics are known exactly.

## Outputs

`tools/ou3_lever_arm_study.py` writes:

- `lever_arm_runs.csv`: one row per sea / axis / distance / modeling arm;
- `lever_arm_summary.csv`: pooled RMS over equal-duration scored windows;
- `lever_arm_report.md`: publication-ready numerical summary and interpretation boundary;
- `ou3_lever_arm_3d.svg`: pooled 3-D displacement degradation relative to the CG baseline;
- `ou3_lever_arm_tilt.svg`: pooled maximum roll/pitch degradation relative to the CG baseline;
- `manifest.json`: source commit, data provenance, model statement, and SHA-256 hashes.

The full study is run in CI by `.github/workflows/ou3-lever-arm-study.yml` against `oceanography-waves-lib` release `v1.1.3`.
