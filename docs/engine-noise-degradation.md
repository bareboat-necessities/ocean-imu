# Engine-noise degradation study

## Purpose

Every other noise source in this repository is electronic: accelerometer and
gyroscope white noise, initial bias, bias random walk, and magnetometer
calibration error.  A cruising sailboat carrying this estimator spends a
substantial share of its passages under power, and an inboard auxiliary diesel
puts mechanical vibration into the same two sensors.

This study measures how far the deployed OU-II, OU-III, and TFG estimators
degrade under that vibration, and identifies which part of the sensor path
causes it.  It is a **sensor-path** study: the vessel's rigid-body response to
the sea, the wave records, and the magnetometer are all unchanged.

## The vibration model

`src/util/W3dSimCommon.{h,cpp}` gains an engine vibration model for the
prescribed auxiliary-vibration stress applied to the 28 ft RAO response: a naturally aspirated
three-cylinder four-stroke diesel on flexible mounts, a 2.6:1 reduction gear,
and a three-blade fixed propeller.  It injects

- crank orders: the half order, the first order, the firing order `n_cyl/2`,
  and the firing harmonics `2x`, `3x`, `4x`;
- driveline lines: shaft rate, propeller blade rate, and twice blade rate;
- a broadband structural floor, high-passed at 5 Hz so the engine contributes
  little directly in the wave band;
- governor hunting (an Ornstein-Uhlenbeck speed deviation of 0.4 % RMS with a
  2 s correlation time, plus a slow periodic hunt) and cycle-to-cycle
  combustion variability, which widen each order from a line into a band;
- accelerometer vibration rectification, specified in mg per g^2;
- gyroscope angular vibration via an effective lever arm, and the gyroscope's
  own linear-acceleration sensitivity.

Amplitudes are physical rather than fitted.  A line of order `k` is an
inertial excitation growing as `f^2`, transmitted through a mount of natural
frequency `f_n` and damping ratio `zeta` with

`T(f) = sqrt((1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2))`, `r = f/f_n`,

so one overall gain, calibrated once so the hull broadband RMS equals
`level_mps2` at `reference_rpm`, fixes the level at every other speed.  The
speed dependence that falls out is the familiar one rather than a plain `rpm^2`
law: a diesel auxiliary is rough at idle, where the low orders sit near the
mount resonance and `T > 1`.  At the default settings the hull RMS runs from
0.381 m/s^2 at 800 rpm to 0.808 m/s^2 at 3200 rpm, a factor of 2.1 over a
factor of 4 in speed.

The sensor is given a finite two-pole anti-alias bandwidth ahead of the 200 Hz
sample rate.  Each line's amplitude is attenuated at its **true** frequency and
its phase is then advanced at the sample rate, so orders above Nyquist fold on
their own; at 2400 rpm the 4.5 and 6.0 orders (180 and 240 Hz) are recorded at
20 and 40 Hz, and at 2000 rpm the 6.0 order lands on the sample rate exactly
and folds to DC.

The model is **off unless `W3D_ENGINE_RPM` is set**, so every existing noise
realization in the repository is bit-identical.  It is installed in the shared
`process_wave_file_for_tracker` runner, so all three families see the same
vibration.  Configuration:

| Variable | Default | Meaning |
| --- | ---: | --- |
| `W3D_ENGINE_RPM` | off | engine speed; 0 or unset disables the model |
| `W3D_ENGINE_CYLINDERS` | 3 | four-stroke cylinder count; firing order is `n/2` |
| `W3D_ENGINE_GEAR_RATIO` | 2.6 | engine rev per shaft rev |
| `W3D_ENGINE_BLADES` | 3 | propeller blade count |
| `W3D_ENGINE_LEVEL_MPS2` | 0.60 | hull broadband RMS at the reference speed |
| `W3D_ENGINE_REFERENCE_RPM` | 2400 | speed at which the level is stated |
| `W3D_ENGINE_MOUNT_HZ` | 10 | mount natural frequency |
| `W3D_ENGINE_MOUNT_ZETA` | 0.12 | mount damping ratio |
| `W3D_ENGINE_BANDWIDTH_HZ` | 80 | sensor anti-alias bandwidth |
| `W3D_ENGINE_GYRO_LEVER_M` | 1.5 | lever arm for hull angular vibration |
| `W3D_ENGINE_GYRO_G_SENS` | 1.78e-4 | gyro g-sensitivity, (rad/s)/(m/s^2) |
| `W3D_ENGINE_VRE_MG_PER_G2` | 1.0 | accelerometer vibration rectification |
| `W3D_ENGINE_STOP_SEC` | 0 | shut the engine down this far into the record; 0 runs it throughout |
| `W3D_ENGINE_SEED` | 20260828 | phase and modulation seed |

The simulator prints an `ENGINE_VIBRATION` banner with the resulting hull and
recorded RMS, followed by one `ENGINE_LINE` per harmonic giving its true
frequency, its aliased frequency, and its recorded amplitude.

## Current results and reproduction

All records use oceanography-waves-lib v1.2.1 vessel-rao-28ft. The 456-run
degradation matrix, 192-run OU-III mitigation matrix, and paired cross-family
guard check use the final 900 s at the selected filter settings, with unchanged quality gates.

Run `tools/engine_noise_degradation.py`, `tools/ou3_engine_noise_mitigation.py`
and `tools/engine_noise_guard_families.py` with `--data-dir plots/kalman_ou_ii`.
Current raw measurements, summaries, manifests and figures are in
`reports/results/engine_noise_degradation/` and
`reports/results/engine_noise_mitigation/`.

The generated reports retain every engine condition and the matched cross-family guard check. Engine-off parity and the engine-on tradeoffs are reported explicitly.

The fixed deployed guard uses two low-pass poles at 14 Hz and a 25 Hz detector,
with 0.03–0.08 m/s² engagement rails. The covariance gain is 0.75. These runs
do not optimize those settings or isolate the cost of quiet-input group delay.
The preset is not a measured 28 ft engine installation, and the study does
not include encounter-frequency changes, propulsion-induced CG motion or
magnetic interference from machinery.
