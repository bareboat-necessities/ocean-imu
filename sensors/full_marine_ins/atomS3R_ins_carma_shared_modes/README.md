# AtomS3R marine INS — shared-modes QMEKF (experimental)

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_carma_shared_modes**

This sketch runs the **experimental** standalone shared-modes QMEKF, [`KalmanCarmaSharedModes`](../../../src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h), on the M5Stack AtomS3R. It is [atomS3R_ins_carma3](../atomS3R_ins_carma3/README.md) with one change: the three independent CARMA(4,1) motion priors are replaced by **shared directional stochastic modes**. Each axis carries one complexified CARMA(4,1) mode in its six-state minimal realization, the axes share its dynamics, and a Hermitian excitation is estimated from the filter's own world-frame acceleration. The model can therefore represent the horizontal-vertical quadrature of waves. A **coupling confidence** withdraws the coupling and falls back toward the independent-axis prior when the coupling cannot be trusted: a prior inconsistent with recent data, little wave energy, manoeuvring, abrupt sea-state changes or inconsistent innovations ([design](../../../docs/kalman-carma-shared-modes.md)).

Open the sketch: [`atomS3R_ins_carma_shared_modes.ino`](atomS3R_ins_carma_shared_modes.ino).

> **Research sketch, not the shipping estimator, and not run on hardware yet.** The evidence is simulation-only: the validation in [`docs/kalman-carma-shared-modes.md`](../../../docs/kalman-carma-shared-modes.md) and [`doc/kalman_carma_shared_modes`](../../../doc/kalman_carma_shared_modes). There is no sea trial, no stability analysis of the adaptive filter, and no measured stack high-water mark. Use [OU-III](../atomS3R_ins_kalman_ou3/README.md) for deployment.

## Differences from atomS3R_ins_carma3

- One line: the filter type is `SeaStateFusion_CarmaSharedModes`, the CARMA3 wrapper with [`KalmanCarmaSharedModes<float>`](../../../src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h) as its inner filter (27 error states instead of 21). With the coupling off it is the CARMA3 filter's law exactly. Startup, magnetic north lock, tuner, wave-excitation evidence, outputs, sentences, display and calibration paths are the CARMA3 sketch's.
- `begin(cfg)` selects the deployable configuration: the coupling is estimated from the filter's own acceleration, with the coupling confidence. No oracle, RAO or truth information is used.

## Memory and timing (ESP32-S3, Xtensa GCC 14.2.0, `-Os`)

| | CARMA3 (A41) | shared modes |
|---|---|---|
| QMEKF object (heap, one allocation at `begin`) | 6.1 kB | 12.5 kB |
| wrapper object (heap) | 3.2 kB | 3.2 kB |
| largest single stack frame | 1.3 kB | 1.9 kB |
| worst static call chain of one update | 4.8 kB | 4.9 kB |

These figures come from `tests/kalman_carma_shared_modes/arduino/csm_size_probe/build.sh`, a static analysis. Calls through function pointers and library frames are not counted, so the call-chain figure is a lower bound; measure the loop task's stack high-water mark on the device before relying on the 8 KB default loop stack. Host update time of the filter core is about 1.6 times that of CARMA3; the embedded per-update time has not been measured.
