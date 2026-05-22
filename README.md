# ocean-imu

Marine IMU and wave-processing algorithms in modern C++ for sensor fusion, sea-state estimation, and simulation-driven validation.

## Motivation

A marine AHRS cannot just reuse typical  popular drone and aerospace IMU filters unchanged. In aerospace, motion is usually modeled as rotation about the center of mass (satellite), and drones often initialize while sitting still before takeoff, so the accelerometer gives a clean gravity direction. On a ship, the system may be turned on while already moving in waves and wind, with heave, roll, pitch, and translational accelerations all mixed into the IMU signals. That means the filter has to learn tilt during motion, avoid trusting wave-distorted acceleration as pure gravity, and keep working across very different sea conditions. In practice, a ship AHRS/INS needs wave-aware initialization, motion compensation, and tuning that can adapt to different sea states dynamically.

The algorithms presented here not only implement tilt-compensated compass and basic roll/pitch/rate-of-turn sensors, they provide corrections for wave induced motion,
and additionally reconstruct 3D displacement of a vessel (heave - strongly observable, surge/sway - weakly) in real-time and estimate apparent (to the vessel being observer) waves direction. 

The code of filters is written in C++ and can run on a microcontroller such as esp32 or on a regular computer. The compass calibration implemented as a part of this library can be run directly on your microcontroller unit. The testing framework uses Stokes/Airy waves with Pierson-Moskowitz or JONSWAP spectrums and a choice of a directional spread model (cosine by default).

## Methods

The articles describing the math behind the methods: 

INS Filters:

- [3D Wave Kalman with OU](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/kalman_ou-w3d.pdf)

- [PII Observer](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/pii_observer-model.pdf)


Frequency tracking:

- [KalmANF Frequency Tracker](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/freq-tracking_adaptive_notch_kalman.pdf)

- [Aranovskiy Frequency Tracker](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/freq-tracking_aranovskiy.pdf)

- [PLL Frequency Tracker](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/freq-tracking_pll.pdf)

- [Zero Crossing Frequency Tracker](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/freq-tracking_zero_crossing.pdf)


IMU Calibration:

- [IMU Calibration Method](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/imu_calibrate-method.pdf)


Wave Models:

- [Wave Models](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/wave_sim-waves.pdf)

- [Fenton Waves](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/wave_sim-fenton.pdf)

- [Spectral Models](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/wave_sim-spectral.pdf)

- [Sea Metrics](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/spectrum-sea_merics.pdf)

- [Vessel RAO](https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/wave_sim-vessel-RAO.pdf)


There are two versions of Kalman INS filters and one filter (PII observer) based on control theory:

- OU_III is using higher-order integral drift correction. It is 21 dimentional state Kalman filter. It seems works better on high precision IMUs.
- OU_II is using more direct integral drift correction and more responsive to sea state changes. It is 18 dimentional state Kalman filter.
- PII observer is based on control theory. It is very computationally light-weight with no matrix operations. It's less accurate than Klaman filters.
- All above filters are adaptive.
- All filters tested to run on esp32s3.
- All filters tested to run on Windows and Linux as well.

Arduino .ino schetches for esp32s3 (on atomS3R):

- Kalman OU_II: https://github.com/bareboat-necessities/ocean-imu/tree/main/sensors/full_marine_ins/atomS3R_ins_kalman_ou2
- Kalman OU_III: https://github.com/bareboat-necessities/ocean-imu/tree/main/sensors/full_marine_ins/atomS3R_ins_kalman_ou3
- PII observer: https://github.com/bareboat-necessities/ocean-imu/tree/main/sensors/full_marine_ins/atomS3R_ins_pii_ovserver


## Overview

`ocean-imu` collects reusable components for marine motion estimation and wave analytics. The codebase is organized so individual modules can be built and exercised from `tests/*` without requiring a monolithic build system.

### Main capabilities

- Attitude and heading workflows (yaw/roll/pitch, AHRS-oriented utilities)
- Wave and frequency-domain analysis
- Wave direction detection
- Kalman-based sea-state models (OU-II and OU-III variants)
- Tuning and filtering helpers for marine sensor pipelines

## Results

Results of the tests, simulations and documentation of the code build from development branch: https://github.com/bareboat-necessities/ocean-imu/releases/tag/vTest

Kalman3D_Wave_OU_II results

<p align="center">
  <img src="./img/samples/w3d_ou2_pmstokes_medium_zkin.svg?raw=true" style="max-width: 50%;">
</p>


<p align="center">
  <img src="./img/samples/w3d_ou2_pmstokes_medium.svg?raw=true" style="max-width: 50%;">
</p>

## Repository layout

```text
src/                 core algorithms and reusable components
  ahrs/              attitude and heading routines
  avg/               averaging and smoothing helpers
  detrend/           detrending helpers
  discrete/          discrete-time utilities
  freq/              frequency-domain utilities
  imu_calibrate/     IMU calibration logic
  kalman_ou_ii/      OU-II Kalman model components
  kalman_ou_iii/     OU-III Kalman model components
  nmea/              NMEA parsing/helpers
  pii_observer/      observer/filter components
  spectrum/          spectral charts
  tuner/             auto-tuning helpers
  util/              shared support code
  wave_dir/          wave direction estimation

tests/               module-level build and validation targets
  ahrs/              AHRS-focused tests and examples
  detrend/           detrending tests
  freq/              builds freq-track
  imu_calibrate/     IMU calibration tests
  kalman_ou_ii/      builds kalman_ou_ii-sim
  kalman_ou_iii/     builds kalman_ou_iii-sim
  pii_observer/      builds pii_observer-adaptive
  wave_sim/          wave simulation programs

sensors/             sensor integration and application examples
  */                 standalone sensor-oriented demos/utilities

doc/                 module documentation and notes
plots/               generated plotting scripts/assets
img/                 images and sample result figures
```

## Prerequisites

- `g++` with C++20 support
- `make`
- Eigen headers (typically from `libeigen3-dev`)

Typical Eigen include location on Linux:

- `/usr/include/eigen3`

## Simulation data dependency

Some validation and simulation workflows depend on data released in:

- https://github.com/bareboat-necessities/oceanography-waves-lib

In CI/docs, this is often referenced as `sim-data-files.zip` from that project’s releases.

You can fetch and unpack this data for local runs with:

```bash
make fetch-sim-data
```

## Build

Run builds from the specific test folder you want to validate:

```bash
cd tests/freq && make all
cd tests/kalman_ou_ii && make all
cd tests/kalman_ou_iii && make all
cd tests/pii_observer && make all
```

If Eigen is not on the default include path for your environment, pass an include override:

```bash
make all CPPFLAGS+='-I/usr/include/eigen3'
```

## Validation

Primary project validation command (when available in your environment):

```bash
make all
```

For module-level validation, run `make all` inside the relevant folder under `tests/`.

