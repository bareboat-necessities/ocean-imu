# ocean-imu

Marine IMU and wave-processing algorithms in modern C++ for sensor fusion, sea-state estimation, and simulation-driven validation.

## Overview

`ocean-imu` collects reusable components for marine motion estimation and wave analytics. The codebase is organized so individual modules can be built and exercised from `tests/*` without requiring a monolithic build system.

### Main capabilities

- Attitude and heading workflows (yaw/roll/pitch, AHRS-oriented utilities)
- Wave and frequency-domain analysis
- Wave direction detection
- Kalman-based sea-state models (OU-II and OU-III variants)
- Tuning and filtering helpers for marine sensor pipelines

## Motivation

A marine AHRS cannot just reuse typical  popular drone and aerospace IMU filters unchanged. In aerospace, motion is usually modeled as rotation about the center of mass (satellite), and drones often initialize while sitting still before takeoff, so the accelerometer gives a clean gravity direction. On a ship, the system may be turned on while already moving in waves and wind, with heave, roll, pitch, and translational accelerations all mixed into the IMU signals. That means the filter has to learn tilt during motion, avoid trusting wave-distorted acceleration as pure gravity, and keep working across very different sea states. In practice, a ship AHRS/INS needs wave-aware initialization, motion compensation, and tuning that can adapt to different sea states dynamically.

## Results

Results of the tests, simulations and documentation of the code build from development branch: https://github.com/bareboat-necessities/ocean-imu/releases/tag/vTest

Kalman3D_Wave_OU_II results

<p align="center">
  <img src="./img/samples/w3d_ou2_pmstokes_medium_zkin.svg?raw=true" style="max-width: 50%;">
</p>


<p align="center">
  <img src="./img/samples/w3d_ou2_pmstokes_medium.svg?raw=true" style="max-width: 50%;">
</p>

## Main Article

This article describes the main method: https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/kalman_ou-w3d.pdf

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

## Notes

- This repository is designed for targeted, module-by-module builds.
- Keep performance-sensitive behavior deterministic and validated via simulation-oriented checks.
