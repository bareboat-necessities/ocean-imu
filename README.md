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

## Results

Results of the tests, simulations and documentation of the code build from development branch: https://github.com/bareboat-necessities/ocean-imu/releases/tag/vTest

Kalman3D_Wave_OU_III results

<p align="center">
  <img src="./img/samples/w3d_ou3_pmstokes_medium_zkin.svg?raw=true" style="max-width: 50%;">
</p>

<p align="center">
  <img src="./img/samples/w3d_ou3_pmstokes_medium.svg?raw=true" style="max-width: 50%;">
</p>

## Main Article

This article describes the main method: https://github.com/bareboat-necessities/ocean-imu/releases/download/vTest/kalman_ou-w3d.pdf

## Repository layout

```text
src/
  ahrs/             AHRS-related components
  freq/             frequency-domain utilities
  kalman_ou_ii/     OU-II Kalman model components
  kalman_ou_iii/    OU-III Kalman model components
  pii_observer/     observer/filter components
  tuner/            auto-tuning helpers
  util/             shared support code (e.g., simulation helpers)
  wave_dir/         wave direction estimation

tests/
  freq/             builds freq-track
  kalman_ou_ii/     builds kalman_ou_ii-sim
  kalman_ou_iii/    builds kalman_ou_iii-sim
  pii_observer/     builds pii_observer-adaptive
  spectrum/         placeholder target
  wave_sim/         placeholder target
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
