# Calibration validation

The reference tree is main at `d2d34116e3a1a1ca7e26cde49e0f2473dcb14a08`.
This report separates synthetic numerical validation, firmware resource
checks, and physical qualification. The procedure and reuse assessment are
in [calibration-architecture.md](calibration-architecture.md).

## Numerical checks

The original calibration suite was built and run before the extensions.
`AccelCalFit.h` and `AccelCalCapture.h` have no changes against the reference
tree. All 17 printed accelerometer campaign/gravity summary rows reproduce
the baseline exactly, including the 30-seed robust-fit campaign. The old
schema assumptions in fixtures were updated without replacing the fitter.

| Check | Result |
| --- | --- |
| Retained workflow and persistence fault tests | 676/676 and 8995/8995 checks |
| Accelerometer campaign | 1627/1627 checks, 30 seeds; ten release-reference scenarios pass |
| Analytic gyro SO(3) sensitivity | Matches independent finite differences |
| Full gyro matrix, exact synthetic rotations | Frobenius matrix error 9.61e-8; independent endpoint RMS 4.09e-8 rad |
| Eight noisy gyro seeds | Worst matrix error 0.0010231 versus 0.0613840 for bias-only identity, about 60-fold smaller |
| Gyro rejection tests | Degenerate axes, corrupted independent rotations, gaps and NaNs fail |
| Final accelerometer model at gyro endpoints | Recovers gyro matrix with distorted/bias-shifted raw accelerometer endpoints |
| Magnetic alignment with wizard-sized pose sets | 16 fit/5 independent poses; independent dip RMS 1.02e-5 degrees in the exact fixture |
| Frame/schema tests | Noncommuting composition, gravity-only yaw ambiguity, reflected rotations, old keys/versions, truncation and CRC corruption |
| Forty-minute white-noise trace | 480001 samples, 2400000 ms; accel N=0.00200198 versus 0.002, gyro N=0.000199806 versus 0.0002 |
| Long-term identifiability | White-only trace reports no bias-instability or walk terms; short traces cannot report a model |
| Mixed flicker-like/random-walk trace | Accel plateau coefficient 0.00237667 versus approximately 0.00225904; gyro walk 3.12645e-5 versus 3e-5 |
| Correlated magnetic Gaussian innovations | 12 gross rejects in 100000 trials at the 0.01% gross-tail threshold |
| Ninety-degree magnetic heading disturbance | All six Kalman implementations reject 100 disturbances with zero heading change in the fixture; clean recovery passes |
| Covariance and injection | Rejected magnetic updates preserve covariance; accepted Joseph updates remain finite, symmetric and PSD |
| Actual qMEKF compass consumer | Dip/norm and 90-degree disturbances cannot redefine the trusted reference; recovery and late first acquisition pass |
| Retained gyro-bias projection regression | OU-II/III still project large accepted corrections; gross magnetic outliers preserve the complete mean and covariance |
| Wizard task dispatch | Actual dispatcher runs off the caller stack, retains the caller's core, publishes save/cancel correctly, handles allocation failure, and services a registered caller watchdog |

The mixed-noise test uses independent stochastic sources: logarithmically
spaced stationary OU components approximate a flicker plateau over a finite
band, while integrated white increments create known gyro random walk. This
does not claim that a plateau identifies flicker outside the measured range.
All coefficients above use the units and Allan convention documented in the
architecture guide. These are synthetic results, not sensor specifications.

Reproduce the focused checks with GCC 13.3 / C++20 and Eigen:

```sh
make -C tests/imu_calibrate build
bash tests/imu_calibrate/run_tests.sh
```

The replay executable uses the exact same bounded-memory analyzer as the
firmware, rather than a second Allan implementation:

```sh
tests/imu_calibrate/sensor_noise-replay capture.log > noise.json
```

It reads calibrated six-axis NOISERAW records with wrapping uint32
microsecond timestamps. Firmware also checks independent magnetic stillness;
that check is not reconstructed by this six-axis replay.

## Firmware resources and reproduction

Toolchain: ESP32 Arduino 3.3.7, M5Unified 0.2.13, M5GFX 0.2.19,
Arduino Eigen 0.3.2. Enable OPI PSRAM and the default_8MB partition scheme:

```sh
arduino-cli compile --library . --jobs 1 --warnings default \
  --build-property 'compiler.cpp.extra_flags=-funroll-loops -fno-finite-math-only -fstack-usage' \
  --fqbn esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc,PSRAM=opi,PartitionScheme=default_8MB \
  sensors/full_marine_ins/atomS3R_ins_kalman_ou3
```

The 680-byte schema retains existing NVS access and verified writes. Target
DWARF type layouts give 3104 bytes for the streaming Allan state (3120 bytes
on the host ABI), 656600 bytes for the 32768-sample rotation fitter including
bounded metadata, and 84824 bytes for the wizard. Rotation samples occupy
20 bytes each; the fitter is explicitly allocated in PSRAM and freed after
fitting. The wizard is heap-owned behind the existing runner entry, rather than
placed on the loop-task stack or retained globally by each compass.
The UI task and fitting worker each have a 32768-byte stack, so a fit can use
65536 bytes of task stacks in addition to caller/RTOS memory. The compiler
reports 64 bytes for the public dispatcher and 48 bytes for its runner. The
4144-byte workflow and 720-byte restore wrapper execute on the UI task;
neither is nested beneath the application's setup/loop stack. ESP-IDF stack sizes and high-water
marks are interpreted as bytes.

Crossbuilds use the documented board options and CI optimization flags:

| Firmware | Flash bytes | Global internal RAM bytes |
| --- | ---: | ---: |
| OU-III INS | 1348359 | 31168 |
| qMEKF compass | 1062811 | 28056 |
| Synthetic benchmark | 695083 | 25624 |

Global RAM excludes runtime heap and task-stack allocations. The benchmark's
compiler-reported individual stack frames include gyro fit 3392 bytes,
normal-equation evaluation 944 bytes, the NIS solve lambda 368 bytes, and
the benchmark entry 6896 bytes. Production worker entry frames are 2512 bytes
for the gyro/alignment fit and 4880 bytes for stationary noise acquisition;
the wizard task entry itself uses 32 bytes. These are individual static frame estimates;
they are neither measured task peaks nor a proof of total stack headroom.
Optimization, inlining and call chains affect the actual peak.

Build and run `sensors/imu_basic/atomS3R_calibration_benchmark` on an AtomS3R
to collect BENCH records: gyro fit microseconds/accuracy, mean and worst
noise-update time, three magnetic-gate timing paths, free internal RAM,
free PSRAM, allocation sizes, and actual task stack headroom. The production
wizard also logs fit time and stack headroom. No board was available in the
implementation environment; no on-device runtime values are claimed.

## Qualification limits

A full native repository run uses `make all` and the pinned simulation
archive from oceanography-waves-lib v1.2.3, SHA256
`87e3cee9505439ec3c0a3ea8c7516789a8f48e3412314c1dcd483894908d433c`.
The focused checks pass through both Make and CMake. The full native replay
run is still in progress; its result will be recorded before the PR is ready.
The separate 1135-test Python run currently reports 18 failures and 8 errors,
all in committed proof/replay provenance and reproduction checks. The CMake
source-path contract and the magnetic effective-covariance source binding pass.
The independent proof check currently fails closed:

```sh
python3 -m tools.stability.ou3_theorem.build_evidence --output /tmp/ocean-stability-evidence.json
```

It exits 1 with `validation_pass=false`: three changed source bindings, an
exact-certificate reproduction mismatch, and nine source-bound diagnostic
failures. The changed magnetic operator requires new qualification; old
fingerprints and finite diagnostic payloads have not been relabeled.
See [the proof research state](ou3-proof-research-state.md) for the retained
structures and the next qualification step. This remains a release blocker.

Build serially on constrained hosts: parallel Eigen-heavy native and MCU
compilers exhausted this environment's 8 GB memory limit. A serial build
completed. The full OU-III image exceeds the default 1310720-byte application
slot; use the documented 3342336-byte OTA slots on the board's 8 MB flash.

Real calibration accuracy requires independently captured device rotations,
fixture registration, and magnetic-environment checks. Real sensor noise
requires a fresh uninterrupted 40-minute stationary acquisition after warmup.
The noise fit does not qualify a marine-motion model or arbitrarily long
physical bias histories. The revised magnetic service policy requires fresh
stability qualification; archived stability evidence is not promoted by
these regression tests.
