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
| Magnetic alignment with wizard-sized pose sets | 9 fit/4 independent poses; the exact fixture passes independent dip validation |
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
| Turn guide | Start/finish holds are distinct; TWIST and TIP show short actions; hand support is allowed and each target is written once |
| Sensor/UI scheduling | Production worker and UI execute on separate host threads; injected 80 ms LCD stalls preserve 200 Hz capture, the accepted turn, and the 20 ms gap gate |
| Screen layout | Real renderer text fits 128 pixels; USB-up and both edge-down rotations match the accelerometer guide; the real UI header preserves pose orientation through next/retry/confirmation menus |
| Guided geometry | Direct face tips expose an information deficiency; alternating twists followed by the named tips recover the full matrix with independent validation, with unchanged conditioning, physical and residual limits |
| Short gyro holds | A 3.5-second target uses 2 seconds of useful samples and two 1.5-second magnetic windows; seventeen guided campaigns retain matrix error below 0.005 and independent validation |
| Controller replay | Seventeen twelve-turn campaigns cover waits, noise, pauses, gaps, hold disturbances, slow motion, late pose jitter, retries, nonidentity scale/cross-axis matrices and eight noisy nonidentity seeds; matrix Frobenius error is below 0.005 |
| USB-up regressions | A late single-reading pose excursion reproduces the old reversion and now retains the finish screen; a rejected mean in turn 2 opens retry while retaining turn 1; a live tap-to-retry also preserves accepted data |
| Reading and final-fit recovery | Instruction cards wait five simulated minutes for input; a failed twelve-turn fit appends four turns without refitting bias or repeating cards; explicit abort after twenty failed turns preserves the captured set until exit |
| Reduced-plan integrity | Eight training rotations and four untouched validation rotations; magnetic registration from actual held endpoints; corrupted validation rotations and undersized sealed records are rejected |
| Turn retries | The actual wrapper retries a timeout at turn 7 while retaining six accepted turns and fitting bias only once; worker allocation failure also preserves accepted data |

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

The final `timeout 180 make all` rerun reached the unchanged CARMA build and
returned exit 124 while compiling `accel_cal-test.o` (recipe reported
`Terminated`). The complete calibration runner and qMEKF crossbuild pass;
this does not establish a repository-wide pass.

## Firmware resources and reproduction

Toolchain: ESP32 Arduino 3.3.7, M5Unified 0.2.13, M5GFX 0.2.19,
Arduino Eigen 0.3.2. Enable OPI PSRAM and the default_8MB partition scheme:

```sh
arduino-cli compile --library . --jobs 1 --warnings default \
  --build-property 'compiler.cpp.extra_flags=-funroll-loops -fno-finite-math-only -fstack-usage' \
  --fqbn esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc,PSRAM=opi,PartitionScheme=default_8MB \
  sensors/full_marine_ins/atomS3R_ins_kalman_ou3
```

The 680-byte schema retains existing NVS access and verified writes. The streaming Allan state uses 3104 bytes on the target (3120 bytes on
the host ABI). The rotation fitter reserves 131072 twenty-byte samples plus
bounded metadata, about 2.5 MiB of PSRAM. The heap-owned wizard also holds the small atomic capture
status context. Rotation samples occupy
20 bytes each; the fitter is explicitly allocated in PSRAM and freed after
fitting. The wizard is heap-owned behind the existing runner entry, rather than
placed on the loop-task stack or retained globally by each compass.
The UI task and each sensor/fit worker have a 32768-byte stack. Sensor capture
and fitting run sequentially, so either operation can use
65536 bytes of task stacks in addition to caller/RTOS memory. The compiler estimates 1536 bytes for the capture worker frame and 160 bytes
for its UI runner; these are not hardware peak measurements. It
reports 64 bytes for the public dispatcher and 48 bytes for its runner. The
4144-byte workflow and 720-byte restore wrapper execute on the UI task;
neither is nested beneath the application's setup/loop stack. ESP-IDF stack sizes and high-water
marks are interpreted as bytes.

The local ESP32 toolchain was removed by workspace maintenance before the
pose-orientation, finish-hold and twelve-turn corrections. Those changes are verified by
host controller/worker tests and by compiling the full UI header against LCD
and clock stubs; a new device crossbuild is not claimed here. The table below
records the last crossbuild before that correction, using the documented
board options and CI optimization flags:

| Firmware | Flash bytes | Global internal RAM bytes |
| --- | ---: | ---: |
| OU-III INS | 1348359 | 31168 |
| qMEKF compass (before pose UI correction) | 1068699 | 28056 |
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
