# Calibration validation

The reference tree is main at `d2d34116e3a1a1ca7e26cde49e0f2473dcb14a08`.
This report separates synthetic numerical validation, firmware resource
checks, and physical qualification. The procedure and reuse assessment are
in [calibration-architecture.md](calibration-architecture.md); requirement status is
in [calibration-requirements.md](calibration-requirements.md).

## Numerical checks

The original calibration suite was built and run before the extensions.
`AccelCalFit.h`, `AccelCalCapture.h`, `MagCalFit.h` and `MagCalCapture.h`
are byte-for-byte unchanged against the reference tree. All 17 printed accelerometer campaign/gravity summary rows reproduce
the baseline exactly, including the 30-seed robust-fit campaign. The old
schema assumptions in fixtures were updated without replacing the fitter.

| Check | Result |
| --- | --- |
| Retained workflow and persistence fault tests | 676/676 and 8995/8995 checks |
| Accelerometer campaign | 1650/1650 checks, 30 seeds; ten release-reference scenarios pass |
| Analytic gyro SO(3) sensitivity | Matches independent finite differences |
| Full gyro matrix, exact synthetic rotations | Frobenius matrix error 9.95e-8; independent endpoint RMS 5.81e-8 rad |
| Eight noisy gyro seeds | Worst matrix error 0.00104775 versus 0.0613840 for bias-only identity, about 59-fold smaller |
| Gyro rejection tests | Degenerate axes, corrupted independent rotations, gaps and NaNs fail |
| Final accelerometer model at gyro endpoints | Recovers gyro matrix with distorted/bias-shifted raw accelerometer endpoints |
| Magnetic alignment with wizard-sized pose sets | 9 fit/4 independent poses; the exact fixture passes independent dip validation |
| Fresh full-calibration admission | Actual entry selects no previous slope or gyro reference; calibrated noise mode and cancellation remain separate |
| Noise worker | Complete simulated 40-minute sensor acquisition, wrapping clock, stopped/missing/frozen samples, temperature and cancellation; no invalid run yields a model |
| Noise runtime wiring | Unified record round-trip retains density units and identified-axis masks; short and invalid records cannot configure a filter |
| Frame/schema tests | Noncommuting composition, gravity-only yaw ambiguity, reflected rotations, old keys/versions, truncation and CRC corruption |
| Forty-minute white-noise trace | 480001 samples, 2400000 ms; accel N=0.00200198 versus 0.002, gyro N=0.000199806 versus 0.0002 |
| Long-term identifiability | White-only trace reports no bias-instability or walk terms; short traces cannot report a model |
| Mixed flicker-like/random-walk trace | Accel plateau coefficient 0.00237667 versus approximately 0.00225904; gyro walk 3.12645e-5 versus 3e-5 |
| Correlated magnetic Gaussian innovations | 12 gross rejects in 100000 trials at the 0.01% gross-tail threshold |
| Ninety-degree magnetic heading disturbance | All six Kalman implementations reject 100 disturbances with zero heading change in the fixture; clean recovery passes |
| Gross field strength with uncertain covariance | All six cores reject 1.6x/0.6x field magnitude even when large R makes NIS small; the test reproduces TFG's previously omitted ratio guard |
| Covariance and injection | Rejected magnetic updates preserve covariance; accepted Joseph updates remain finite, symmetric and PSD |
| Actual qMEKF compass consumer | Dip/norm and 90-degree disturbances cannot redefine the trusted reference; recovery and late first acquisition pass |
| Retained gyro-bias projection regression | OU-II/III still project large accepted corrections; gross magnetic outliers preserve the complete mean and covariance |
| Wizard task dispatch | Actual dispatcher runs off the caller stack, retains the caller's core, publishes save/cancel correctly, handles allocation failure, and services a registered caller watchdog |
| Turn guide | Start/finish holds are distinct; TWIST and TIP show short actions; hand support is allowed and each target is written once |
| Sensor/UI scheduling | Production worker and UI execute on separate host threads; injected 80 ms LCD stalls preserve 200 Hz capture, the accepted turn, and the 20 ms gap gate |
| Screen layout | Real renderer text fits 128 pixels; USB-up and both edge-down rotations match the accelerometer guide; the real UI header preserves pose orientation through next/retry/confirmation menus |
| Guided geometry | Direct face tips expose an information deficiency; alternating twists followed by the named tips recover the full matrix with independent validation, with unchanged conditioning, physical and residual limits |
| Human-supported endpoints | Bias stillness profile stays strict; a distinct endpoint profile admits 0.1-degree rocking at 4/6/8 Hz with coupled gravity/rate/lever-arm acceleration. All three finish within 4.2 s of placement; worst matrix error 0.002822 versus the 0.005 test bound |
| Slow placement | Twenty seconds of constant slow yaw remains unqualified and waiting in place; no false acceptance or premature retry |
| Short gyro holds | A 3.5-second target uses 2 seconds of useful samples and two 1.5-second magnetic windows; twenty guided campaigns retain matrix error below 0.005 and independent validation |
| Controller replay | Twenty twelve-turn campaigns cover waits, noise, pauses, gaps, hold disturbances, slow motion, late pose jitter, retries, nonidentity scale/cross-axis matrices and eight noisy nonidentity seeds; matrix Frobenius error is below 0.005 |
| USB-up regressions | A late single-reading pose excursion reproduces the old reversion and now retains the finish screen; a misplaced finish pose recovers in place without discarding earlier turns; a live tap-to-retry also preserves accepted data |
| Reading and final-fit recovery | Instruction cards wait five simulated minutes for input; a failed twelve-turn fit appends four turns without refitting bias or repeating cards; explicit abort after twenty failed turns preserves the captured set until exit |
| Reduced-plan integrity | Eight training rotations and four untouched validation rotations; magnetic registration from actual held endpoints; corrupted validation rotations and undersized sealed records are rejected |
| Turn retries | The actual wrapper retries a timeout at turn 7 while retaining six accepted turns and fitting bias only once; worker allocation failure also preserves accepted data |

The mixed-noise test uses independent stochastic sources: logarithmically
spaced stationary OU components approximate a flicker plateau over a finite
band, while integrated white increments create known gyro random walk. This
does not claim that a plateau identifies flicker outside the measured range.
All coefficients above use the units and Allan convention documented in the
architecture guide. These are synthetic results, not sensor specifications.

Reproduce the focused checks with GCC 13.3 / C++20 and Eigen 3.4.0:

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

The focused calibration build and complete runner pass. The new worker/admission regressions and TFG strength-anomaly fixture fail
against the preceding source for the intended behavioral reasons: stopped-stream
handling, previous-coefficient reuse and omission of the common field-ratio guard.

## Firmware resources and reproduction

Toolchain: Arduino CLI 1.3.1, ESP32 Arduino 3.3.7, M5Unified 0.2.13,
M5GFX 0.2.19, Arduino Eigen 0.3.2 and ArxTypeTraits 0.3.2. Enable OPI PSRAM and the default_8MB partition scheme:

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
65536 bytes of task stacks in addition to caller/RTOS memory. The compiler estimates 1552 bytes for the capture worker frame and 128 bytes
for its UI runner; these are not hardware peak measurements. It
reports 64 bytes for the public dispatcher and 48 bytes for its runner. The
4144-byte workflow and 720-byte restore wrapper execute on the UI task;
neither is nested beneath the application's setup/loop stack. ESP-IDF stack sizes and high-water
marks are interpreted as bytes.

The architecture audit crossbuilds the qMEKF compass, OU-III INS, TFG INS and
benchmark with the documented board options and CI optimization flags. The
TFG image was rebuilt again after the final human-supported endpoint changes;
its row below includes the latest shared wizard/capture translation unit.
The other rows record the architecture-audit builds before that last UI/capture
correction. Compilation and linking do not establish operation on a board.

| Firmware | Flash bytes | Global internal RAM bytes |
| --- | ---: | ---: |
| OU-III INS | 1354811 | 31168 |
| qMEKF compass | 1069683 | 28056 |
| TFG INS (latest endpoint profile) | 1274739 | 40904 |
| Synthetic benchmark | 695803 | 25624 |

Global RAM excludes runtime heap and task-stack allocations. Production
worker entry frames are 2512 bytes for the gyro/alignment fit and 4880 bytes
for stationary noise acquisition; the wizard task entry uses 32 bytes.
These are individual static frame estimates, not measured task peaks or a
proof of total stack headroom. Optimization, inlining and call chains affect
the actual peak.

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
All native targets build successfully when compiled serially. The default
`make all` algorithm suites run through `tests/util` and stop in validation:
`python3 ../../tools/tfg_comparison.py --check --require-ci` raises
`ValueError: TFG comparison source provenance is stale` (outer make exit 2).
After the final changes, all directories were rebuilt and the same known
validation blocker checked first, without removing any build/test directory.
The full TFG native suite was rerun after the strength-guard correction and
passes, as do the TFG-CARMA3 core (199 checks), measurement/covariance contracts
(44 checks), and the wave-direction and wave-simulation suites after that
barrier.
This does not establish a repository-wide green result.
The separate 1135-test Python run reports 14 failures, 6 errors and 1 skip,
all in committed proof/replay provenance and reproduction checks. The CMake
source-path contract and the magnetic effective-covariance source binding pass.
The independent proof check currently fails closed:

```sh
python3 -m tools.stability.ou3_theorem.build_evidence --output /tmp/ocean-stability-evidence.json
```

It exits 1 with `validation_pass=false`: seven changed source bindings, two
exact-certificate reproduction mismatches, and nine source-bound diagnostic
failures. The changed magnetic operator requires new qualification; old
fingerprints and finite diagnostic payloads have not been relabeled.
See [the proof research state](ou3-proof-research-state.md) for the retained
structures and the next qualification step. This remains a release blocker.

Run Eigen-heavy native and MCU compilations serially on constrained hosts.
Concurrent single-job builds hit this host's 8 GB limit while compiling
`tfg_carma3-test.cpp`: `g++: fatal error: Killed signal terminated program cc1plus`. The full OU-III image exceeds the default 1310720-byte application
slot; use the documented 3342336-byte OTA slots on the board's 8 MB flash.

Real calibration accuracy requires independently captured device rotations,
fixture registration, and magnetic-environment checks. Real sensor noise
requires a fresh uninterrupted 40-minute stationary acquisition after warmup.
The noise fit does not qualify a marine-motion model or arbitrarily long
physical bias histories. The revised magnetic service policy requires fresh
stability qualification; archived stability evidence is not promoted by
these regression tests.
