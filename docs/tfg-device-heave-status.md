# TFG device-heave investigation status

## Scope and unresolved result

The shipping estimator is unchanged. The reported real-device tens-of-metres
heave error has not been assigned a validated root cause or fixed. The tests
below distinguish the zero-state handoff from later nonlinear transients; they
are not a certificate of device stability. No output clipping, detrender
workaround, coefficient retuning, or evidence-gate relaxation is included.

The investigated source is PR #587 at
`b59ce712476e2a87c2f6aace420474fa1a9b2bf6`, based on #586 / `59118d19`.
Local source identity was checked against GitHub blob identifiers:

- `Kalman3D_Wave_TFG.h`: `d8d9d72c591df15161177164a9012b009b82afeb`.
- `SeaStateFusionFilter_TFG.h`: `66cda600dd624fc91a4b42260fe4e3f5f335db6b`.

## Executed observations

The original 600 s, 200 Hz frozen-bias handoff regression produced:

| Vertical residual (m/s^2) | Final p.z (m) | Peak absolute p.z (m) | Final a_w.z (m/s^2) |
| --- | --- | --- | --- |
| 0 | 0 | 0 | 0 |
| 0.01 | 0.0695634 | 0.0802305 | 0.00950194 |
| 0.03 | 0.208686 | 0.240682 | 0.0285046 |
| 0.05 | 0.347809 | 0.401142 | 0.0475083 |

Accelerometer bias remained zero. The strengthened test additionally checks
all state/covariance values on every sample and prints S.z and its maximum.
Its assertion is 0.5 m, not the former 20 m. This assertion does not change
estimator output.

The full device-wrapper fixture runs 13 ten-minute histories, including
startup, adaptation, magnetic refinement, held/unlocked accelerometer bias,
periodic wave-acceleration covariance synchronization, tilted orientation,
noise/motion, and 10--200 Hz loop rates. The original focused CI run passed
all cases with peak raw displacement below 1.20 m. Its regression assertion
is now 2 m. These are synthetic histories, not recordings from the AtomS3R.

A separate stationary stress replay used body gyro residual
`(0.2, -0.1, 0.06) rad/s`, body accelerometer residual `(0,0,0.05) m/s^2`,
orientation `Rz(0.5) Ry(0.4) Rx(0.2)`, and world field `(20,0,43)`.
The paired replay produced peak raw p.z of 70.121559 m for TFG and
2.675621 m for OU-III, without detrending. The gyro residual is deliberately
large; no available device recording establishes that it represents the
reported device failure. This is a reproducing stress history, not a
calibrated sensor-performance claim.

In a separately scheduled TFG-only causal replay, peak p.z was 72.67 m.
Holding accelerometer bias gave 72.5674 m; disabling periodic covariance sync
gave 43.6439 m; disabling magnetic refinement gave 20.8682 m; supplying the
true gyro bias at handoff gave 66.5424 m. Therefore none of those isolated
changes is a validated fix. Applying adaptation during Cold worsened the
stress peak to 109.35 m and was discarded. Candidate estimator changes were
not retained. Different replay clocks/schedules must not be compared as
bit-identical runs.

## Test/build repairs

The Makefile test runner now actually executes the process-covariance,
axis-factor, handoff, and device regressions. The shared log-frequency test
is included in the detrend runner. CMake builds the missing test targets.
Detrend object dependencies include headers, and its shell runner explicitly
sets `-e` even when invoked using `bash run_tests.sh`.

Device trace CSVs are opt-in through `OCEAN_IMU_TRACE_DIR`, which must name an
existing directory. Normal simulator runs do not receive these diagnostic
CSVs as accidental dataset inputs. The focused workflow retains the traces.
The #584 yaw-free startup magnetic-compass path is unchanged.

## Remaining blockers and evidence boundary

The unresolved software question is the source of the disproportionate
nonlinear TFG correction transient in the stress history, and whether that
same mechanism occurs for the real device's calibrated sample history.
Passing the small-disturbance regression does not answer either question.
A proposed fix must reduce the raw-state failure before detrending, preserve
normal wave performance, and pass the same-history comparison; changing
coefficients merely to suppress this example is not sufficient.

Local `make all` was executed and stopped while compiling
`tests/kalman_ou_ii/kalman_ou_ii-sim.cpp` with:
`g++: fatal error: Killed signal terminated program cc1plus`.
The nested make errors were Makefile:68, Makefile:59, and Makefile:45.
This is a compiler/resource failure, not a passing full validation.

The committed validation and robustness bundles are stale for
`src/detrend/AdaptiveWaveDetrender3D.h`; the evidence contract correctly
rejects both. Full source-bound regeneration and green end-to-end CI remain
required. No refreshed provenance or publication result is claimed here.
