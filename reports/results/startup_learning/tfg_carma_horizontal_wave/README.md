# Transverse wave correction in the TFG-CARMA filters

Paired native replays compare the PR before the correction (`00f1829351418f13768d8c5c01245244714407e2`) with the shared implementation in `4b6d7fc999521713018358ff454b350ae07ea960`. Each replay covers a 1,200 s record at 200 Hz. Scores use its final 900 s. Both filters use the same wave-direction and gyro inputs in each pair.

## Irregular head and following seas

Six draws on 16 JONSWAP and PM–Stokes records give 96 pairs per filter. Each table entry is the square root of the mean of the squared record RMS values.

| Filter | 3D displacement (m) | Y displacement (m) | 3D accelerometer bias (m/s²) | Yaw (°) |
|---|---:|---:|---:|---:|
| Independent-axis | 0.22223 → 0.21884 | 0.13696 → 0.13362 | 0.03667 → 0.03628 | 0.6818 → 0.6804 |
| Shared-mode | 0.20187 → 0.19852 | 0.13065 → 0.12719 | 0.03639 → 0.03600 | 0.6782 → 0.6769 |

| Filter | Following-sea 3D RMS change | Head-sea 3D RMS change |
|---|---:|---:|
| Independent-axis | -1.410% | -1.642% |
| Shared-mode | -1.757% | -1.580% |

| Filter | X accelerometer-bias RMS (m/s²) | Y accelerometer-bias RMS (m/s²) | Z accelerometer-bias RMS (m/s²) |
|---|---:|---:|---:|
| Independent-axis | 0.025290 → 0.025105 | 0.026485 → 0.026121 | 0.001911 → 0.001911 |
| Shared-mode | 0.024924 → 0.024737 | 0.026441 → 0.026076 | 0.001939 → 0.001938 |

On the fixed draw all 64 filter, gyro-profile, direction and wave-record combinations passed their pre-existing gates. On the 8.5 m PM–Stokes head sea with default gyro timing, Y RMS changed as follows.

| Filter | Prior (m) | With correction (m) | Prior ceiling (m) |
|---|---:|---:|---:|
| Independent-axis | 0.29204 | 0.28868 | 0.29351 |
| Shared-mode | 0.28691 | 0.28343 | 0.28835 |

The effect develops during live tracking. The same 96 pairs per filter gave the following displacement RMS by time after the first sample.

| Filter | 40–60 s (m) | 60–150 s (m) |
|---|---:|---:|
| Independent-axis | 1.30863 → 1.30863 | 2.60699 → 2.60694 |
| Shared-mode | 1.33509 → 1.33509 | 2.85139 → 2.85135 |

## Held-out noise and motion

Changes below use the same paired calculation. Seeds 203 and 509 are also shown separately from the six-draw average. After freezing the release source, we selected two further draws, seeds 1063 and 3559, on all 16 irregular records. Regular waves use seed 203 on 24 cnoidal, Gerstner and Fenton records per filter. The vibration runs use seed 203 on the four 8.5 m irregular records at each of 0.60 and 1.20 m/s², 2,400 RPM and 80 Hz bandwidth.

| Filter | Set | Pairs | 3D displacement | Accelerometer bias | Yaw |
|---|---|---:|---:|---:|---:|
| Independent-axis | Seeds 203 and 509 | 32 | -1.602% | -1.133% | +0.105% |
| Independent-axis | Post-freeze draws | 32 | -1.462% | -0.832% | -0.139% |
| Independent-axis | Regular waves | 24 | -0.932% | -0.327% | +0.114% |
| Independent-axis | Engine vibration | 8 | -0.496% | -0.806% | -0.081% |
| Shared-mode | Seeds 203 and 509 | 32 | -1.733% | -1.210% | +0.088% |
| Shared-mode | Post-freeze draws | 32 | -1.560% | -0.804% | -0.147% |
| Shared-mode | Regular waves | 24 | -1.780% | -0.333% | +0.087% |
| Shared-mode | Engine vibration | 8 | -0.451% | -1.147% | -0.099% |

The simulator changes initialization and sensor-noise draws with `W3D_SEED`. Draws on the same wave trace therefore test sensor and startup variability. The regular-wave records and vibration cases test different inputs. These measurements come from simulation. [Per-record scores and input hashes](paired-rms.csv)

## Reproduction

Use the v1.4.0 inputs locked in `tools/sim_dataset_manifest.json`, x86-64-v3, GCC 14.2, Eigen 3.4 and the native Ubuntu math library from the strict replay workflow. At each commit, build `tests/kalman_tfg_carma/tfg_carma3-sim` and `tests/kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim` with their Makefiles. Run each executable separately for each input.

```sh
OCEAN_IMU_SIM_CASE=head W3D_VALIDATION_WINDOW_SEC=900 W3D_WRITE_TIMESERIES=0 \
  tests/kalman_tfg_carma/tfg_carma3-sim --input "$wave_csv"
```

Repeat for `head` and `following`. Leave `W3D_SEED` unset for the default draw, then use 7, 41, 101, 203 and 509 individually. The additional held-out draws use 1063 and 3559. To emit the two earlier windows, set `W3D_VALIDATION_SEGMENTS=early:40:60,release:60:150`. Pass `--regression` for the additional fixed-draw legacy gyro profile. Set `W3D_ALL_WAVE_TYPES=1` for regular waves. For vibration set `W3D_ENGINE_RPM=2400`, `W3D_ENGINE_BANDWIDTH_HZ=80` and `W3D_ENGINE_LEVEL_MPS2` to the stated level. [Filter equations](../../../../docs/tfg-carma3-research.md#horizontal-wave-correction)

The native simulator prints per-axis accelerometer-bias RMS separately from the aggregate validation metrics. The CSV retains their reported precision. [Source, binaries and score hash](receipt.json)

## Fixed-draw limits

All 544 cells and 2,768 gated channels passed the earlier limits before the ordinary 0.5% rounding margin was reapplied. Of the 576 per-record thresholds for the two TFG-CARMA filters, 404 tightened, 135 increased and 37 held. The largest increase was 0.202%.

| Default global limit | Independent-axis, before → after | Shared-mode, before → after |
|---|---:|---:|
| Y displacement (m) | 0.29351 → 0.29013 | 0.28835 → 0.28486 |
| 3D displacement, JONSWAP (% reference maximum) | 6.65939 → 6.65676 | 6.26656 → 6.26406 |
| 3D displacement, PM–Stokes (% reference maximum) | 7.21953 → 7.22067 | 6.66428 → 6.66433 |
| 3D accelerometer bias (%) | 57.12069 → 56.72944 | 57.51143 → 57.01154 |
| Pitch (°) | 0.18233 → 0.18075 | 0.17526 → 0.17400 |
| Yaw (°) | 0.94453 → 0.94456 | 0.93490 → 0.93488 |

Outside the two TFG-CARMA filters, 66 numeric fields in the full matrix differed from the previous published replay. All were in QMEKF shared-mode CARMA. The largest relative difference was 0.00333%, from 27.0459% to 27.0450% on legacy Y accelerometer bias in the 1.5 m following PM–Stokes record.

[Current global limits](../../../../tests/common/CurrentReplayLimits.h) · [Per-record limits and source receipts](../../filter_sim_gates/per_record_limits.json)
