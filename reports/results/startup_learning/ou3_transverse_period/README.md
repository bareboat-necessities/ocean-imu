# Testing the measured transverse period in OU-III

OU-III has an explicit integral of displacement, $S$. The TFG–CARMA filters reconstruct this quantity from their wave states. This experiment tests whether the same cross-wave observation improves the deployed OU-III filter on 1,200 s head- and following-sea replays sampled at 200 Hz.

The measured-period candidate lowered six-draw 3D displacement RMS by 6.70% and passed 10 of 32 fixed replay checks. Individual pitch, yaw and accelerometer-bias increases kept the shipping OU-III drift observation in place.

## Observation

The shipping drift observation targets $S=0$. The candidate replaces its cross-wave component with $n^T(S+\beta v)=0$ and retains the original along-wave and vertical components. The unit vector $n$ comes from the measured wave axis and estimated boat heading. The candidate estimates $\beta$ from centered transverse displacement and velocity variances over eight encounter periods. It starts after four periods of history and requires live bias learning and a reliable axis. The Jacobian includes the effect of attitude on the wave axis in NED coordinates. The principal candidate uses the same 0.625 transverse residual-variance fraction as TFG-CARMA.

## Late tracking errors

The main score covers seconds 300–1,200. Values pool the squared per-record RMS errors before taking their square root. Six draws cover the default noise sequence and seeds 7, 41, 101, 203 and 509 across all 16 irregular-sea records. The reserved draws use seeds 1063 and 3559. Engine tests use seed 203 on four 8.5 m records at each of two vibration levels.

| Measured transverse period | Pairs | 3D displacement (m) | Y displacement (m) | 3D accelerometer bias (m/s²) | Pitch (°) | Yaw (°) |
|---|---:|---:|---:|---:|---:|---:|
| Six noise draws | 96 | 0.42728 → 0.39865 | 0.24322 → 0.22278 | 0.03834 → 0.03822 | 0.16559 → 0.16428 | 0.73747 → 0.74737 |
| Reserved noise draws | 32 | 0.43288 → 0.40317 | 0.25033 → 0.22931 | 0.05738 → 0.05802 | 0.19164 → 0.19856 | 0.76604 → 0.77268 |
| Engine vibration | 8 | 0.87929 → 0.82700 | 0.52958 → 0.50375 | 0.13799 → 0.13993 | 0.75115 → 0.76174 | 4.89836 → 4.96197 |

### Record limits and attitude feedback

The along-wave and vertical variant withholds direct transverse attitude and gyro-bias updates. The final two variants keep the original drift innovation for the named inertial states. Each state consumes one innovation. Their Joseph covariance update includes the cross-covariance between the two residuals.

| Default gyro profile, 16 records | Passed | 3D displacement change | Accel-bias change | Pitch change | Yaw change |
|---|---:|---:|---:|---:|---:|
| Disabled, same compiled experiment | 16/16 | +0.00% | +0.00% | +0.00% | +0.00% |
| Measured transverse period | 5/16 | -7.03% | +2.35% | +3.88% | +3.21% |
| Measured period, original variance | 4/16 | -4.67% | +3.23% | +5.02% | +3.14% |
| Vertical period on the transverse axis | 6/16 | -6.14% | +0.24% | +0.12% | +1.28% |
| Along-wave and vertical attitude updates | 0/16 | +18.95% | +15.39% | +43.81% | +11.38% |
| Original drift observation for attitude, biases and acceleration | 0/16 | +12.46% | +9.59% | +21.80% | +4.81% |
| Original drift observation for attitude and biases | 0/16 | +12.20% | +9.59% | +21.44% | +4.68% |

The principal candidate passed 5/16 default and 5/16 legacy-gyro checks. Its largest default-draw pitch increase was 0.06910° to 0.09252° on the following-sea record `wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv`. The largest default-draw yaw increase was 0.80380° to 0.91668°. On reserved noise draws the largest pitch increase was 0.11967° to 0.15135°.

Both variants that retained the original inertial-state drift observation failed all 16 default replay checks. The disabled variant reproduced all 624 comparison fields across the three windows.

## Reproduction

Start from `e663f0e6b7cf4df9f23727ed80bcf8d81404d3e2` and the v1.4.0 input archives. Apply `01-transverse.patch` and build `tests/kalman_ou_iii/kalman_ou_iii-sim` with the Makefile and the recorded GCC, Eigen and math-library versions. Use this build for arms 0–4. Apply `02-correlated-block.patch` on top for arms 5–6. The experiment selects its observation with `OU3_TRANSVERSE_ARM`.

```sh
OCEAN_IMU_SIM_CASE=head OU3_TRANSVERSE_ARM=1 W3D_WRITE_TIMESERIES=0 \
  W3D_VALIDATION_WINDOW_SEC=900 W3D_VALIDATION_SEGMENTS=early:40:60,release:60:150 \
  tests/kalman_ou_iii/kalman_ou_iii-sim --input "$wave_csv"
```

Repeat with following-sea inputs. Set `W3D_SEED` for additional sensor draws. The fixed legacy profile adds `--regression` and uses default sensor noise. For vibration, set `W3D_ENGINE_RPM=2400`, `W3D_ENGINE_BANDWIDTH_HZ=80`, and `W3D_ENGINE_LEVEL_MPS2` to 0.60 or 1.20. Fixed replay limits apply to the default noise draw. The other draws measure paired physical errors.

[Per-record RMS for all three windows](paired-rms.csv.gz) · [Checksums and per-sea measurements](receipt.json) · [Exact patches and original logs](reproduction.tar.gz)
