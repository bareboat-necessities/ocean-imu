# Stationary versus moving-boat measurements

Both arms are fresh replays of oceanography-waves-lib v1.3.1, with exact 200 Hz clocks. Before is the zero-speed 28 ft vessel archive; after is the 5 kn archive with ±3° heading oscillation every 120 s. Each pair uses the identical rebuilt executable, default sensor and initialization draws, and trailing 900 s scoring. No estimator retuning is applied. Hashes and per-record controls are retained in `stationary_control.json`, `measurements.json`, and `dataset_comparison.json`.

| Filter | Worst heave before → after (% Hs) | Relative increase | Worst 3D RMS before → after (m) | Accel bias RMS before → after (mm/s²) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 6.64 → 8.49 | +28% | 0.922 → 1.072 | 52.05 → 45.27 |
| OU-III | 4.25 → 8.24 | +94% | 0.736 → 0.900 | 58.01 → 63.50 |
| TFG | 3.23 → 5.72 | +77% | 0.712 → 0.862 | 56.79 → 91.84 |
| CARMA | 2.46 → 3.73 | +51% | 0.396 → 0.420 | 68.71 → 47.53 |
| Shared QMEKF | 2.36 → 3.53 | +50% | 0.360 → 0.354 | 68.98 → 46.06 |
| TFG-CARMA | 2.46 → 3.73 | +51% | 0.396 → 0.420 | 68.78 → 47.59 |
| Shared TFG | 2.36 → 3.53 | +50% | 0.360 → 0.354 | 68.97 → 45.60 |
| TFG cubic | 2.97 → 4.83 | +63% | 1.153 → 1.537 | 101.77 → 210.33 |
| Adaptive PII | 8.75 → 10.49 | +20% | — | — |
| NLO | 6.98 → 9.43 | +35% | — | — |

Each column reports its own worst record; these maxima need not occur in the same sea. These deterministic comparisons do not provide population confidence intervals.

| Filter | Roll before → after (deg) | Pitch before → after (deg) | Yaw before → after (deg) |
| --- | ---: | ---: | ---: |
| OU-II | 0.224 → 0.213 | 0.212 → 0.188 | 0.742 → 0.899 |
| OU-III | 0.321 → 0.345 | 0.137 → 0.155 | 0.824 → 0.871 |
| TFG | 0.310 → 0.529 | 0.210 → 0.246 | 0.884 → 1.069 |
| CARMA | 0.349 → 0.211 | 0.203 → 0.194 | 0.763 → 0.938 |
| Shared QMEKF | 0.352 → 0.202 | 0.204 → 0.206 | 0.769 → 0.977 |
| TFG-CARMA | 0.349 → 0.211 | 0.203 → 0.191 | 0.770 → 0.927 |
| Shared TFG | 0.352 → 0.199 | 0.201 → 0.199 | 0.777 → 0.947 |
| TFG cubic | 0.439 → 0.401 | 0.497 → 1.166 | 1.130 → 3.063 |
| Adaptive PII | 5.394 → 4.397 | 4.519 → 3.581 | 10.303 → 8.233 |
| AHRS (last 60 s) | 0.933 → 0.643 | 5.332 → 2.996 | 12.658 → 7.080 |

The largest relative increases in worst heave are, in order: OU-III, TFG, TFG cubic. OU-II and the CARMA wrappers do not let still-water detection attenuate wave variance, which a 0.27 m following sea at 5 kn triggers. TFG has a substantial accelerometer-bias and attitude penalty; the shared-mode CARMA variants retain nearly unchanged worst absolute 3D RMS despite worse normalized heave in the smallest sea. These paired inputs change encounter frequency and steering together. The [regression diagnosis](moving-vessel-regression-diagnosis.md) explains the stillness fix and the limits of attributing the remaining differences.

The [engine-noise study](../reports/results/engine_noise_degradation/engine_noise_report.md) separately tests machinery stress on the moving records. At 800 rpm, OU-III with the guard plus covariance inflation reaches 2.417 times its engine-off 3D RMS; the old zero-speed 1.7-times claim does not carry over. Its pitch offset also worsens when covariance inflation is added. This is a measured operating limitation, not a change to deployed tuning.

Position and velocity truth are expressed about the declared cruise track. All specific force, turning acceleration and angular rates remain intact. See the [coordinate convention](moving-vessel-replays.md) and [current results](moving-vessel-results.md).
