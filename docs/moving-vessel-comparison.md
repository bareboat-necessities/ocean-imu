# Stationary versus moving-boat measurements

Both arms are fresh replays of oceanography-waves-lib v1.3.1, with exact 200 Hz clocks. Before is the zero-speed 28 ft vessel archive; after is the 5 kn archive with ±3° heading oscillation every 120 s. Each pair uses the identical rebuilt executable, default sensor and initialization draws, and trailing 900 s scoring. No estimator retuning is applied. Hashes and per-record controls are retained in `stationary_control.json`, `measurements.json`, and `dataset_comparison.json`.

| Filter | Worst heave before → after (% Hs) | Relative increase | Worst 3D RMS before → after (m) | Accel bias RMS before → after (mm/s²) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 6.63 → 20.84 | +214% | 0.922 → 1.070 | 52.04 → 131.91 |
| OU-III | 4.25 → 8.24 | +94% | 0.736 → 0.900 | 58.01 → 63.50 |
| TFG | 4.38 → 7.82 | +78% | 0.745 → 0.911 | 57.36 → 156.77 |
| CARMA | 3.21 → 4.67 | +45% | 0.394 → 0.449 | 68.62 → 60.53 |
| Shared QMEKF | 2.45 → 3.70 | +51% | 0.366 → 0.364 | 69.91 → 47.48 |
| TFG-CARMA | 3.21 → 4.67 | +45% | 0.394 → 0.449 | 68.69 → 60.69 |
| Shared TFG | 2.45 → 3.70 | +51% | 0.366 → 0.364 | 69.41 → 47.60 |
| TFG cubic | 4.32 → 8.25 | +91% | 1.185 → 1.579 | 101.47 → 210.45 |
| Adaptive PII | 8.75 → 15.39 | +76% | — | — |
| NLO | 7.32 → 23.52 | +221% | — | — |

Each column reports its own worst record; these maxima need not occur in the same sea. These deterministic comparisons do not provide population confidence intervals.

| Filter | Roll before → after (deg) | Pitch before → after (deg) | Yaw before → after (deg) |
| --- | ---: | ---: | ---: |
| OU-II | 0.226 → 0.586 | 0.212 → 0.503 | 0.726 → 0.893 |
| OU-III | 0.321 → 0.345 | 0.137 → 0.155 | 0.824 → 0.871 |
| TFG | 0.351 → 0.667 | 0.093 → 0.656 | 0.950 → 1.633 |
| CARMA | 0.349 → 0.310 | 0.202 → 0.195 | 0.761 → 0.937 |
| Shared QMEKF | 0.351 → 0.212 | 0.212 → 0.209 | 0.754 → 0.986 |
| TFG-CARMA | 0.349 → 0.309 | 0.203 → 0.193 | 0.770 → 0.935 |
| Shared TFG | 0.350 → 0.212 | 0.208 → 0.201 | 0.758 → 0.952 |
| TFG cubic | 0.433 → 0.401 | 0.493 → 1.160 | 1.057 → 3.049 |
| Adaptive PII | 5.394 → 4.397 | 4.519 → 3.581 | 10.303 → 8.233 |
| AHRS (last 60 s) | 0.933 → 0.643 | 5.332 → 2.996 | 12.658 → 7.080 |

OU-II and NLO have the largest relative heave deterioration. TFG has a substantial accelerometer-bias and attitude penalty; the shared-mode CARMA variants retain nearly unchanged worst absolute 3D RMS despite worse normalized heave in the smallest sea. These paired inputs change encounter frequency and steering together, so this experiment does not identify their separate causal contributions.

The [engine-noise study](../reports/results/engine_noise_degradation/engine_noise_report.md) separately tests machinery stress on the moving records. At 800 rpm, OU-III with the guard plus covariance inflation reaches 2.417 times its engine-off 3D RMS; the old zero-speed 1.7-times claim does not carry over. Its pitch offset also worsens when covariance inflation is added. This is a measured operating limitation, not a change to deployed tuning.

Position and velocity truth are expressed about the declared cruise track. All specific force, turning acceleration and angular rates remain intact. See the [coordinate convention](moving-vessel-replays.md) and [current results](moving-vessel-results.md).
