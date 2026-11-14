# Stationary versus moving-boat measurements

Both arms are fresh replays of oceanography-waves-lib v1.3.1, with exact 200 Hz clocks. Before is the zero-speed 28 ft vessel archive; after is the 5 kn archive with ±3° heading oscillation every 120 s. Each pair uses the identical rebuilt executable, default sensor and initialization draws, and trailing 900 s scoring. No estimator retuning is applied. Hashes and per-record controls are retained in `stationary_control.json`, `measurements.json`, and `dataset_comparison.json`.

| Filter | Worst heave before → after (% Hs) | Relative increase | Worst 3D RMS before → after (m) | Accel bias RMS before → after (mm/s²) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 6.65 → 8.48 | +28% | 0.922 → 1.072 | 52.05 → 45.26 |
| OU-III | 4.25 → 8.24 | +94% | 0.736 → 0.900 | 58.01 → 63.50 |
| TFG | 4.38 → 7.82 | +78% | 0.745 → 0.911 | 57.36 → 156.77 |
| CARMA | 2.46 → 3.73 | +51% | 0.396 → 0.420 | 68.71 → 47.53 |
| Shared QMEKF | 2.45 → 3.70 | +51% | 0.366 → 0.364 | 69.91 → 47.48 |
| TFG-CARMA | 2.46 → 3.73 | +51% | 0.396 → 0.420 | 68.78 → 47.59 |
| Shared TFG | 2.45 → 3.70 | +51% | 0.366 → 0.364 | 69.41 → 47.60 |
| TFG cubic | 4.32 → 8.25 | +91% | 1.185 → 1.579 | 101.47 → 210.45 |
| Adaptive PII | 8.75 → 15.39 | +76% | — | — |
| NLO | 7.32 → 23.52 | +221% | — | — |

Each column reports its own worst record; these maxima need not occur in the same sea. These deterministic comparisons do not provide population confidence intervals.

| Filter | Roll before → after (deg) | Pitch before → after (deg) | Yaw before → after (deg) |
| --- | ---: | ---: | ---: |
| OU-II | 0.224 → 0.213 | 0.212 → 0.188 | 0.742 → 0.899 |
| OU-III | 0.321 → 0.345 | 0.137 → 0.155 | 0.824 → 0.871 |
| TFG | 0.351 → 0.667 | 0.093 → 0.656 | 0.950 → 1.633 |
| CARMA | 0.349 → 0.211 | 0.203 → 0.194 | 0.763 → 0.939 |
| Shared QMEKF | 0.351 → 0.212 | 0.212 → 0.209 | 0.754 → 0.986 |
| TFG-CARMA | 0.349 → 0.211 | 0.203 → 0.191 | 0.770 → 0.928 |
| Shared TFG | 0.350 → 0.212 | 0.208 → 0.201 | 0.758 → 0.952 |
| TFG cubic | 0.433 → 0.401 | 0.493 → 1.160 | 1.057 → 3.049 |
| Adaptive PII | 5.394 → 4.397 | 4.519 → 3.581 | 10.303 → 8.233 |
| AHRS (last 60 s) | 0.933 → 0.643 | 5.332 → 2.996 | 12.658 → 7.080 |

NLO has the largest relative heave deterioration, followed by OU-III and cubic TFG. OU-II and the CARMA wrappers do not let still-water detection attenuate wave variance, which a 0.27 m following sea at 5 kn triggers. TFG has a substantial accelerometer-bias and attitude penalty; the shared-mode CARMA variants retain nearly unchanged worst absolute 3D RMS despite worse normalized heave in the smallest sea. These paired inputs change encounter frequency and steering together; the [regression diagnosis](moving-vessel-regression-diagnosis.md) separates them with counterfactual records and attributes the change to the encounter frequency.

<!-- BEGIN current-peak-records -->
The following records maximize each reported error in the current moving arm; these are different endpoints, not a single worst sea.

| Filter | Normalized heave | Accelerometer bias | Yaw |
| --- | --- | --- | --- |
| OU-II | JONSWAP 0.27 m | JONSWAP 0.27 m | JONSWAP 1.5 m |
| OU-III | JONSWAP 0.27 m | PM-Stokes 0.27 m | JONSWAP 4 m |
| TFG | JONSWAP 0.27 m | PM-Stokes 8.5 m | PM-Stokes 8.5 m |
| CARMA | PM-Stokes 0.27 m | PM-Stokes 0.27 m | JONSWAP 1.5 m |
| Shared QMEKF | PM-Stokes 0.27 m | PM-Stokes 0.27 m | JONSWAP 1.5 m |
| TFG-CARMA | PM-Stokes 0.27 m | PM-Stokes 0.27 m | JONSWAP 1.5 m |
| Shared TFG | PM-Stokes 0.27 m | PM-Stokes 0.27 m | JONSWAP 1.5 m |
| TFG cubic | JONSWAP 0.27 m | JONSWAP 8.5 m | JONSWAP 8.5 m |
| Adaptive PII | JONSWAP 1.5 m | — | PM-Stokes 8.5 m |

NLO normalized heave also peaks in 0.27 m PM-Stokes. TFG’s worst bias and yaw errors occur in the 8.5 m PM-Stokes record; its degradation is not confined to small waves.
<!-- END current-peak-records -->

The [engine-noise study](../reports/results/engine_noise_degradation/engine_noise_report.md) separately tests machinery stress on the moving records. At 800 rpm, OU-III with the guard plus covariance inflation reaches 2.417 times its engine-off 3D RMS; the old zero-speed 1.7-times claim does not carry over. Its pitch offset also worsens when covariance inflation is added. This is a measured operating limitation, not a change to deployed tuning.

Position and velocity truth are expressed about the declared cruise track. All specific force, turning acceleration and angular rates remain intact. See the [coordinate convention](moving-vessel-replays.md) and [current results](moving-vessel-results.md).
