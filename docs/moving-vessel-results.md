# Underway replay results

The current external-data baseline uses oceanography-waves-lib v1.4.0: a 28 ft estimated hull RAO at 5 kn, following seas northbound and head/bow-quartering seas southbound, both with a ±3° heading sine every 20 s. Each case contains eight 1200 s, 200 Hz records. The final 900 s are scored, except AHRS attitude (final 60 s). Displacement and velocity truth are relative to the declared cruise track; total acceleration and delivered IMU samples are retained. These are deterministic regression measurements, not confidence intervals or sea trials.

## Northbound, following seas

| Filter | JONSWAP heave (% Hs) | PM-Stokes heave (% Hs) | 3D RMS (m) | Yaw RMS (deg) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 8.551 | 8.089 | 1.057 | 0.910 |
| OU-III | 8.405 | 7.973 | 0.879 | 0.804 |
| TFG | 5.868 | 5.563 | 0.855 | 0.915 |
| TFG cubic | 4.794 | 4.909 | 1.430 | 3.808 |
| CARMA | 3.541 | 3.650 | 0.416 | 0.904 |
| Shared QMEKF | 3.445 | 3.479 | 0.350 | 0.902 |
| TFG-CARMA | 3.541 | 3.650 | 0.416 | 0.894 |
| Shared TFG | 3.445 | 3.479 | 0.350 | 0.888 |
| Shared tangent | 3.445 | 3.479 | 0.350 | 0.887 |
| AHRS | — | — | — | 1.426 |
| Adaptive PII | 10.145 | 9.284 | — | 2.422 |
| TVG-NLO | 9.583 | 9.145 | — | — |
| Heave PII | 10.194 | 9.352 | — | — |
| Godhavn | 382.948 | 296.548 | — | — |
| Godhavn extended | 17.407 | 18.196 | — | — |
| Richter | 16.237 | 16.918 | — | — |
| Kuchler | 13.419 | 14.422 | — | — |

[Full metrics](../reports/results/sim_v1.4.0/following/metrics.csv) · [Charts and provenance](../reports/results/sim_v1.4.0/following/README.md)

## Southbound, head / bow-quartering seas

| Filter | JONSWAP heave (% Hs) | PM-Stokes heave (% Hs) | 3D RMS (m) | Yaw RMS (deg) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 6.076 | 5.862 | 0.874 | 1.038 |
| OU-III | 3.621 | 3.558 | 0.694 | 1.088 |
| TFG | 2.773 | 2.782 | 0.684 | 1.047 |
| TFG cubic | 2.705 | 2.682 | 0.990 | 1.264 |
| CARMA | 2.166 | 2.533 | 0.446 | 0.945 |
| Shared QMEKF | 2.097 | 2.362 | 0.425 | 0.935 |
| TFG-CARMA | 2.166 | 2.533 | 0.446 | 0.940 |
| Shared TFG | 2.097 | 2.363 | 0.425 | 0.930 |
| Shared tangent | 2.097 | 2.362 | 0.425 | 0.930 |
| AHRS | — | — | — | 1.654 |
| Adaptive PII | 10.327 | 13.747 | — | 3.011 |
| TVG-NLO | 7.171 | 7.824 | — | — |
| Heave PII | 10.483 | 13.855 | — | — |
| Godhavn | 395.715 | 489.737 | — | — |
| Godhavn extended | 26.438 | 31.342 | — | — |
| Richter | 22.413 | 25.517 | — | — |
| Kuchler | 20.910 | 21.765 | — | — |

[Full metrics](../reports/results/sim_v1.4.0/head/metrics.csv) · [Charts and provenance](../reports/results/sim_v1.4.0/head/README.md)

Each column is the worst of eight records; maxima need not occur together. Each case/sea/profile has its own fixed limit. All profiles, roll/pitch, bias metrics and input/executable hashes remain in the case-specific JSON/CSV reports.

See [the scoring convention](moving-vessel-replays.md) and [gate protocol](filter-simulator-gates.md).
