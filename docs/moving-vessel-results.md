# Underway replay results

The current external-data baseline uses oceanography-waves-lib v1.3.1: a 28 ft estimated hull RAO moving at 5 kn with heading oscillating plus/minus 3 degrees about north every 120 s. Eight fixed JONSWAP and PM-Stokes records run for 1200 s at 200 Hz; the last 900 s are scored with default sensor and initialization seeds. Displacement and velocity truth are relative to the declared cruise track. All recorded specific force, turning acceleration, changing attitude and body rates reach the estimator unchanged. The 3D denominator uses wave displacement, not forward travel. These deterministic regression measurements are not confidence intervals or sea trials. A dataset change is not an estimator improvement.

| Filter | JONSWAP heave (% Hs) | PM-Stokes heave (% Hs) | 3D RMS (m) | Yaw RMS (deg) | Accel bias RMS (mm/s²) |
| --- | ---: | ---: | ---: | ---: | ---: |
| OU-II | 20.841 | 11.315 | 1.0704 | 0.893 | 131.91 |
| OU-III | 8.236 | 7.839 | 0.9001 | 0.871 | 63.50 |
| TFG | 7.819 | 7.476 | 0.9107 | 1.633 | 156.77 |
| CARMA | 4.192 | 4.667 | 0.4488 | 0.937 | 60.53 |
| Shared QMEKF | 3.550 | 3.702 | 0.3644 | 0.986 | 47.48 |
| TFG-CARMA | 4.192 | 4.669 | 0.4489 | 0.935 | 60.69 |
| Shared TFG | 3.551 | 3.702 | 0.3645 | 0.952 | 47.60 |
| TFG cubic | 8.247 | 8.016 | 1.5789 | 3.049 | 210.45 |
| Adaptive PII | 15.387 | 14.733 | — | 8.233 | — |

Each column is the worst of eight records; maxima need not occur together.
The compact per-record measurements and both sampling profiles are in `reports/results/filter_sim_gates/measurements.json` and `underway_summary.csv`.
See [the scoring and surrogate convention](moving-vessel-replays.md).
