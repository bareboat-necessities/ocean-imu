# Underway replay results

The current external-data baseline uses oceanography-waves-lib v1.3.1: a 28 ft estimated hull RAO moving at 5 kn with heading oscillating plus/minus 3 degrees about north every 120 s. Eight fixed JONSWAP and PM-Stokes records run for 1200 s at 200 Hz; the last 900 s are scored with default sensor and initialization seeds. Displacement and velocity truth are relative to the declared cruise track. All recorded specific force, turning acceleration, changing attitude and body rates reach the estimator unchanged. The 3D denominator uses wave displacement, not forward travel. These deterministic regression measurements are not confidence intervals or sea trials. A dataset change is not an estimator improvement.

| Filter | JONSWAP heave (% Hs) | PM-Stokes heave (% Hs) | 3D RMS (m) | Yaw RMS (deg) | Accel bias RMS (mm/s²) |
| --- | ---: | ---: | ---: | ---: | ---: |
| OU-II | 8.486 | 8.027 | 1.0716 | 0.899 | 45.27 |
| OU-III | 8.236 | 7.839 | 0.9001 | 0.871 | 63.50 |
| TFG | 5.716 | 5.459 | 0.8621 | 1.069 | 91.84 |
| CARMA | 3.572 | 3.726 | 0.4201 | 0.938 | 47.53 |
| Shared QMEKF | 3.456 | 3.532 | 0.3537 | 0.977 | 46.06 |
| TFG-CARMA | 3.572 | 3.727 | 0.4201 | 0.927 | 47.59 |
| Shared TFG | 3.456 | 3.533 | 0.3537 | 0.947 | 45.60 |
| TFG cubic | 4.689 | 4.832 | 1.5367 | 3.063 | 210.33 |
| Adaptive PII | 10.487 | 9.623 | — | 8.233 | — |
| TVG-NLO | 9.430 | 9.119 | — | — | — |

Each column is the worst of eight records; maxima need not occur together.
The compact per-record measurements and both sampling profiles are in `reports/results/filter_sim_gates/measurements.json` and `underway_summary.csv`.
See [the scoring and surrogate convention](moving-vessel-replays.md).
