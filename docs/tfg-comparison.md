# OU-III and TFG on the moving-boat records

The full companion study uses oceanography-waves-lib v1.3.1: a 28 ft analytical hull response at 5 kn, with ±3° heading over 120 s. It retains 180 replay rows: nine conditions, ten paired wave/noise/initialization seed triplets, and two deployed estimators. Each history lasts 1200 s at 200 Hz; the final 900 s are scored. Every complete finite replay remains in the statistics, including deterministic regression-gate misses.

The JONSWAP aggregate averages the four heights and the controlled crossfade within each seed. These are ten paired seed differences, not fifty independent conditions. TFG uses the wave-matched vertical integral observation and the spectral-law constrained attitude/gyro gain.

| Ensemble | OU-III heave RMS (% Hs) | TFG | TFG minus OU-III, points [bootstrap 95% interval] |
| --- | ---: | ---: | ---: |
| JONSWAP plus crossfade | 6.312 | 4.515 | -1.798 [-1.848, -1.753] |
| PM-Stokes | 5.989 | 4.315 | -1.674 [-1.712, -1.638] |
| Controlled crossfade | 6.067 | 4.149 | -1.918 [-2.026, -1.811] |

TFG has lower mean heave error in these declared ensembles. The other channels show tradeoffs:

| Metric, JONSWAP plus crossfade | OU-III | TFG |
| --- | ---: | ---: |
| 3D displacement RMS (m) | 0.4149 | 0.38331 |
| Roll RMS (deg) | 0.24237 | 0.29956 |
| Pitch RMS (deg) | 0.25258 | 0.29237 |
| Yaw RMS (deg) | 1.9343 | 1.8758 |
| Accelerometer-bias RMS (m/s²) | 0.060425 | 0.07488 |
| Gyroscope-bias RMS (rad/s) | 0.0001286 | 0.00010928 |

These secondary channels are exploratory without multiplicity adjustment. The comparison concerns complete estimators; it does not isolate group geometry from motion regularization. The crossfade is a controlled transition, not a physically evolving random sea. This is a development RAO ensemble, not qualification on an unseen vessel.

Reproduce with `python3 tools/tfg_comparison.py --mode full`. Raw rows, paired summaries, gate diagnostics, input/binary hashes and producer workflow receipts are in `reports/results/tfg_comparison/`. The current bundle was produced by [CI run 43661880261](https://github.com/bareboat-necessities/ocean-imu/actions/runs/43661880261) at source `6552e14`. Publication requires `--check --require-ci`; source changes, including changes to the cruise-track helper, invalidate reuse.

The [stationary versus moving comparison](moving-vessel-comparison.md) separately measures the effect of changing physical input with the same executable.
