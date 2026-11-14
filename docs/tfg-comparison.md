# OU-III and TFG on the moving-boat records

The full companion study uses oceanography-waves-lib v1.3.1: a 28 ft analytical hull response at 5 kn, with ±3° heading over 120 s. It retains 180 replay rows: nine conditions, ten paired wave/noise/initialization seed triplets, and two deployed estimators. Each history lasts 1200 s at 200 Hz; the final 900 s are scored. Every complete finite replay remains in the statistics, including deterministic regression-gate misses.

The JONSWAP aggregate averages the four heights and the controlled crossfade within each seed. These are ten paired seed differences, not fifty independent conditions.

| Ensemble | OU-III heave RMS (% Hs) | TFG | TFG minus OU-III, points [bootstrap 95% interval] |
| --- | ---: | ---: | ---: |
| JONSWAP plus crossfade | 6.312 | 6.087 | -0.226 [-0.266, -0.186] |
| PM-Stokes | 5.989 | 5.831 | -0.157 [-0.185, -0.127] |
| Controlled crossfade | 6.067 | 5.807 | -0.260 [-0.368, -0.138] |

TFG has lower mean heave error in these declared ensembles. The other channels show tradeoffs:

| Metric, JONSWAP plus crossfade | OU-III | TFG |
| --- | ---: | ---: |
| 3D displacement RMS (m) | 0.41489 | 0.40521 |
| Roll RMS (deg) | 0.24237 | 0.29599 |
| Pitch RMS (deg) | 0.25258 | 0.26661 |
| Yaw RMS (deg) | 1.9343 | 1.8554 |
| Accelerometer-bias RMS (m/s²) | 0.060425 | 0.071343 |
| Gyroscope-bias RMS (rad/s) | 0.0001286 | 0.00010923 |

These secondary channels are exploratory without multiplicity adjustment. The comparison concerns complete estimators.

Reproduce with `python3 tools/tfg_comparison.py --mode full`. Raw rows, paired summaries, gate diagnostics, input/binary hashes and producer workflow receipts are in `reports/results/tfg_comparison/`. Publication requires `--check --require-ci`; source changes, including changes to the cruise-track helper, invalidate reuse.

The [stationary versus moving comparison](moving-vessel-comparison.md) separately measures the effect of changing physical input with the same executable.
