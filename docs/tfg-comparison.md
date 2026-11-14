# OU-III and TFG on the moving-boat records

The full companion study uses oceanography-waves-lib v1.3.0: a 28 ft analytical hull response at 5 kn, with ±3° heading over 120 s. It retains 180 replay rows: nine conditions, ten paired wave/noise/initialization seed triplets, and two deployed estimators. Each record lasts 1200 s at 200 Hz; the final 900 s are scored. All complete finite replays remain in the statistics, including deterministic regression-gate misses.

The JONSWAP companion aggregate averages the four heights **and the controlled crossfade** within each seed. It therefore differs from the four-height primary OU-II/OU-III endpoint. These are ten paired seed differences, not fifty independent conditions.

| Ensemble | OU-III heave RMS (% Hs) | TFG | TFG minus OU-III, points [bootstrap 95% interval] |
| --- | ---: | ---: | ---: |
| JONSWAP plus crossfade | 6.31 | 6.09 | −0.225 [−0.266, −0.185] |
| PM-Stokes | 5.99 | 5.83 | −0.16 [−0.18, −0.13] |
| Controlled crossfade | 6.07 | 5.81 | −0.26 [−0.37, −0.14] |

TFG has lower mean heave error in these declared ensembles. The other channels show tradeoffs: in the JONSWAP-plus-crossfade aggregate, accelerometer-bias RMS rises from 0.06042 m/s² (OU-III) to 0.07157 m/s² (TFG), roll RMS from 0.2424° to 0.2963°, and pitch from 0.2525° to 0.2687°. TFG has lower mean 3D displacement, yaw, and gyro-bias errors. These secondary channels are exploratory without multiplicity adjustment; the comparison concerns complete estimators.

Reproduce with `python3 tools/tfg_comparison.py --mode full`. Raw rows, paired summaries, gate diagnostics, input/binary hashes, and the publication table are in `reports/results/tfg_comparison/`. The current local full replay has no GitHub workflow receipt. The CI publication gate still requires its own replay provenance; it is not bypassed by the local results.

The deterministic [stationary versus moving comparison](moving-vessel-comparison.md) answers a different question: what changes when the identical executable sees the new physical input.
