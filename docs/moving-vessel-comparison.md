# Stationary versus moving-boat measurements

Before is a fresh replay of the v1.2.3 stationary-vessel records; after is
the v1.3.0 5 kn dataset with heading oscillating ±3° every 120 s. Each native
filter uses the **identical rebuilt executable**, default sensor/initialization
draws and trailing 900 s windows in both arms. Input and executable hashes
are recorded in `stationary_control.json` and `measurements.json`. The
change is physical input, with no estimator retuning. The fresh stationary
controls also reproduce the archived regression baseline to the displayed precision.

| Filter | Worst heave before → after (% Hs) | Relative increase | Worst absolute 3D RMS before → after (m) | Accel bias RMS before → after (mm/s²) |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 6.63 → 20.86 | +215% | 0.922 → 1.070 | 52.04 → 131.96 |
| OU-III | 4.25 → 8.24 | +94% | 0.736 → 0.900 | 58.01 → 63.50 |
| TFG | 4.38 → 7.82 | +78% | 0.745 → 0.911 | 57.38 → 156.83 |
| CARMA | 3.21 → 4.67 | +45% | 0.394 → 0.449 | 68.62 → 60.53 |
| Shared QMEKF | 2.45 → 3.70 | +51% | 0.366 → 0.364 | 69.91 → 47.48 |
| TFG-CARMA | 3.21 → 4.67 | +45% | 0.394 → 0.449 | 68.69 → 60.69 |
| Shared TFG | 2.45 → 3.70 | +51% | 0.366 → 0.364 | 69.41 → 47.60 |

The worst normalized heave error is in a **0.27 m sea** for every filter:
JONSWAP for OU-II, OU-III and TFG; PM-Stokes for the CARMA families.
Shared-mode filters retain essentially unchanged worst absolute 3D RMS,
while OU-II and TFG show the largest accelerometer-bias deterioration.

Against the recorded auxiliary baselines, NLO worst raw heave rises from 7.32% to 23.52% Hs. AHRS worst roll/pitch/yaw
improves from 0.93°/5.33°/12.66° to 0.64°/3.00°/7.08°.

A raw-data spectral check on the 0.27 m JONSWAP record finds a heave PSD peak
period of 2.98 s before and 5.96 s after; vertical-acceleration RMS falls
from 0.246 to 0.106 m/s², while true heave RMS changes only 0.0534 to
0.0488 m. This supports the interpretation that lower encounter frequency
and reduced inertial excitation make the low-sea problem harder. It does
not isolate the effects of steering, adaptation and bias observability.

Position/velocity reference truth is expressed about the declared cruise
track. All specific force, turn acceleration and angular rates remain intact.
See [the coordinate convention](moving-vessel-replays.md) and
[current results](moving-vessel-results.md). Columns report their own worst
record and need not peak in the same sea.
