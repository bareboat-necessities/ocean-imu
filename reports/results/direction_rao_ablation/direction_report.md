# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Axis and directed angular errors use available estimates only; axis availability and unresolved travel fractions retain the full scoring window as their denominator. An unavailable axis is not a zero-degree bearing. A dash means at least one record has no finite angular score, so the eight-record mean is undefined; missing records are not silently excluded.

| Family | RPM | Arm | Axis error ° | Axis available % | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 11.70 | 76.16 | 19.55 | 71.81 | 0.41 | 27.77 |
| OU-II | 0 | guard-only | 11.70 | 76.16 | 19.55 | 71.81 | 0.41 | 27.77 |
| OU-II | 0 | guard-and-rao | 2.36 | 73.02 | 16.35 | 69.60 | 0.23 | 30.17 |
| OU-II | 2400 | unconditioned | 52.46 | 85.09 | 91.44 | 6.63 | 4.37 | 89.00 |
| OU-II | 2400 | guard-only | 11.24 | 75.31 | 25.61 | 64.08 | 0.73 | 35.18 |
| OU-II | 2400 | guard-and-rao | 17.76 | 72.63 | -- | 68.23 | 0.24 | 31.54 |
| OU-III | 0 | unconditioned | 11.74 | 75.98 | 19.64 | 71.79 | 0.42 | 27.79 |
| OU-III | 0 | guard-only | 11.74 | 75.98 | 19.64 | 71.79 | 0.42 | 27.79 |
| OU-III | 0 | guard-and-rao | 2.16 | 72.95 | 15.93 | 69.46 | 0.21 | 30.33 |
| OU-III | 2400 | unconditioned | 47.65 | 77.50 | 88.50 | 5.79 | 2.38 | 91.84 |
| OU-III | 2400 | guard-only | 10.98 | 75.09 | 22.99 | 64.33 | 0.81 | 34.86 |
| OU-III | 2400 | guard-and-rao | 18.70 | 72.55 | -- | 68.13 | 0.23 | 31.65 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. Angular accuracy is conditional on usable motion. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
