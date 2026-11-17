# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Axis and directed angular errors use available estimates only; axis availability and unresolved travel fractions retain the full scoring window as their denominator. An unavailable axis is not a zero-degree bearing. A dash means at least one record has no finite angular score, so the eight-record mean is undefined; missing records are not silently excluded.

| Family | RPM | Arm | Axis error ° | Axis available % | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 9.83 | 72.57 | 56.01 | 60.20 | 1.38 | 38.42 |
| OU-II | 0 | guard-only | 9.83 | 72.57 | 56.01 | 60.20 | 1.38 | 38.42 |
| OU-II | 0 | guard-and-rao | 11.29 | 71.49 | -- | 64.13 | 0.23 | 35.64 |
| OU-II | 2400 | unconditioned | 57.83 | 99.30 | 98.70 | 0.39 | 0.30 | 99.31 |
| OU-II | 2400 | guard-only | 11.52 | 71.43 | 67.69 | 30.70 | 2.61 | 66.69 |
| OU-II | 2400 | guard-and-rao | 12.92 | 70.39 | -- | 61.81 | 0.29 | 37.90 |
| OU-III | 0 | unconditioned | 9.73 | 72.43 | 58.47 | 60.35 | 1.38 | 38.27 |
| OU-III | 0 | guard-only | 9.73 | 72.43 | 58.47 | 60.35 | 1.38 | 38.27 |
| OU-III | 0 | guard-and-rao | 11.25 | 71.28 | -- | 63.86 | 0.23 | 35.91 |
| OU-III | 2400 | unconditioned | 46.30 | 75.66 | 94.91 | 1.61 | 1.69 | 96.70 |
| OU-III | 2400 | guard-only | 11.44 | 71.25 | 59.38 | 31.26 | 2.66 | 66.09 |
| OU-III | 2400 | guard-and-rao | 13.13 | 70.17 | -- | 61.36 | 0.34 | 38.31 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. Angular accuracy is conditional on usable motion. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
