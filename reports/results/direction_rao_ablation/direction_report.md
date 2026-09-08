# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Axis and directed angular errors use available estimates only; axis availability and unresolved travel fractions retain the full scoring window as their denominator. An unavailable axis is not a zero-degree bearing. A dash means at least one record has no finite angular score, so the eight-record mean is undefined; missing records are not silently excluded.

| Family | RPM | Arm | Axis error ° | Axis available % | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 11.62 | 76.10 | 19.60 | 71.76 | 0.42 | 27.82 |
| OU-II | 0 | guard-only | 11.62 | 76.10 | 19.60 | 71.76 | 0.42 | 27.82 |
| OU-II | 0 | guard-and-rao | 1.78 | 73.03 | 16.11 | 69.59 | 0.21 | 30.20 |
| OU-II | 2400 | unconditioned | 45.69 | 78.89 | 78.35 | 5.35 | 1.39 | 93.26 |
| OU-II | 2400 | guard-only | 10.80 | 75.07 | 23.25 | 64.20 | 0.83 | 34.97 |
| OU-II | 2400 | guard-and-rao | 17.19 | 72.61 | -- | 68.25 | 0.22 | 31.53 |
| OU-III | 0 | unconditioned | 11.79 | 76.00 | 19.75 | 71.78 | 0.43 | 27.79 |
| OU-III | 0 | guard-only | 11.79 | 76.00 | 19.75 | 71.78 | 0.43 | 27.79 |
| OU-III | 0 | guard-and-rao | 2.23 | 72.96 | 15.91 | 69.50 | 0.21 | 30.29 |
| OU-III | 2400 | unconditioned | 47.35 | 77.60 | 91.86 | 5.12 | 2.58 | 92.30 |
| OU-III | 2400 | guard-only | 11.03 | 75.11 | 23.35 | 64.31 | 0.79 | 34.90 |
| OU-III | 2400 | guard-and-rao | 18.74 | 72.55 | -- | 68.08 | 0.24 | 31.69 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. Angular accuracy is conditional on usable motion. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
