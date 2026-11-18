# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Axis and directed angular errors use available estimates only; axis availability and unresolved travel fractions retain the full scoring window as their denominator. An unavailable axis is not a zero-degree bearing. A dash means at least one record has no finite angular score, so the eight-record mean is undefined; missing records are not silently excluded.

| Family | RPM | Arm | Axis error ° | Axis available % | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 12.50 | 78.70 | 19.49 | 76.28 | 0.46 | 23.27 |
| OU-II | 0 | guard-only | 12.50 | 78.70 | 19.49 | 76.28 | 0.46 | 23.27 |
| OU-II | 0 | guard-and-rao | 7.84 | 73.79 | 17.67 | 70.34 | 0.23 | 29.43 |
| OU-II | 2400 | unconditioned | 37.50 | 98.05 | -- | 32.87 | 2.26 | 64.87 |
| OU-II | 2400 | guard-only | 12.53 | 79.01 | 19.38 | 76.71 | 0.35 | 22.93 |
| OU-II | 2400 | guard-and-rao | 9.11 | 73.72 | 18.34 | 70.24 | 0.14 | 29.62 |
| OU-III | 0 | unconditioned | 12.43 | 78.63 | 19.12 | 76.26 | 0.40 | 23.33 |
| OU-III | 0 | guard-only | 12.43 | 78.63 | 19.12 | 76.26 | 0.40 | 23.33 |
| OU-III | 0 | guard-and-rao | 7.71 | 73.77 | 17.50 | 70.43 | 0.22 | 29.35 |
| OU-III | 2400 | unconditioned | 32.20 | 84.25 | 58.65 | 31.69 | 3.16 | 65.15 |
| OU-III | 2400 | guard-only | 12.34 | 78.84 | 19.74 | 76.50 | 0.32 | 23.18 |
| OU-III | 2400 | guard-and-rao | 8.69 | 73.66 | 18.36 | 70.18 | 0.17 | 29.65 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. Angular accuracy is conditional on usable motion. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
