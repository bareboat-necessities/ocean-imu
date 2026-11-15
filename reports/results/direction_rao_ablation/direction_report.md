# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Axis and directed angular errors use available estimates only; axis availability and unresolved travel fractions retain the full scoring window as their denominator. An unavailable axis is not a zero-degree bearing. A dash means at least one record has no finite angular score, so the eight-record mean is undefined; missing records are not silently excluded.

| Family | RPM | Arm | Axis error ° | Axis available % | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 4.47 | 71.82 | -- | 61.10 | 1.09 | 37.80 |
| OU-II | 0 | guard-only | 4.47 | 71.82 | -- | 61.10 | 1.09 | 37.80 |
| OU-II | 0 | guard-and-rao | -- | 71.17 | -- | 64.55 | 0.24 | 35.20 |
| OU-II | 2400 | unconditioned | 46.77 | 75.04 | 104.97 | 0.28 | 0.29 | 99.43 |
| OU-II | 2400 | guard-only | -- | 71.25 | -- | 31.30 | 2.09 | 66.61 |
| OU-II | 2400 | guard-and-rao | -- | 70.22 | -- | 62.55 | 0.14 | 37.30 |
| OU-III | 0 | unconditioned | 4.37 | 71.68 | -- | 61.31 | 0.97 | 37.72 |
| OU-III | 0 | guard-only | 4.37 | 71.68 | -- | 61.31 | 0.97 | 37.72 |
| OU-III | 0 | guard-and-rao | -- | 70.99 | -- | 64.67 | 0.22 | 35.11 |
| OU-III | 2400 | unconditioned | 45.72 | 74.43 | 106.70 | 0.77 | 0.79 | 98.44 |
| OU-III | 2400 | guard-only | -- | 71.10 | -- | 31.60 | 2.18 | 66.22 |
| OU-III | 2400 | guard-and-rao | -- | 70.00 | -- | 61.93 | 0.17 | 37.90 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. Angular accuracy is conditional on usable motion. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
