# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Directed RMS uses samples where direction is resolved; unresolved fractions must be considered alongside it. A dash means at least one record has no finite angular score, so the eight-record mean is undefined; missing records are not silently excluded.

| Family | RPM | Arm | Axis error ° | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 14.94 | 28.57 | 71.91 | 0.60 | 27.50 |
| OU-II | 0 | guard-only | 14.94 | 28.57 | 71.91 | 0.60 | 27.50 |
| OU-II | 0 | guard-and-rao | 9.00 | 16.73 | 69.51 | 0.23 | 30.26 |
| OU-II | 2400 | unconditioned | 58.58 | -- | 3.44 | 0.55 | 96.01 |
| OU-II | 2400 | guard-only | 16.49 | 50.29 | 59.43 | 2.52 | 38.05 |
| OU-II | 2400 | guard-and-rao | 13.05 | 33.20 | 66.24 | 0.87 | 32.89 |
| OU-III | 0 | unconditioned | 15.76 | 23.57 | 71.20 | 0.43 | 28.36 |
| OU-III | 0 | guard-only | 15.76 | 23.57 | 71.20 | 0.43 | 28.36 |
| OU-III | 0 | guard-and-rao | 9.98 | 17.74 | 69.32 | 0.26 | 30.42 |
| OU-III | 2400 | unconditioned | 60.19 | -- | 2.87 | 1.29 | 95.84 |
| OU-III | 2400 | guard-only | 17.74 | 49.39 | 59.96 | 2.79 | 37.25 |
| OU-III | 2400 | guard-and-rao | 14.20 | 31.26 | 66.50 | 0.64 | 32.86 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. The equalizer improves average angular accuracy, while quiet-water unresolved fractions increase. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
