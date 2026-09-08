# Direction RAO and engine-guard ablation

Each row is an unweighted mean across eight full-record replays, scored over the final 900 seconds. Directed RMS uses samples where direction is resolved; unresolved fractions must be considered alongside it.

| Family | RPM | Arm | Axis error ° | Directed RMS ° | Correct % | Wrong % | Unresolved % |
|---|---:|---|---:|---:|---:|---:|---:|
| OU-II | 0 | unconditioned | 14.37 | 29.81 | 72.59 | 0.72 | 26.69 |
| OU-II | 0 | guard-only | 14.37 | 29.81 | 72.59 | 0.72 | 26.69 |
| OU-II | 0 | guard-and-rao | 8.36 | 16.93 | 69.55 | 0.23 | 30.22 |
| OU-II | 2400 | unconditioned | 58.37 | nan | 3.41 | 0.64 | 95.94 |
| OU-II | 2400 | guard-only | 16.29 | 49.33 | 59.58 | 2.38 | 38.05 |
| OU-II | 2400 | guard-and-rao | 12.83 | 34.48 | 66.22 | 0.86 | 32.93 |
| OU-III | 0 | unconditioned | 15.75 | 23.02 | 71.21 | 0.43 | 28.36 |
| OU-III | 0 | guard-only | 15.75 | 23.02 | 71.21 | 0.43 | 28.36 |
| OU-III | 0 | guard-and-rao | 10.00 | 17.73 | 69.33 | 0.27 | 30.41 |
| OU-III | 2400 | unconditioned | 58.87 | nan | 2.46 | 0.95 | 96.59 |
| OU-III | 2400 | guard-only | 17.70 | 49.30 | 59.79 | 2.72 | 37.49 |
| OU-III | 2400 | guard-and-rao | 14.19 | 31.16 | 66.43 | 0.66 | 32.92 |

Quiet guard-on/off measurements are identical. RAO-on/off attitude and displacement measurements are identical in every matched pair. The equalizer improves average angular accuracy, while quiet-water unresolved fractions increase. Engine vibration substantially reduces availability without the guard; correction cannot recover direction from absent or weak motion. Quality and confidence thresholds are unchanged. The unconditioned arm disables both the guard and its covariance inflation, so that comparison measures their combined effect.
