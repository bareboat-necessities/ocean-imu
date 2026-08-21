# OU-III joint period/statistics confirmation

seeds=10
baseline=baseline (K_T=8, K_log=0.50, K_sigma=2)
candidates=joint_log_005 (K_T=4, K_log=0.05, K_sigma=4); joint_log_010 (K_T=4, K_log=0.10, K_sigma=4)

Paired geometric-mean ratios vs production baseline; negative is better.

| candidate | stationary Z | stationary 3D | blend Z | recover Z | mean Tz [s] | mean sigma [m/s2] |
|---|---:|---:|---:|---:|---:|---:|
| baseline (K_T=8, K_log=0.50, K_sigma=2) | +0.000% [+0.000,+0.000] | +0.000% [+0.000,+0.000] | +0.000% [+0.000,+0.000] | +0.000% [+0.000,+0.000] | 5.4721 | 0.9849 |
| joint_log_005 (K_T=4, K_log=0.05, K_sigma=4) | -1.063% [-1.257,-0.868] | -0.162% [-0.295,-0.030] | -12.274% [-13.621,-10.905] | -0.267% [-1.408,+0.887] | 5.5417 | 0.9643 |
| joint_log_010 (K_T=4, K_log=0.10, K_sigma=4) | -0.936% [-1.127,-0.744] | -0.065% [-0.202,+0.073] | -12.566% [-13.947,-11.164] | -0.229% [-1.397,+0.953] | 5.5499 | 0.9639 |

Direct paired ratio: joint_log_010 versus joint_log_005; negative favors 0.10.

| comparison | stationary Z | stationary 3D | blend Z | recover Z |
|---|---:|---:|---:|---:|
| 0.10 / 0.05 | +0.128% [+0.105,+0.152] | +0.098% [+0.079,+0.116] | -0.334% [-0.431,-0.236] | +0.038% [-0.062,+0.138] |
