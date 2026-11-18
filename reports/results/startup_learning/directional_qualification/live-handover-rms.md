# Live OU-III handover components

These replays change the state supplied at handover and run the production OU-III filter on the same records and sensor draws. Each component restoration starts from the deployed handover. The candidate accelerometer-bias seed uses only sensors. The early cohort has 32 matched starts per policy and the long cohort has 48. Each includes both sea directions.

## First second after handover

| Handover | Tilt (°) | Yaw (°) | Gyro bias (rad/s) | Accelerometer bias (m/s²) |
|---|---:|---:|---:|---:|
| Previous handover | 0.8439 | 4.8111 | 0.0007818 | 0.078587 |
| Deployed handover | 0.5436 | 3.7367 | 0.0004233 | 0.078587 |
| Restore previous tilt | 0.7148 | 3.8927 | 0.0004235 | 0.078587 |
| Restore previous yaw | 0.6044 | 4.5959 | 0.0004232 | 0.078587 |
| Restore previous field | 0.5496 | 3.7350 | 0.0004233 | 0.078587 |
| Restore previous gyro bias | 0.5444 | 3.7374 | 0.0007815 | 0.078587 |
| Add gravity-axis accel-bias seed | 0.5459 | 3.7363 | 0.0004233 | 0.065268 |

## Late RMS from 300 to 1,200 seconds

| Handover | X (m) | Y (m) | Z (m) | 3-D (m) | Accelerometer bias (m/s²) | Yaw (°) |
|---|---:|---:|---:|---:|---:|---:|
| Previous handover | 0.216155 | 0.180696 | 0.197872 | 0.344278 | 0.039893 | 0.817777 |
| Deployed handover | 0.216192 | 0.180659 | 0.197870 | 0.344280 | 0.039885 | 0.816055 |
| Restore previous tilt | 0.216190 | 0.180666 | 0.197870 | 0.344283 | 0.039882 | 0.816249 |
| Restore previous yaw | 0.216170 | 0.180686 | 0.197871 | 0.344282 | 0.039886 | 0.818338 |
| Restore previous field | 0.216186 | 0.180652 | 0.197870 | 0.344273 | 0.039883 | 0.815908 |
| Restore previous gyro bias | 0.216186 | 0.180667 | 0.197871 | 0.344282 | 0.039896 | 0.815175 |
| Add gravity-axis accel-bias seed | 0.216129 | 0.180692 | 0.197949 | 0.344305 | 0.039834 | 0.817760 |

## Coverage and trade-offs

In the early cohort, 13 of 16 following-sea starts and all 16 head-sea starts met 35 seconds. The corresponding long-cohort counts were 20 of 24 and 24 of 24. All seven policies used the same handover epoch and deadline outcome within each matched draw. In the original 960-start qualification, all 72 late following-sea cases used the proxy fallback and went live by 40.270 seconds. [All 960 handover epochs and corrected gyro checks](../field_steps/all-handover-gyro-continuity.json)

The candidate accelerometer-bias seed lowered first-second pooled bias RMS from 0.078587 to 0.065268 m/s². Late pooled bias RMS went from 0.039885 to 0.039834 m/s², while Y, 3-D and yaw RMS increased. The calibrated prior remains the production default.

| Candidate at 300–1,200 s | 3-D higher / pairs | Bias higher / pairs | Yaw higher / pairs |
|---|---:|---:|---:|
| Following | 7 / 24 | 10 / 24 | 16 / 24 |
| Head | 7 / 24 | 6 / 24 | 9 / 24 |

The horizontal wave readout reduced long-run X, Y and 3-D RMS in each of the 336 matched policy replays against the preceding estimator. The early record tables matched exactly. The first-second and late tables above compare the handover policies with that readout enabled.

[All early and late RMS by sea and policy](live-handover-arms.json) · [Paired candidate loss counts and worst draws](gravity-axis-bias-tradeoffs.json) · [Original six-job workflow](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44373741328)
