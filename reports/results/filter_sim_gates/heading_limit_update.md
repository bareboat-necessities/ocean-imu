# Heading-corrected replay limits

The limits qualify the complete 544-cell two-direction CI measurement run on `a02d9fd4c3a882dd50649ed4b1cdd27b93b859c7`. Both the steering model and yaw corrections are active. The existing 0.5% allowance and upward 1e-5 rounding are unchanged. All 272 horizontal cells independently gate geographic wave Y RMS in metres.

Qualifying `baseline_*.json` receipts retain their failures against the earlier limits. The independent enforced run [43994230010](https://github.com/bareboat-necessities/ocean-imu/actions/runs/43994230010) passes all 544 cells using these caps; the numerical records match the qualifying measurements exactly. Accuracy improvements and adverse changes, including held-out draws, are reported in [yaw-regression.md](../heading_kinematics/yaw-regression.md).

## Global default yaw ceilings

Before is the steering-only branch. Each global cap covers both directions; per-record gates are usually lower. Interval-profile caps are also regenerated.

| Family | Before (deg) | After (deg) |
| --- | ---: | ---: |
| carma | 2.52379 | 0.94918 |
| csm | 2.51374 | 0.93871 |
| ou2 | 2.59181 | 1.04316 |
| ou3 | 2.55470 | 1.09866 |
| tangent | 2.50632 | 0.93448 |
| tcarma | 2.51554 | 0.94413 |
| tcsm | 2.50625 | 0.93448 |
| tfg | 2.62044 | 1.05168 |
| tfg_legacy | 3.86696 | 3.75715 |

## Per-record changes

All gated channels are included. Increases reflect measured tradeoffs rather than a larger allowance. The JSON lists the largest increases; the versioned per-record tables retain every cap.

| Metric | Tightened | Increased | Unchanged | New |
| --- | ---: | ---: | ---: | ---: |
| acc_3d | 201 | 71 | 0 | 0 |
| acc_z | 197 | 72 | 3 | 0 |
| displacement | 181 | 91 | 0 | 0 |
| displacement_y | 153 | 106 | 13 | 0 |
| gyro_3d | 162 | 110 | 0 | 0 |
| heave | 161 | 112 | 255 | 0 |
| pitch | 200 | 86 | 2 | 0 |
| roll | 213 | 69 | 6 | 0 |
| yaw | 143 | 161 | 0 | 0 |

## Independent Y caps, 0.27 m JONSWAP

Metres RMS over the final 900 s, production defaults.

| Family | Following | Head |
| --- | ---: | ---: |
| carma | 0.00452 | 0.00234 |
| csm | 0.00427 | 0.00231 |
| ou2 | 0.01262 | 0.00632 |
| ou3 | 0.01310 | 0.00833 |
| tangent | 0.00428 | 0.00231 |
| tcarma | 0.00452 | 0.00234 |
| tcsm | 0.00428 | 0.00231 |
| tfg | 0.01708 | 0.00838 |
| tfg_legacy | 0.03915 | 0.02561 |
