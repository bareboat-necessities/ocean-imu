# Engine-noise mitigation: the OU-III front-end vibration guard

The degradation study showed that machinery vibration costs OU-III a
large and almost entirely systematic error, that the size of it tracks
recorded out-of-band accelerometer power rather than the placement of the
aliased crank orders, and that the gyroscope path contributes nothing.
That points at a single remedy: keep the out-of-band content out of the
accelerometer before anything reads it.

`AccelVibrationGuard` sits at the one point in `updateCore_` where raw
measurements arrive, so the Mahony proxy, the MEKF, and the tilt watchdog
all see the same conditioned signal.  It low-passes the accelerometer in
the empty decade between the wave band and the machinery band
(**2 poles at 14 Hz**), and engages only when a separate
high-pass detector says there is machinery to remove.  The third arm adds
the vibration-aware measurement covariance: the commanded accelerometer
sigma is raised to `sqrt(sigma^2 + (0.75 * excess)^2)` from the same gated
excess, so the covariance and the measurement describe the same conditions.

Scoring uses the trailing **900 s** of each 1200 s record,
pooled over the eight stationary records as `sqrt(mean(record_RMS^2))`.

Source commit used for the replay: `02adc1b0d12b34ff91b1c30952aef2f53d68db4f`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.4047 | 0.117 | 0.51 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.4047 | 0.117 | 0.51 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0405 | 0.4047 | 0.117 | 0.51 | 1.000 |
| 800 rpm | off | — | — | — | 3.9014 | 2.017 | 99.29 | 9.640 |
|  | guard | 0.0867 | 1.000 | — | 0.5833 | 1.514 | 24.34 | 1.441 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 0.5727 | 1.519 | 19.59 | 1.415 |
| 1600 rpm | off | — | — | — | 12.7705 | 3.564 | 90.28 | 31.555 |
|  | guard | 0.1182 | 1.000 | — | 0.5004 | 1.287 | 12.57 | 1.236 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.4615 | 0.885 | 6.64 | 1.140 |
| 2400 rpm | off | — | — | — | 19.7530 | 1.390 | 90.30 | 48.809 |
|  | guard | 0.1443 | 1.000 | — | 0.4997 | 1.163 | 11.60 | 1.235 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.4506 | 0.497 | 4.66 | 1.113 |
| 3200 rpm | off | — | — | — | 20.3989 | 3.396 | 91.15 | 50.405 |
|  | guard | 0.1803 | 1.000 | — | 0.6231 | 1.507 | 22.78 | 1.540 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.4651 | 0.569 | 5.58 | 1.149 |
| 2400 rpm, quiet mount | off | — | — | — | 3.5172 | 1.801 | 90.09 | 8.691 |
|  | guard | 0.0724 | 0.831 | — | 0.4569 | 0.689 | 7.68 | 1.129 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.4406 | 0.443 | 4.81 | 1.089 |
| 2400 rpm, engine bed | off | — | — | — | 127.9366 | 25.594 | 103.92 | 316.127 |
|  | guard | 0.2884 | 1.000 | — | 0.9638 | 1.414 | 64.58 | 2.381 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.4922 | 0.554 | 5.20 | 1.216 |
| 2400 rpm, wide sensor | off | — | — | — | 21.0117 | 1.157 | 90.20 | 51.919 |
|  | guard | 0.1935 | 1.000 | — | 0.4569 | 0.662 | 6.19 | 1.129 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.4432 | 0.058 | 1.35 | 1.095 |

## Transparency with no engine running

The engine-off rows are **identical to every digit**, because the detector never crosses its lower rail: the guard leaves
the measurement path untouched and returns its input unchanged.  This is
the property that lets the guard ship without re-cutting any fitted gate
or invalidating a committed replay.

The per-record output retains detector readings and engagement for all
eight vessel-response cases; these are checked independently of height.

## Residual and scope

Conditioning introduces group delay. The fixed-gain arms measure its
combined effect with vibration attenuation and covariance inflation.
This run does not isolate delay on a forced-on quiet input and does not
optimize the deployed covariance gain 0.75.

## Figure

- `ou_engine_noise_guard.svg`: pooled 3-D error and standing tilt offset across the
  three configurations, against the engine-off baseline.

Mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.

## What this does not do

Group delay is the price of conditioning and is paid whether or not there
is anything left to remove. No front-end filter separates machinery whose
orders reach into the wave band, since there is nothing there to separate
them from vessel motion. Compare each recorded condition and arm directly.

Mechanical isolation and a tighter sensor anti-alias filter still act on
the quantity that matters, and are the only things that reduce the input
rather than manage it.
