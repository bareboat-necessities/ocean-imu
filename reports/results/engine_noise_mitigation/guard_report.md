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

Source commit used for the replay: `9de0543507509cabb478b751b83bf5c74711764d`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.4052 | 0.146 | 0.58 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.4052 | 0.146 | 0.58 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0405 | 0.4052 | 0.146 | 0.58 | 1.000 |
| 800 rpm | off | — | — | — | 4.3703 | 2.055 | 99.76 | 10.786 |
|  | guard | 0.0867 | 1.000 | — | 0.5987 | 1.511 | 24.23 | 1.478 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 0.5898 | 1.517 | 19.50 | 1.456 |
| 1600 rpm | off | — | — | — | 15.7470 | 3.781 | 90.24 | 38.864 |
|  | guard | 0.1182 | 1.000 | — | 0.5168 | 1.283 | 12.48 | 1.276 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.4727 | 0.871 | 6.53 | 1.167 |
| 2400 rpm | off | — | — | — | 22.1857 | 1.411 | 90.29 | 54.755 |
|  | guard | 0.1443 | 1.000 | — | 0.5167 | 1.159 | 11.51 | 1.275 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.4564 | 0.474 | 4.54 | 1.127 |
| 3200 rpm | off | — | — | — | 26.3297 | 4.665 | 91.24 | 64.982 |
|  | guard | 0.1803 | 1.000 | — | 0.6525 | 1.510 | 22.70 | 1.610 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.4702 | 0.545 | 5.45 | 1.160 |
| 2400 rpm, quiet mount | off | — | — | — | 3.8638 | 1.828 | 90.15 | 9.536 |
|  | guard | 0.0724 | 0.831 | — | 0.4674 | 0.674 | 7.56 | 1.154 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.4488 | 0.423 | 4.69 | 1.108 |
| 2400 rpm, engine bed | off | — | — | — | 125.2212 | 25.126 | 103.68 | 309.048 |
|  | guard | 0.2884 | 1.000 | — | 0.9946 | 1.419 | 64.71 | 2.455 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.4883 | 0.528 | 5.07 | 1.205 |
| 2400 rpm, wide sensor | off | — | — | — | 23.5226 | 1.129 | 90.15 | 58.054 |
|  | guard | 0.1935 | 1.000 | — | 0.4687 | 0.648 | 6.09 | 1.157 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.4405 | 0.061 | 1.26 | 1.087 |

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
