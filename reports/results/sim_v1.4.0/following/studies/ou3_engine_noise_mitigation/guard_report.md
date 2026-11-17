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

Source commit used for the replay: `37d4ddb7a2d8431f91d2807c526f6adc4e22b1f8`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.4853 | 0.091 | 0.66 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.4853 | 0.091 | 0.66 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0478 | 0.4853 | 0.091 | 0.66 | 1.000 |
| 800 rpm | off | — | — | — | 13.5467 | 1.521 | 90.26 | 27.913 |
|  | guard | 0.0867 | 1.000 | — | 2.2667 | 1.097 | 91.56 | 4.671 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 1.1247 | 1.811 | 36.54 | 2.317 |
| 1600 rpm | off | — | — | — | 26.4915 | 1.976 | 90.10 | 54.586 |
|  | guard | 0.1182 | 1.000 | — | 0.8571 | 1.677 | 23.90 | 1.766 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.5909 | 1.096 | 8.00 | 1.218 |
| 2400 rpm | off | — | — | — | 100.2985 | 5.996 | 89.81 | 206.667 |
|  | guard | 0.1443 | 1.000 | — | 0.8144 | 1.607 | 24.23 | 1.678 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.5238 | 0.530 | 4.42 | 1.079 |
| 3200 rpm | off | — | — | — | 145.3359 | 10.496 | 90.53 | 299.468 |
|  | guard | 0.1803 | 1.000 | — | 1.7244 | 1.367 | 71.13 | 3.553 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.5271 | 0.438 | 3.83 | 1.086 |
| 2400 rpm, quiet mount | off | — | — | — | 5.4088 | 0.476 | 90.09 | 11.145 |
|  | guard | 0.0724 | 0.831 | — | 0.5813 | 0.982 | 9.74 | 1.198 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.5195 | 0.597 | 5.23 | 1.070 |
| 2400 rpm, engine bed | off | — | — | — | 197.8592 | 6.252 | 91.98 | 407.693 |
|  | guard | 0.2884 | 1.000 | — | 2.3274 | 0.394 | 89.56 | 4.796 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.5609 | 0.406 | 3.34 | 1.156 |
| 2400 rpm, wide sensor | off | — | — | — | 34.5222 | 0.227 | 89.95 | 71.134 |
|  | guard | 0.1935 | 1.000 | — | 0.5590 | 0.924 | 10.05 | 1.152 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.5080 | 0.057 | 0.93 | 1.047 |

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
