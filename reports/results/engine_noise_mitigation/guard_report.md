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

Source commit used for the replay: `6b20d608bca82f6487fd6d4c1d8344f4af829eef`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.4052 | 0.146 | 0.57 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.4052 | 0.146 | 0.57 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0000 | 0.4052 | 0.146 | 0.57 | 1.000 |
| 800 rpm | off | — | — | — | 4.4617 | 2.208 | 104.13 | 11.010 |
|  | guard | 0.0867 | 1.000 | — | 0.5992 | 1.378 | 57.36 | 1.479 |
|  | guard+R | 0.0867 | 1.000 | 0.0517 | 0.5902 | 1.656 | 37.08 | 1.456 |
| 1600 rpm | off | — | — | — | 15.7608 | 3.896 | 103.91 | 38.892 |
|  | guard | 0.1182 | 1.000 | — | 0.5172 | 1.454 | 28.60 | 1.276 |
|  | guard+R | 0.1182 | 1.000 | 0.0724 | 0.4727 | 1.086 | 8.31 | 1.166 |
| 2400 rpm | off | — | — | — | 22.1951 | 1.424 | 103.90 | 54.769 |
|  | guard | 0.1443 | 1.000 | — | 0.5170 | 1.501 | 28.45 | 1.276 |
|  | guard+R | 0.1443 | 1.000 | 0.0906 | 0.4565 | 0.548 | 5.08 | 1.126 |
| 3200 rpm | off | — | — | — | 26.3465 | 4.739 | 103.94 | 65.013 |
|  | guard | 0.1803 | 1.000 | — | 0.6535 | 1.611 | 55.78 | 1.613 |
|  | guard+R | 0.1803 | 1.000 | 0.1165 | 0.4702 | 0.571 | 5.66 | 1.160 |
| 2400 rpm, quiet mount | off | — | — | — | 3.8651 | 1.904 | 103.96 | 9.538 |
|  | guard | 0.0724 | 0.831 | — | 0.4675 | 1.009 | 10.93 | 1.154 |
|  | guard+R | 0.0724 | 0.831 | 0.0433 | 0.4489 | 0.616 | 6.02 | 1.108 |
| 2400 rpm, engine bed | off | — | — | — | 125.2986 | 26.237 | 103.42 | 309.191 |
|  | guard | 0.2884 | 1.000 | — | 0.9943 | 1.485 | 82.58 | 2.454 |
|  | guard+R | 0.2884 | 1.000 | 0.1960 | 0.4883 | 0.541 | 5.17 | 1.205 |
| 2400 rpm, wide sensor | off | — | — | — | 23.5255 | 1.245 | 103.83 | 58.052 |
|  | guard | 0.1935 | 1.000 | — | 0.4688 | 1.030 | 12.53 | 1.157 |
|  | guard+R | 0.1935 | 1.000 | 0.1261 | 0.4405 | 0.036 | 1.31 | 1.087 |

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
