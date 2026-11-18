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

Source commit used for the replay: `fee96e2b3be2d3d9d67b654f76ebea92b291df57`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.3065 | 0.161 | 0.88 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.3065 | 0.161 | 0.88 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0331 | 0.3065 | 0.161 | 0.88 | 1.000 |
| 800 rpm | off | — | — | — | 2.5511 | 0.547 | 107.29 | 8.323 |
|  | guard | 0.0867 | 1.000 | — | 0.4028 | 0.957 | 9.08 | 1.314 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 0.3976 | 0.941 | 8.58 | 1.297 |
| 1600 rpm | off | — | — | — | 5.6502 | 1.881 | 99.11 | 18.433 |
|  | guard | 0.1182 | 1.000 | — | 0.3681 | 0.646 | 5.28 | 1.201 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.3489 | 0.384 | 3.38 | 1.138 |
| 2400 rpm | off | — | — | — | 8.5609 | 0.639 | 98.82 | 27.929 |
|  | guard | 0.1443 | 1.000 | — | 0.3588 | 0.512 | 4.78 | 1.171 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.3454 | 0.348 | 2.85 | 1.127 |
| 3200 rpm | off | — | — | — | 10.0197 | 1.143 | 103.56 | 32.688 |
|  | guard | 0.1803 | 1.000 | — | 0.3861 | 0.786 | 8.17 | 1.260 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.3590 | 0.582 | 4.47 | 1.171 |
| 2400 rpm, quiet mount | off | — | — | — | 3.2844 | 0.714 | 90.44 | 10.715 |
|  | guard | 0.0724 | 0.831 | — | 0.3582 | 0.649 | 4.70 | 1.168 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.3464 | 0.491 | 3.49 | 1.130 |
| 2400 rpm, engine bed | off | — | — | — | 97.1712 | 22.879 | 102.88 | 317.005 |
|  | guard | 0.2884 | 1.000 | — | 0.5014 | 1.238 | 17.23 | 1.636 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.3772 | 0.643 | 4.74 | 1.231 |
| 2400 rpm, wide sensor | off | — | — | — | 8.9115 | 0.529 | 98.41 | 29.072 |
|  | guard | 0.1935 | 1.000 | — | 0.3532 | 0.453 | 3.60 | 1.152 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.3453 | 0.228 | 1.45 | 1.126 |

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
