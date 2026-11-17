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
| engine off | off | — | — | — | 0.3751 | 0.161 | 0.89 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.3751 | 0.161 | 0.89 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0331 | 0.3751 | 0.161 | 0.89 | 1.000 |
| 800 rpm | off | — | — | — | 2.6143 | 0.545 | 105.88 | 6.970 |
|  | guard | 0.0867 | 1.000 | — | 0.4770 | 0.957 | 9.08 | 1.272 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 0.4726 | 0.941 | 8.59 | 1.260 |
| 1600 rpm | off | — | — | — | 5.6553 | 1.940 | 94.10 | 15.077 |
|  | guard | 0.1182 | 1.000 | — | 0.4410 | 0.646 | 5.28 | 1.176 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.4218 | 0.384 | 3.39 | 1.125 |
| 2400 rpm | off | — | — | — | 8.5337 | 0.609 | 98.91 | 22.752 |
|  | guard | 0.1443 | 1.000 | — | 0.4311 | 0.512 | 4.78 | 1.149 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.4190 | 0.348 | 2.85 | 1.117 |
| 3200 rpm | off | — | — | — | 10.1786 | 1.229 | 103.82 | 27.137 |
|  | guard | 0.1803 | 1.000 | — | 0.4598 | 0.786 | 8.17 | 1.226 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.4356 | 0.582 | 4.47 | 1.161 |
| 2400 rpm, quiet mount | off | — | — | — | 3.3272 | 0.718 | 90.36 | 8.871 |
|  | guard | 0.0724 | 0.831 | — | 0.4309 | 0.650 | 4.70 | 1.149 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.4186 | 0.491 | 3.50 | 1.116 |
| 2400 rpm, engine bed | off | — | — | — | 89.6585 | 19.259 | 103.33 | 239.036 |
|  | guard | 0.2884 | 1.000 | — | 0.5762 | 1.241 | 17.27 | 1.536 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.4594 | 0.643 | 4.74 | 1.225 |
| 2400 rpm, wide sensor | off | — | — | — | 8.9598 | 0.529 | 98.61 | 23.887 |
|  | guard | 0.1935 | 1.000 | — | 0.4250 | 0.453 | 3.60 | 1.133 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.4199 | 0.228 | 1.45 | 1.120 |

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
