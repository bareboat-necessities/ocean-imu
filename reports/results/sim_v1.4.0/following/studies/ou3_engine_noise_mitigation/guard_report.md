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
| engine off | off | — | — | — | 0.3876 | 0.091 | 0.66 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.3876 | 0.091 | 0.66 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0478 | 0.3876 | 0.091 | 0.66 | 1.000 |
| 800 rpm | off | — | — | — | 13.2706 | 1.521 | 90.30 | 34.235 |
|  | guard | 0.0867 | 1.000 | — | 2.1219 | 1.093 | 91.91 | 5.474 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 0.9886 | 1.810 | 36.48 | 2.550 |
| 1600 rpm | off | — | — | — | 29.7766 | 1.773 | 90.06 | 76.815 |
|  | guard | 0.1182 | 1.000 | — | 0.7416 | 1.676 | 23.88 | 1.913 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.4827 | 1.095 | 8.00 | 1.245 |
| 2400 rpm | off | — | — | — | 135.8790 | 6.307 | 90.34 | 350.530 |
|  | guard | 0.1443 | 1.000 | — | 0.6991 | 1.606 | 24.21 | 1.804 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.4146 | 0.530 | 4.41 | 1.070 |
| 3200 rpm | off | — | — | — | 114.7363 | 11.957 | 90.59 | 295.988 |
|  | guard | 0.1803 | 1.000 | — | 1.5269 | 1.387 | 69.93 | 3.939 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.4161 | 0.437 | 3.83 | 1.073 |
| 2400 rpm, quiet mount | off | — | — | — | 5.2467 | 0.475 | 90.11 | 13.535 |
|  | guard | 0.0724 | 0.831 | — | 0.4737 | 0.981 | 9.74 | 1.222 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.4137 | 0.596 | 5.22 | 1.067 |
| 2400 rpm, engine bed | off | — | — | — | 193.8073 | 6.972 | 91.51 | 499.968 |
|  | guard | 0.2884 | 1.000 | — | 2.2532 | 0.397 | 89.52 | 5.813 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.4416 | 0.405 | 3.34 | 1.139 |
| 2400 rpm, wide sensor | off | — | — | — | 34.6177 | 0.217 | 89.92 | 89.304 |
|  | guard | 0.1935 | 1.000 | — | 0.4505 | 0.921 | 10.03 | 1.162 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.4007 | 0.057 | 0.92 | 1.034 |

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
