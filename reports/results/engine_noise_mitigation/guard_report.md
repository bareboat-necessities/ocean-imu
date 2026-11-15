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

Source commit used for the replay: `f7cd0fa140a2f9647ef33d16bf40c25258442158`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.4939 | 0.133 | 0.62 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.4939 | 0.133 | 0.62 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0515 | 0.4939 | 0.133 | 0.62 | 1.000 |
| 800 rpm | off | — | — | — | 18.1502 | 2.497 | 90.17 | 36.751 |
|  | guard | 0.0867 | 1.000 | — | 2.0790 | 1.265 | 91.38 | 4.209 |
|  | guard+R | 0.0867 | 1.000 | 0.0701 | 1.1938 | 1.875 | 40.54 | 2.417 |
| 1600 rpm | off | — | — | — | 113.5495 | 15.826 | 89.11 | 229.915 |
|  | guard | 0.1182 | 1.000 | — | 0.8993 | 1.746 | 25.75 | 1.821 |
|  | guard+R | 0.1182 | 1.000 | 0.0881 | 0.6249 | 1.171 | 9.09 | 1.265 |
| 2400 rpm | off | — | — | — | 58.8752 | 5.925 | 90.00 | 119.210 |
|  | guard | 0.1443 | 1.000 | — | 0.8560 | 1.633 | 26.07 | 1.733 |
|  | guard+R | 0.1443 | 1.000 | 0.1044 | 0.5428 | 0.537 | 5.10 | 1.099 |
| 3200 rpm | off | — | — | — | 114.3476 | 18.241 | 90.29 | 231.531 |
|  | guard | 0.1803 | 1.000 | — | 1.9660 | 1.326 | 75.36 | 3.981 |
|  | guard+R | 0.1803 | 1.000 | 0.1282 | 0.5434 | 0.428 | 4.40 | 1.100 |
| 2400 rpm, quiet mount | off | — | — | — | 5.8892 | 1.197 | 89.83 | 11.925 |
|  | guard | 0.0724 | 0.831 | — | 0.6126 | 1.072 | 11.04 | 1.240 |
|  | guard+R | 0.0724 | 0.831 | 0.0630 | 0.5392 | 0.633 | 6.10 | 1.092 |
| 2400 rpm, engine bed | off | — | — | — | 135.4851 | 15.610 | 91.99 | 274.330 |
|  | guard | 0.2884 | 1.000 | — | 2.1605 | 0.723 | 89.62 | 4.375 |
|  | guard+R | 0.2884 | 1.000 | 0.2037 | 0.5764 | 0.399 | 3.92 | 1.167 |
| 2400 rpm, wide sensor | off | — | — | — | 40.0888 | 0.603 | 89.91 | 81.172 |
|  | guard | 0.1935 | 1.000 | — | 0.5826 | 1.023 | 11.10 | 1.180 |
|  | guard+R | 0.1935 | 1.000 | 0.1371 | 0.5191 | 0.079 | 1.31 | 1.051 |

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
