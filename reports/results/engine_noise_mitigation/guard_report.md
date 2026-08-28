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
high-pass detector says there is machinery to remove.

Scoring uses the trailing **900 s** of each 1200 s record,
pooled over the eight stationary records as `sqrt(mean(record_RMS^2))`.

Source commit used for the replay: `d0ce9cbe5ee9701ce3cca18cff9d4a96ff5b99b9`.

## Result

| Condition | Guard | Detector [m/s²] | Engaged | 3-D [m] | 3-D offset [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | 0.5224 | 0.1396 | 0.144 | 0.63 | 1.00 |
| engine off | on | 0.0080 | 0.000 | 0.5224 | 0.1396 | 0.144 | 0.63 | 1.00 |
| 800 rpm | off | — | — | 0.8770 | 0.5115 | 2.014 | 12.61 | 1.68 |
| 800 rpm | on | 0.0867 | 1.000 | 0.5982 | 0.2247 | 1.327 | 5.10 | 1.15 |
| 1600 rpm | off | — | — | 1.4021 | 0.9051 | 1.211 | 17.85 | 2.68 |
| 1600 rpm | on | 0.1182 | 1.000 | 0.5889 | 0.2093 | 1.125 | 4.07 | 1.13 |
| 2400 rpm | off | — | — | 4.2318 | 3.8746 | 3.075 | 45.38 | 8.10 |
| 2400 rpm | on | 0.1443 | 1.000 | 0.5961 | 0.2226 | 1.180 | 4.17 | 1.14 |
| 3200 rpm | off | — | — | 8.6509 | 8.2819 | 3.283 | 66.55 | 16.56 |
| 3200 rpm | on | 0.1803 | 1.000 | 0.6746 | 0.3313 | 1.773 | 6.27 | 1.29 |
| 2400 rpm, quiet mount | off | — | — | 1.3724 | 1.1528 | 2.744 | 15.21 | 2.63 |
| 2400 rpm, quiet mount | on | 0.0724 | 0.831 | 0.5673 | 0.1849 | 0.584 | 2.38 | 1.09 |
| 2400 rpm, engine bed | off | — | — | 37.5310 | 10.6370 | 7.705 | 103.88 | 71.84 |
| 2400 rpm, engine bed | on | 0.2884 | 1.000 | 0.8345 | 0.5241 | 2.260 | 8.36 | 1.60 |
| 2400 rpm, wide sensor | off | — | — | 7.5696 | 6.7859 | 2.360 | 85.70 | 14.49 |
| 2400 rpm, wide sensor | on | 0.1935 | 1.000 | 0.5795 | 0.1892 | 0.520 | 2.08 | 1.11 |

## Transparency with no engine running

The engine-off rows are **identical to every digit**, because the detector never crosses its lower rail: the guard leaves
the measurement path untouched and returns its input unchanged.  This is
the property that lets the guard ship without re-cutting any fitted gate
or invalidating a committed replay.

It holds because the detector is placed above the sea rather than at the
conditioning corner.  Across a 31:1 range of significant wave height the
clean detector reading varies by about one percent, so a big sea does not
look like machinery to it.

## Figure

- `ou_engine_noise_guard.svg`: pooled 3-D error and standing tilt offset with the
  guard off and on, against the engine-off baseline.

Mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.

## What this does not do

The guard removes vibration from the measurement path.  It does not make
the estimator vibration-aware: the accelerometer measurement covariance
is unchanged, so the filter still believes a conditioned sample is as
good as a quiet one.  It also cannot help with machinery whose orders
reach into the wave band, which no front-end filter can separate from the
sea, and it does not touch the residual that survives its own stopband.
