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

Source commit used for the replay: `456430d4d2cb2449dd0ae3c7781f71c5eedf599e`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.4046 | 0.157 | 0.57 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.4046 | 0.157 | 0.57 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0000 | 0.4046 | 0.157 | 0.57 | 1.000 |
| 800 rpm | off | — | — | — | 5.2136 | 2.699 | 104.08 | 12.885 |
|  | guard | 0.0867 | 1.000 | — | 0.5982 | 1.382 | 57.58 | 1.478 |
|  | guard+R | 0.0867 | 1.000 | 0.0517 | 0.5915 | 1.680 | 37.15 | 1.462 |
| 1600 rpm | off | — | — | — | 19.7699 | 4.837 | 103.96 | 48.858 |
|  | guard | 0.1182 | 1.000 | — | 0.5122 | 1.479 | 28.62 | 1.266 |
|  | guard+R | 0.1182 | 1.000 | 0.0724 | 0.4695 | 1.104 | 8.35 | 1.160 |
| 2400 rpm | off | — | — | — | 27.2576 | 1.166 | 103.88 | 67.363 |
|  | guard | 0.1443 | 1.000 | — | 0.5140 | 1.526 | 28.54 | 1.270 |
|  | guard+R | 0.1443 | 1.000 | 0.0906 | 0.4534 | 0.557 | 5.10 | 1.121 |
| 3200 rpm | off | — | — | — | 31.3423 | 6.295 | 103.89 | 77.457 |
|  | guard | 0.1803 | 1.000 | — | 0.6695 | 1.634 | 55.78 | 1.655 |
|  | guard+R | 0.1803 | 1.000 | 0.1165 | 0.4712 | 0.579 | 5.68 | 1.164 |
| 2400 rpm, quiet mount | off | — | — | — | 4.5247 | 2.129 | 103.92 | 11.182 |
|  | guard | 0.0724 | 0.831 | — | 0.4603 | 1.036 | 10.98 | 1.137 |
|  | guard+R | 0.0724 | 0.831 | 0.0433 | 0.4415 | 0.634 | 6.05 | 1.091 |
| 2400 rpm, engine bed | off | — | — | — | 112.3932 | 27.927 | 103.60 | 277.761 |
|  | guard | 0.2884 | 1.000 | — | 1.0629 | 1.549 | 82.84 | 2.627 |
|  | guard+R | 0.2884 | 1.000 | 0.1960 | 0.4950 | 0.546 | 5.18 | 1.223 |
| 2400 rpm, wide sensor | off | — | — | — | 27.7936 | 1.224 | 103.73 | 68.687 |
|  | guard | 0.1935 | 1.000 | — | 0.4582 | 1.055 | 12.58 | 1.132 |
|  | guard+R | 0.1935 | 1.000 | 0.1261 | 0.4397 | 0.036 | 1.30 | 1.087 |

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
