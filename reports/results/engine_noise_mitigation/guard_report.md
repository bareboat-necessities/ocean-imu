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

Source commit used for the replay: `cfb9c8067b5b5f0fb4e04bb022fb76b1fee38053`.

## Result

| Condition | Arm | Detector [m/s²] | Engaged | Racc σ [m/s²] | 3-D [m] | Pitch offset [deg] | Yaw [deg] | vs engine-off |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| engine off | off | — | — | — | 0.5224 | 0.144 | 0.63 | 1.000 |
|  | guard | 0.0080 | 0.000 | — | 0.5224 | 0.144 | 0.63 | 1.000 |
|  | guard+R | 0.0080 | 0.000 | 0.0000 | 0.5224 | 0.144 | 0.63 | 1.000 |
| 800 rpm | off | — | — | — | 0.8770 | 2.014 | 12.61 | 1.679 |
|  | guard | 0.0867 | 1.000 | — | 0.5982 | 1.327 | 5.10 | 1.145 |
|  | guard+R | 0.0867 | 1.000 | 0.0517 | 0.6024 | 1.369 | 5.00 | 1.153 |
| 1600 rpm | off | — | — | — | 1.4021 | 1.211 | 17.85 | 2.684 |
|  | guard | 0.1182 | 1.000 | — | 0.5889 | 1.125 | 4.07 | 1.127 |
|  | guard+R | 0.1182 | 1.000 | 0.0724 | 0.5827 | 0.814 | 2.97 | 1.115 |
| 2400 rpm | off | — | — | — | 4.2318 | 3.075 | 45.38 | 8.100 |
|  | guard | 0.1443 | 1.000 | — | 0.5961 | 1.180 | 4.17 | 1.141 |
|  | guard+R | 0.1443 | 1.000 | 0.0906 | 0.5768 | 0.444 | 1.84 | 1.104 |
| 3200 rpm | off | — | — | — | 8.6509 | 3.283 | 66.55 | 16.559 |
|  | guard | 0.1803 | 1.000 | — | 0.6746 | 1.773 | 6.27 | 1.291 |
|  | guard+R | 0.1803 | 1.000 | 0.1165 | 0.5876 | 0.551 | 2.16 | 1.125 |
| 2400 rpm, quiet mount | off | — | — | — | 1.3724 | 2.744 | 15.21 | 2.627 |
|  | guard | 0.0724 | 0.831 | — | 0.5673 | 0.584 | 2.38 | 1.086 |
|  | guard+R | 0.0724 | 0.831 | 0.0433 | 0.5631 | 0.369 | 1.67 | 1.078 |
| 2400 rpm, engine bed | off | — | — | — | 37.5310 | 7.705 | 103.88 | 71.839 |
|  | guard | 0.2884 | 1.000 | — | 0.8345 | 2.260 | 8.36 | 1.597 |
|  | guard+R | 0.2884 | 1.000 | 0.1960 | 0.6110 | 0.567 | 2.16 | 1.170 |
| 2400 rpm, wide sensor | off | — | — | — | 7.5696 | 2.360 | 85.70 | 14.489 |
|  | guard | 0.1935 | 1.000 | — | 0.5795 | 0.520 | 2.08 | 1.109 |
|  | guard+R | 0.1935 | 1.000 | 0.1261 | 0.5716 | 0.047 | 0.69 | 1.094 |

## Transparency with no engine running

The engine-off rows are **identical to every digit**, because the detector never crosses its lower rail: the guard leaves
the measurement path untouched and returns its input unchanged.  This is
the property that lets the guard ship without re-cutting any fitted gate
or invalidating a committed replay.

It holds because the detector is placed above the sea rather than at the
conditioning corner.  Across a 31:1 range of significant wave height the
clean detector reading varies by about one percent, so a big sea does not
look like machinery to it.

## Why the residual does not go to one

Conditioning costs group delay, and that cost is still there when there is
nothing left to remove.  Forcing the guard on over a *quiet* input isolates
it: on the two JONSWAP records at Hs 1.5 and 8.5 the delay alone accounts
for 1.021x and 1.063x, against deployed residuals of 1.074x and 1.151x.
So roughly half the remaining gap is the guard's own delay, which no
covariance change can touch, and the covariance arm attacks the other half.

That is also why the two channels disagree about the best gain.  Attitude
keeps improving as the accelerometer is de-weighted further, but the
accelerometer is the only wave measurement there is, so past a gain of
about 1.25 displacement turns back up as the estimate leans on the OU
prior instead.  0.75 sits at the displacement optimum with margin.

## Figure

- `ou_engine_noise_guard.svg`: pooled 3-D error and standing tilt offset across the
  three configurations, against the engine-off baseline.

Mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.

## What this does not do

Group delay is the price of conditioning and is paid whether or not there
is anything left to remove, so the residual cannot reach 1.00 while the
guard is engaged.  And no front-end filter helps with machinery whose
orders reach into the wave band, since there is nothing there to separate
them from the sea: the 800 rpm row is that limit showing itself early,
and it is the one condition where the covariance stage does not pay.

Mechanical isolation and a tighter sensor anti-alias filter still act on
the quantity that matters, and are the only things that reduce the input
rather than manage it.
