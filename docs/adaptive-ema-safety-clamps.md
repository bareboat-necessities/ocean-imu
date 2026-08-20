# Universal safety clamps for sea-scaled adaptation EMAs

The OU-III, OU-II and TFG adaptation paths now share one safety envelope for
**dynamically estimated EMA time scales**.  This is a guard, not another fitted
tuning coefficient.

## Bounds

- dynamic sea/model time base used by an EMA: **0.5 to 6.0 s**
- final dynamically derived EMA horizon: **0.05 to 30.0 s**

The time-base clamp is applied only when a measured `T_sea` or dynamically
estimated OU `tau` is used to construct an EMA horizon.  It does **not** clamp
the physical wave-period estimate or the OU process time constant itself.
Explicit fixed-second ablation setters keep their requested values.

## Why these limits are outside the calibrated envelope

The versioned eight-sea reference family spans approximately `T_z = 2.3..8.4 s`,
therefore `T_sea = T_z/2 = 1.15..4.2 s`.  The deployed dynamic horizons are:

| channel | law | reference range |
| --- | --- | ---: |
| tuner frequency EMA | `0.10 T_sea` | 0.115..0.42 s |
| common `tau/sigma_aw` EMA | `0.40 T_sea` | 0.46..1.68 s |
| tuner variance EMA (`K=2`) | `2 T_z = 4 T_sea` | 4.6..16.8 s |
| OU-III `r_S` EMA | `1.5 tau` | about 1.9..6.3 s |
| OU-II `r_p0/r_v0` EMA | `3 tau` | about 3.8..12.7 s |
| TFG `r_S` EMA | `3 tau` | same order, below 13 s on the reference family |

Thus every deployed reference trajectory is strictly interior to both guards.
The clamps act only on startup/transient estimator excursions or on inputs well
outside the evidenced marine envelope.

## Why clamping only `T_sea` is almost, but not quite, enough

Clamping the dynamic sea time automatically protects the tuner frequency EMA,
the common `tau/sigma_aw` EMA and the K-period variance EMA.  The drift-channel
smoothers need one additional final guard because their helper can optionally
shorten its horizon with the discrepancy/slew term *after* the sea/model time
has been scaled.  The shared 0.05..30 s outer guard catches that path too.

All three filter families call the same helpers, so the limits cannot silently
diverge between OU-III, OU-II and TFG.
