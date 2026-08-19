# Tuning the dynamically scaled EMA smoothers on the drift-correction channels

Three channels in the two OU filters smooth their adaptation target with an
exponential moving average whose time constant is proportional to the current
OU operating point rather than fixed in seconds:

| family | channel | target law | horizon |
| --- | --- | --- | --- |
| OU-III | `r_S` | `R_S_coeff * sigma_a * tau^3` | `ADAPT_RS_MULT * tau_target` |
| OU-II | `r_p0` | `R_p0_coeff * sigma_a * tau^2` | `ADAPT_R_p0_MULT * tau_target` |
| OU-II | `r_v0` | `R_v0_coeff * sigma_a * tau` | `ADAPT_R_v0_MULT * tau_target` |

All three multipliers shipped at `5.0`. That value had never been measured
against an alternative. This document is the measurement.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_ou_ii/kalman_ou_ii-sim` and
`tests/kalman_ou_iii/kalman_ou_iii-sim`, scoring the trailing 900 s.
`tools/ou_ema_adapt_study.py` reproduces every table.

**Result.** All three multipliers move from `5.0` to `3.0`. The stationary
records are nearly indifferent; the gain is on sea-state transitions, where the
shipped horizon lagged the target by 10 to 20 s.

**Result, second pass.** OU-III's `ADAPT_RS_MULT` moves again, from `3.0` to
`1.5`, when the transition instrument is sharpened: the crossfade the study
blends over went from 360 s to 120 s, which is fast enough for the smoothing
horizon to actually lag it. Section 7 is that measurement. OU-II's two
multipliers and TFG's are untouched and stay at `3.0`; the same instrument is
available for them, but this pass measured OU-III.

## 1. Why the stationary records cannot answer the question

A smoothing horizon is only observable while its target moves. On a stationary
record the target settles after warmup and every horizon converges to the same
value, so the whole sweep collapses into the noise:

| `ADAPT_RS_MULT` | 3D RMS | Z RMS | roll/pitch | accel bias |
| --- | --- | --- | --- | --- |
| 1.0 | 0.9992 | 0.9966 | 1.0033 | 1.0003 |
| 2.0 | 0.9996 | 0.9990 | 1.0064 | 1.0052 |
| 3.0 | 0.9999 | 0.9997 | 1.0039 | 1.0033 |
| 4.0 | 1.0000 | 0.9999 | 1.0018 | 1.0016 |
| 5.0 (shipped) | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 8.0 | 0.9998 | 0.9999 | 1.0003 | 1.0014 |

Ratios to the shipped multiplier on the same record and seed, pooled by
geometric mean over the four JONSWAP and four PM-Stokes records at three IMU
seed triplets (n = 24). Every entry is inside ±0.6%. A 16x change in the
smoothing horizon is worth less than a percent, which is the expected answer
for an instrument that holds the target still.

## 2. The instrument: synthesized sea-state transitions

`tools/ou_ema_adapt_study.py transition` builds non-stationary records with the
same C2 quintic crossfade `tools/ou_validation.py` uses for its transition
scenario, with displacement, velocity and acceleration all carrying the blend's
derivatives so the record stays kinematically consistent. Sections 3-6 below
were measured when that crossfade ran over 420-780 s of a 1200 s record; it now
runs over 540-660 s, and section 7 is the re-measurement. The 900 s scoring
window opens at 300 s, so it contains 240 s of the start sea, the whole 120 s
crossfade, and 540 s of the endpoint sea; the intervals are scored separately
through `W3D_VALIDATION_SEGMENTS`.

Two endpoint pairs are used, each run in both directions (sea building and sea
decaying) at five wave-phase seeds:

| pair | start sea | endpoint sea | mean `tau_applied` |
| --- | --- | --- | --- |
| `large` | JONSWAP H_s = 1.5 m | H_s = 8.5 m record rescaled to 4.0 m | 3.15 s |
| `small` | JONSWAP H_s = 0.27 m | H_s = 1.5 m record rescaled to 0.7 m | 1.77 s |

Across the `large` transition `tau` roughly doubles, and `r_S ~ sigma_a tau^3`
moves by more than an order of magnitude. That is the excursion the smoother
has to follow.

## 3. What the shipped horizon costs

Sweeping `ADAPT_RS_MULT` on the `large` pair, scoring the crossfade interval
only (10 records x 3 IMU seed triplets, n = 30, ratios to the shipped value):

| `ADAPT_RS_MULT` | 3D RMS | Z RMS | Z RMS worst | roll/pitch | accel bias |
| --- | --- | --- | --- | --- | --- |
| 1.0 | 0.9850 | 0.9796 | 1.0130 | 1.0025 | 1.0077 |
| 1.5 | 0.9872 | 0.9832 | 1.0107 | 0.9917 | 0.9926 |
| 2.0 | 0.9891 | 0.9865 | 1.0079 | 0.9976 | 0.9993 |
| 2.5 | 0.9909 | 0.9892 | 1.0052 | 0.9975 | 0.9986 |
| 3.0 | 0.9927 | 0.9916 | 1.0030 | 0.9976 | 0.9980 |
| 4.0 | 0.9963 | 0.9960 | 1.0004 | 0.9994 | 0.9995 |
| 5.0 (shipped) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 8.0 | 1.0185 | 1.0201 | 1.0422 | 1.0099 | 1.0098 |
| 20 | 1.0596 | 1.0621 | 1.1268 | 1.0262 | 1.0246 |
| 50 | 1.1176 | 1.1242 | 1.2809 | 1.0305 | 1.0209 |

The displacement error is monotone in the horizon over two decades of
multiplier. It flattens below about 1.5, where roll/pitch and the accelerometer
bias start to degrade instead: at that point the smoother has stopped rejecting
tuner jitter and is passing it into `r_S`.

OU-II behaves the same way, with the two channels separating cleanly. Sweeping
one multiplier at a time on the `large` pair, crossfade interval, one IMU seed
triplet (n = 10):

| multiplier | 3D RMS (`r_p0` swept) | 3D RMS (`r_v0` swept) | accel bias (`r_v0` swept) |
| --- | --- | --- | --- |
| 0.35 | 0.9827 | 0.9998 | 0.9523 |
| 1.0 | 0.9877 | 1.0005 | 0.9028 |
| 2.0 | 0.9922 | 1.0003 | 0.9547 |
| 5.0 (shipped) | 1.0000 | 1.0000 | 1.0000 |
| 12 | 1.0152 | 1.0013 | 1.0105 |
| 20 | 1.0297 | 1.0032 | 0.9818 |

`r_p0` owns the displacement error and moves little else; `r_v0` does not move
displacement at all but owns roll/pitch and the accelerometer bias, with an
optimum near 1-2 rather than at the short end. That split follows
the two laws: `r_p0 ~ tau^2` prices position drift, `r_v0 ~ tau` prices velocity
drift, and the velocity channel is what the bias estimator competes with.

## 4. The scaling law is right; the constant was not

The horizon is proportional to `tau`, i.e. it is a fixed number of wave periods
rather than a fixed number of seconds. Repeating the OU-III sweep on the
`small` pair, whose `tau` is 1.8x smaller, separates the two candidate laws: a
fixed-seconds law would want a multiplier 1.8x larger there.

Crossfade interval, one IMU seed triplet, n = 10 per pair:

| `ADAPT_RS_MULT` | blend 3D RMS, `large` pair | blend 3D RMS, `small` pair |
| --- | --- | --- |
| 0.5 | 0.9808 | 0.9920 |
| 1.0 | 0.9837 | 0.9917 |
| 1.5 | 0.9861 | 0.9927 |
| 2.0 | 0.9880 | 0.9938 |
| 3.0 | 0.9920 | 0.9960 |
| 5.0 (shipped) | 1.0000 | 1.0000 |
| 20 | 1.0544 | 1.0223 |

Both pairs are minimized by the same multiplier to within the resolution of the
sweep, and the optimum does not shift toward larger multipliers in the shorter
sea. The tau-proportional law is not refuted; only its constant was wrong.

At the deployed operating points the constant is what decides whether the
smoother is a filter or a lag:

| record | `tau_applied` | horizon at 5.0 | horizon at 3.0 |
| --- | --- | --- | --- |
| JONSWAP H_s = 0.27 m | 1.28 s | 6.4 s | 3.8 s |
| JONSWAP H_s = 1.5 m | 2.18 s | 10.9 s | 6.5 s |
| JONSWAP H_s = 4.0 m | 3.59 s | 18.0 s | 10.8 s |
| JONSWAP H_s = 8.5 m | 4.23 s | 21.1 s | 12.7 s |

In a developed sea the shipped horizon was 21 s, a substantial fraction of the
360 s crossfade it had to follow, so `r_S` spent most of the transition on a
stale operating point. That is the lag the sweep is measuring.

## 5. A discrepancy-gated horizon does not help

The two instruments want opposite things - a stationary sea wants a long
horizon for noise rejection, a moving sea wants a short one - which suggests
making the horizon depend on how far the target has drifted from the applied
value rather than picking one constant. `seastate::common::adaptiveSmoothingHorizonSec`
implements that: with a threshold `slew_log > 0` the horizon is divided by
`1 + (|ln(target/applied)| / slew_log)^2`, so jitter inside the threshold leaves
the long horizon intact while a sustained move shortens it quadratically.

Measured on the `large` pair at the shipped multiplier, crossfade interval,
one IMU seed triplet (n = 10):

| `slew_log` | 3D RMS | roll/pitch | accel bias |
| --- | --- | --- | --- |
| 0.1 | 0.9860 | 1.0726 | 1.1214 |
| 0.2 | 0.9900 | 1.0628 | 1.1056 |
| 0.35 | 0.9940 | 1.0433 | 1.0796 |
| 0.75 | 0.9980 | 0.9924 | 1.0052 |
| 1.5 | 0.9993 | 1.0090 | 1.0184 |
| off (deployed) | 1.0000 | 1.0000 | 1.0000 |

It reaches the same displacement as a short fixed horizon and pays far more for
it: `slew_log = 0.1` buys 3D RMS 0.986, which `ADAPT_RS_MULT = 1.5` also buys,
but at 7.3% worse roll/pitch and 12.1% worse accelerometer bias instead of
better. The reason is that the gate keys on the *instantaneous* discrepancy,
which on a stationary sea is dominated by tuner jitter rather than by any real
move, so the horizon collapses on every excursion. Distinguishing a real move
from jitter needs persistence, not magnitude.

The mechanism is kept, disabled (`ADAPT_RS_SLEW_LOG = 0`,
`ADAPT_R_SLEW_LOG = 0`) and reachable through
`OU_ADAPT_RS_SLEW_LOG` / `OU_ADAPT_R_SLEW_LOG`, so the claim stays
reproducible. It is not a recommended configuration.

## 6. Adopted values

`ADAPT_RS_MULT`, `ADAPT_R_p0_MULT` and `ADAPT_R_v0_MULT` all move from `5.0` to
`3.0`. Ratios to the shipped configuration, pooled by geometric mean:

| | OU-III stationary | OU-III blend | OU-II stationary | OU-II blend |
| --- | --- | --- | --- | --- |
| n | 24 | 30 | 24 | 30 |
| 3D RMS | 0.9999 | **0.9927** | 0.9989 | **0.9957** |
| Z RMS | 0.9997 | **0.9916** | 0.9998 | **0.9944** |
| Z RMS, worst record | 1.0066 | 1.0030 | 1.0073 | 1.0022 |
| roll/pitch | 1.0039 | 0.9976 | 0.9974 | 0.9885 |
| yaw | 1.0000 | 0.9995 | 0.9991 | 1.0012 |
| accel bias | 1.0033 | 0.9980 | 0.9990 | 0.9841 |

Over the whole 900 s transition window - the convention the committed
validation scores, three quarters of which is stationary - the same change is
worth 0.31% of 3D RMS for OU-III and 0.21% for OU-II.

The shorter horizons were not taken further than 3.0 even though the crossfade
score keeps improving down to about 1.5, because the stationary worst-record
vertical error rises monotonically the other way and the accelerometer-bias
sentinel moves with it. 3.0 is where the two curves cross at a cost worth
paying.

### Sentinel re-derivation

`bias_3d_percent` is the only gate in either simulator that the change moves
past its limit; every other sentinel still passes on its existing value.

| family | before | after | binding record |
| --- | --- | --- | --- |
| OU-III | 106.25 (limit 106.8) | 108.88 (limit 109.4) | JONSWAP H_s = 1.5 m, accel |
| OU-II | 77.38 (limit 77.8) | 81.75 (limit 82.3) | JONSWAP H_s = 1.5 m, accel |

These sentinels carry about half a percent of headroom by construction, so any
change to the operating point moves them; even `ADAPT_RS_MULT = 4.0` exceeds
the old OU-III limit (107.35). What the number reports is the horizontal
accelerometer bias on the H_s = 1.5 m record, where the estimate is already
worse than predicting zero - 130% of the maximum true bias for OU-II and 182%
for OU-III before the change. The shorter horizon lets the pseudo-measurements
absorb low-frequency content the bias state was absorbing, which is why the
displacement error on that same record is unchanged (OU-II 20.779% -> 20.776%)
while the bias aggregate rises. It is a redistribution between two states, not
a loss of displacement accuracy.

## 7. Second pass: the same question on a crossfade fast enough to lag

Sections 3-6 were measured against a 360 s crossfade. That is a long ramp
compared with the horizon being tuned: at the deployed operating points
`3.0 * tau` is 4-13 s, so the crossfade was 27 to 90 smoothing horizons long and
the smoother had time to settle at every point along it. What the blend
interval mostly scored was therefore how well the schedule sat on a slowly
moving target, not how fast it followed. `tools/ou_validation.py` now crossfades over 120 s, which is about
three times the longest averaging memory anywhere in the adaptation path, and
the endpoint sea is scored twice -- a run-on interval one crossfade long, then
the settled remainder -- so a horizon that keeps carrying the sea it just left
is charged for it separately from one that is simply wrong.

Sweeping `ADAPT_RS_MULT` on that record. Ratios to the shipped `3.0` on the
identical realization; transition columns pool 2 endpoint pairs x 2 directions
x 3 wave-phase seeds x 3 IMU seed triplets (n = 36), the stationary column the
four JONSWAP and four PM-Stokes records at the same three seed triplets
(n = 24). Vertical RMS throughout:

| `ADAPT_RS_MULT` | blend | blend worst | run-on | settled | whole window | stationary | stationary worst |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.35 | 0.9761 | 1.0627 | 0.9842 | 0.9900 | 0.9865 | 0.9907 | 0.9962 |
| 0.5 | 0.9761 | 1.0550 | 0.9859 | 0.9910 | 0.9875 | 0.9923 | 0.9978 |
| 0.75 | 0.9778 | 1.0413 | 0.9885 | 0.9927 | 0.9895 | 0.9944 | 0.9994 |
| 1.0 | 0.9801 | 1.0300 | 0.9905 | 0.9942 | 0.9913 | 0.9958 | 1.0003 |
| **1.5** | **0.9850** | **1.0156** | **0.9937** | **0.9964** | **0.9943** | **0.9977** | **1.0013** |
| 2.0 | 0.9900 | 1.0077 | 0.9961 | 0.9980 | 0.9966 | 0.9988 | 1.0013 |
| 3.0 (was shipped) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 4.0 | 1.0094 | 1.0323 | 1.0033 | 1.0013 | 1.0026 | 1.0005 | 1.0032 |
| 5.0 | 1.0180 | 1.0618 | 1.0066 | 1.0022 | 1.0047 | 1.0007 | 1.0056 |
| 8.0 | 1.0394 | 1.1341 | 1.0177 | 1.0036 | 1.0095 | 1.0007 | 1.0098 |

Three things changed relative to the first pass, and only the first of them is
about the constant:

- **The run-on cost is now scored on its own.** The run-on interval moves
  monotonically with the horizon (0.9905 at 1.0, 1.0177 at 8.0), in the same
  direction as the blend and about five times harder than the settled interval
  beside it (1.0036 at 8.0). Under the 360 s crossfade that interval was the
  first quarter of a single 420 s endpoint segment, so the same cost was
  averaged down by the settled sea around it.
- **The stationary ensemble now agrees with the transition ensemble.** In the
  first pass the stationary worst-record vertical error rose as the horizon
  shortened, and that is what stopped the constant at 3.0. It no longer does:
  at 1.5 the stationary mean is 0.9977 and the worst record 1.0013, and the
  curve is still falling at 0.35. What changed in between is the rest of the
  filter -- the operating point moved to the wave band, the `r_S` law was
  re-gauged, and `MIN_R_S` came down -- so this is a re-measurement at a
  different operating point rather than a contradiction of the old one.
- **The cost is in the tail, not the mean.** The mean keeps improving down to
  0.35 while the worst single realization degrades monotonically the other way:
  +0.8% at 2.0, +1.6% at 1.5, +3.0% at 1.0, +6.3% at 0.35, all on the blend.
  That is the smoother passing tuner jitter into `r_S` on the record where the
  jitter happens to be worst.

`1.5` is where the mean gain and the worst-record loss cross on the blend:
-1.50% against +1.56%. That is the same shape of rule section 6 stopped at
`3.0` with -- take the shortest horizon whose worst-record cost has not yet
outrun its mean gain -- applied to the instrument that can now see the cost.
It is not the minimum of the mean curve, and it is not meant to be: below it
every further tenth of a percent of average gain is bought with two to three
tenths of worst-case loss.

Attitude and the accelerometer bias do not object. On the transition ensemble
roll/pitch is 0.9992 [0.9984, 1.0001] and the bias 0.9999 [0.9991, 1.0006] in
the blend; on the stationary ensemble 0.9998 [0.9996, 1.0000] and 0.9998
[0.9996, 1.0000]. The first pass's finding that these degrade below about 1.5
was measured with the tuner jitter of the pre-wave-band operating point; at the
deployed one they are flat over the whole grid.

### The other two averaging horizons

The `r_S` channel is not the only averaging horizon in the adaptation path, so
the same instrument was pointed at the other two. Neither moves.

`OU_ADAPT_TAU_SEC`, the fixed-seconds EMA on applied `tau` and `sigma_aw`
(shipped 1.8 s), vertical RMS ratios on the same ensemble:

| `ADAPT_TAU_SEC` | blend | run-on | settled | whole window |
| --- | --- | --- | --- | --- |
| 0.5 | 0.9983 | 0.9992 | 0.9993 | 0.9991 |
| 1.0 | 0.9988 | 0.9995 | 0.9996 | 0.9994 |
| 1.8 (shipped) | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 3.0 | 1.0017 | 1.0007 | 1.0005 | 1.0008 |
| 6.0 | 1.0051 | 1.0019 | 1.0011 | 1.0019 |
| 12.0 | 1.0102 | 1.0038 | 1.0015 | 1.0030 |

The sign is the same as `r_S`'s -- shorter follows a transition better -- but
the whole grid spans 0.4% against `r_S`'s 4%, and the available gain below the
shipped value is 0.06% of the window. `tau` enters `r_S` as `tau^3`, so
shortening this horizon is the expensive way to buy what `ADAPT_RS_MULT` buys
directly. It stays at 1.8 s.

`OU_SIGMA_VAR_K_PERIODS`, the `sigma_a` variance horizon, is measured in
docs/ou-sigma-horizon.md and stays at 2.0 for a different reason: its effect
reverses sign between a building and a decaying sea, so there is no
direction-symmetric optimum to move it to.

### Sentinel re-derivation, second pass

All nine OU-III deterministic gates were re-cut with `tools/ou_regauge_gates.py`
against the deployed filter, as docs/quality-gate-regauge.md requires. Eight of
the nine move down, i.e. the filter got better on their binding records; pitch
is the only gate the old bar no longer held, at 0.1940 -> 0.1952 deg on
PM-Stokes H_s = 4.0 m (+0.6%, against a bar cut with half a percent of margin),
which is a re-draw of the same kind section 8 of docs/ou-sigma-horizon.md
describes -- roll/pitch improves on average over both ensembles while this one
record moves inside its own scatter. The full before/after table is in the
comment above `FAIL_LIMITS` in `tests/kalman_ou_iii/kalman_ou_iii-sim.cpp`.

## 8. Reproducing

```sh
# stationary sweep, one knob at a time
tools/ou_ema_adapt_study.py stage1 --family OU_III --records jonswap,pmstokes

# crossfade sweep, both transition directions, five wave-phase seeds
tools/ou_ema_adapt_study.py transition --family OU_III \
    --transition-pair large --transition-seeds 11,12,13,14,15

# section 7: the whole grid on both endpoint pairs and both directions
tools/ou_ema_adapt_study.py transition --family OU_III \
    --grid 0.35,0.5,0.75,1,1.5,2,4,5,8 \
    --transition-pair large,small --transition-dir up,down \
    --transition-seeds 11,12,13 --seeds default,7,99

# the same grid on the stationary ensemble, which is what bounds it from below
tools/ou_ema_adapt_study.py stage1 --family OU_III \
    --grid 0.35,0.5,0.75,1,1.5,2,4,5,8 \
    --records jonswap,pmstokes --seeds default,7,99

# the fixed-seconds tau/sigma horizon, on the same records
tools/ou_ema_adapt_study.py transition --family OU_III --transition-stage confirm \
    --point 'atau=1:OU_ADAPT_TAU_SEC=1' --point 'atau=6:OU_ADAPT_TAU_SEC=6' \
    --transition-pair large,small --transition-dir up,down \
    --transition-seeds 11,12,13 --seeds default,7,99

# scaling-law check: same sweep one octave down in wave period
tools/ou_ema_adapt_study.py transition --family OU_III --transition-pair small

# joint (r_p0, r_v0) grid for OU-II
tools/ou_ema_adapt_study.py transition --family OU_II --transition-stage grid2d

# a named point across three IMU seed triplets
tools/ou_ema_adapt_study.py transition --family OU_II --transition-stage confirm \
    --point 'p0=3,v0=3:OU_ADAPT_R_P0_MULT=3:OU_ADAPT_R_V0_MULT=3' \
    --seeds default,7,99
```

The multipliers are also reachable at runtime through
`setRSAdaptMult()` (OU-III) and `setR_p0_AdaptMult()` / `setR_v0_AdaptMult()`
(OU-II), and are part of the `tools/ou_tuning_sweep.py` search space.
