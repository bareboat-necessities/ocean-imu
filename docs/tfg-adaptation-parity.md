# Bringing TFG to OU parity: adaptation, magnetic init, hard iron — and the retune

`SeaStateFusionFilter_TFG` was written as the core of `SeaStateFusionFilter_OU_III`
and then drifted from it. This records what was brought across, what each piece
is worth on its own, and how the tuner coefficients were re-fitted afterwards.

Everything below is the eight-record set the quality gates use — four JONSWAP
and four PM-Stokes significant wave heights — scored over the final 900 s of a
twenty-minute replay with the default sensor seeds. Columns are: vertical RMS as
a percentage of `H_s`, 3D RMS as a percentage of peak 3D displacement, attitude
RMS in degrees, and accelerometer-bias RMS as a percentage of the true bias.
"worst" is the worst single record, which is what the gates are cut against.

Reproduce any row with:

```
make -C tests/kalman_tfg build
cd tests/kalman_tfg
W3D_COLLECT_ALL_GATES=1 W3D_WRITE_TIMESERIES=0 <VARS> ./kalman_tfg-sim
```

## 1. Headline

| channel | before | after | change |
|---|---|---|---|
| vertical RMS %H_s, mean | 4.638 | 4.369 | **−5.8%** |
| vertical RMS %H_s, worst | 5.209 | 4.778 | **−8.3%** |
| 3D RMS %, mean | 21.04 | 19.74 | **−6.2%** |
| 3D RMS %, worst | 25.78 | 21.03 | **−18.4%** |
| yaw RMS deg, mean | 2.149 | 1.077 | **−49.9%** |
| yaw RMS deg, worst | 2.923 | 1.528 | **−47.7%** |
| roll RMS deg, mean | 0.698 | 0.308 | **−55.9%** |
| pitch RMS deg, mean | 0.306 | 0.340 | **+11.3%** |
| acc-bias RMS %, mean | 144.1 | 92.7 | **−35.7%** |
| acc-bias RMS %, worst | 398.2 | 166.7 | **−58.1%** |

Per record:

| record | Z before | Z after | 3D before | 3D after | yaw before | yaw after | bias before | bias after |
|---|---|---|---|---|---|---|---|---|
| JONSWAP 0.27 | 5.209 | 4.778 | 20.64 | 19.96 | 2.122 | 1.465 | 134.1 | 79.1 |
| JONSWAP 1.50 | 4.641 | 4.307 | 20.99 | 21.03 | 2.154 | 1.151 | 88.9 | 43.9 |
| JONSWAP 4.00 | 4.781 | 4.239 | 19.67 | 19.89 | 1.793 | 0.704 | 398.2 | 84.5 |
| JONSWAP 8.50 | 3.940 | 3.980 | 17.25 | 17.73 | 1.712 | 1.032 | 49.6 | 81.0 |
| PM-Stokes 0.27 | 5.100 | 4.685 | 21.19 | 20.36 | 2.320 | 0.938 | 88.4 | 68.5 |
| PM-Stokes 1.50 | 4.392 | 4.016 | 19.11 | 18.47 | 2.094 | 0.774 | 171.3 | 73.8 |
| PM-Stokes 4.00 | 4.447 | 4.452 | 25.78 | 20.60 | 2.923 | 1.024 | 86.3 | 166.7 |
| PM-Stokes 8.50 | 4.594 | 4.491 | 23.68 | 19.87 | 2.072 | 1.528 | 136.2 | 143.8 |

Yaw improves on 8 of 8 records, vertical on 6 of 8, 3D on 5 of 8 — and the two
records that dominated the old worst case (PM-Stokes 4.0 and 8.5 on 3D) come
down by 5.2 and 3.8 points. Accelerometer bias is the mixed one: four records
improve substantially, two get worse, and the worst case more than halves.

Pitch is the only channel whose *mean* moves the wrong way, and it does so only
at the default seed. Under `W3D_SEED=20260813` it goes 0.558° → 0.369° and under
`W3D_SEED=777` it goes 0.674° → 0.345°; under both alternative realizations every
channel in the headline table improves. The loss at the default seed is
concentrated on the two large PM-Stokes records and arrives with the magnetic
refinement, which halves roll on the same set — combined attitude RMS
`sqrt(roll² + pitch²)` goes 0.762° → 0.459°.

## 2. What each piece is worth

One feature turned off at a time, at the deployed coefficients.

| arm | Z mean | Z worst | 3D mean | 3D worst | yaw mean | yaw worst | roll | pitch | bias mean | bias worst |
|---|---|---|---|---|---|---|---|---|---|---|
| **deployed** | 4.388 | 4.787 | 18.82 | 20.32 | 0.777 | 1.345 | 0.296 | 0.366 | 97.1 | 164.1 |
| `TFG_WAVE_BAND=0` | 4.432 | 4.880 | 18.97 | 20.64 | 0.816 | 1.513 | 0.290 | 0.396 | 102.9 | 166.6 |
| `TFG_AW_COV_SYNC=0` | 4.389 | 4.793 | 18.83 | 20.33 | 0.772 | 1.310 | 0.298 | 0.357 | 95.9 | 163.8 |
| `TFG_MAG_REFINE=0` | 4.301 | 4.738 | 18.73 | 20.25 | 0.688 | 1.534 | 0.600 | 0.241 | 118.0 | 204.3 |
| `TFG_MAG_HARD_IRON=0` | 4.412 | 4.779 | 19.12 | 20.50 | **2.486** | **3.256** | 0.307 | 0.373 | 94.3 | 159.3 |
| `TFG_STARTUP_INIT=staged` | 4.431 | 4.814 | 19.84 | **23.42** | 2.265 | 2.906 | 0.697 | 0.303 | 145.7 | **398.3** |
| old coefficients | 4.447 | 5.006 | 19.89 | 21.18 | 1.097 | 1.759 | 0.329 | 0.390 | 99.5 | 168.5 |

The six live arms were re-measured on the current tree, after the hard-iron
ridge floor was re-cut; the yaw columns in particular are not the ones this
section first reported, and `doc/kalman_ou_iii/ins-startup.tex` quotes the
re-measured values. The `old coefficients` row keeps its original measurement:
that coefficient set no longer exists in the tree, so it cannot be replayed.

Reading it:

**Continuous hard iron is the largest single effect, and it is entirely in yaw.**
Turning it off costs 3.2× on mean yaw and 2.4× on worst yaw and moves nothing
else by more than 4%. This is the feature both OU families carried and this one
did not, and it is the whole of the previously-unexplained yaw gap.

**The proxy startup policy is what fixes the horizontal channel and the bias.**
Worst 3D goes 20.32 → 23.42 and worst bias 164.1 → 398.3 without it. Note that
the staged arm here still has the band-passed sigma channel and the retuned
coefficients, so it is not the old filter — it is the old *startup*, which is
worth about 15% of worst-case 3D and 143% of worst-case bias on its own.

**The magnetic refinement is a genuine trade.** It costs 2.0% of vertical and
13% of mean yaw, and buys a halving of roll (0.600 → 0.296), 18% of mean
accelerometer bias (118.0 → 97.1) and 12% of worst yaw (1.534 → 1.345). It stays on for the same reason OU-III keeps
it: the channels it improves are the weak ones, and roll and bias are the pair
that the provisional reference's tilt error was being parked in.

**The band-passed sigma channel** is worth 1–2% on vertical, 11% on worst yaw
and 6% on mean bias.

**The periodic a_w covariance sync** is within noise on everything: −0.1% on 3D
mean, +1.2% the wrong way on bias mean. It is kept for the consistency argument
in `docs/tfg-design.md` section 8, not for the number, and that is stated rather
than dressed up.

## 3. The retune

Every item above changes what the tuner coefficients multiply, so they were
re-fitted afterwards rather than inherited.

| coefficient | was | now |
|---|---|---|
| `tau_coeff` | 1.0 | 1.0 |
| `sigma_coeff` | 1.0 | 1.0 |
| `R_S_coeff` | 0.35 | **0.28** |
| `S_factor` | 1.87 | **1.20** |
| `R_S_xy_factor` | 1.0 | **1.15** |
| `adapt_tau_sec`, `adapt_RS_mult`, noise floor, band ratios | — | unchanged |

Every sweep in this section scores one deterministic realization per record.
That is adequate for the displacement channels and is not adequate for yaw or
accelerometer bias, which on the largest JONSWAP record span a factor of four
across IMU seeds. Re-run paired over three seed triplets,
[`docs/tfg-adaptation-refit.md`](tfg-adaptation-refit.md) confirms `tau_coeff`,
`sigma_coeff`, `R_S_coeff` and `R_S_xy_factor` where this section left them and
moves `S_factor` from 1.20 to **1.00**, which is the horizontal-to-vertical
acceleration ratio the records themselves carry. It also re-sweeps
`R_S_xy_factor` below 1.15, which the grid below never did.

### tau_coeff — clean minimum at 1.0

| `TFG_TAU_COEFF` | 0.85 | 0.92 | **1.00** | 1.10 | 1.20 |
|---|---|---|---|---|---|
| Z mean | 4.735 | 4.479 | **4.374** | 4.482 | 4.803 |
| 3D mean | 21.74 | 20.53 | **19.89** | 20.00 | 21.08 |

`tau = tau_coeff · 0.5 / f_tune` is the law, and the sea's own zero-crossing
period is what it is. A coefficient other than 1 would be saying the schedule's
own statement of the operating point is wrong.

### sigma_coeff and R_S_coeff — swept jointly, because `r_S ~ sigma·tau³`

Sweeping either alone moves both the OU stationary scale and the integral
regularizer, so they were gridded.

| `sigma_coeff` \\ `R_S_coeff` | 0.22 | 0.28 | 0.35 |
|---|---|---|---|
| 0.80 | 4.653 / 21.63 | 4.435 / 20.54 | 4.373 / 20.03 |
| 0.90 | 4.541 / 21.00 | 4.385 / 20.13 | 4.394 / 19.89 |
| **1.00** | 4.466 / 20.52 | **4.374 / 19.89** | 4.447 / 19.89 |

(Z mean / 3D mean.) The vertical minimum tracks the product `sigma_coeff ·
R_S_coeff ≈ 0.28`, which is the r_S the filter actually receives — but the two
are not interchangeable, because `sigma_coeff` also sets `Sigma_aw` directly.
Holding the product and moving `sigma_coeff` down to 0.8 costs 4% on 3D and 11%
on bias. So `sigma_coeff` stays at 1.0 and `R_S_coeff` takes the reduction.

OU-III runs `sigma_coeff = 0.9`. That difference is now genuine rather than a
units mismatch — both filters measure sigma through the same band — and 0.9 was
tried here: it buys 0.4% of vertical at 4% on 3D and 4% on bias.

The `R_S_coeff` minimum is flat over 0.26–0.30 (Z mean 4.384 / 4.374 / 4.381)
and 0.28 is where 3D stops improving. Why it had to come down at all: the
band-passed sigma channel reads lower than the broadband one it replaced, and
`r_S ~ sigma·tau³` inherits that, so 0.35 now over-regularizes the `S = 0`
constraint.

### S_factor and R_S_xy_factor — the horizontal pair

Both act on the horizontal channel: `S_factor` scales horizontal `sigma_aw` up,
`R_S_xy_factor` loosens the horizontal `S = 0` constraint.

| S / xy | 1.15 | 1.30 | 1.50 |
|---|---|---|---|
| 1.00 | 4.382 / 19.70 / 96.4 | 4.374 / 19.79 / 94.9 | 4.365 / 20.12 / 92.8 |
| **1.20** | **4.369 / 19.74 / 92.7** | 4.362 / 19.84 / 91.5 | 4.353 / 20.19 / 89.8 |
| 1.50 | 4.365 / 19.78 / 93.9 | 4.359 / 19.89 / 93.1 | 4.352 / 20.19 / 92.0 |

(Z mean / 3D mean / bias mean.) `R_S_xy_factor` above about 1.2 starts trading
3D for yaw and bias; 1.15 is where it is still free on all four channels.
`S_factor` at 1.20 takes 5% off bias and 0.05° off pitch for no vertical cost —
the horizontal channel no longer needs the extra `a_w` headroom now that it is
not absorbing a standing heading error.

### Four things deliberately left at OU-III's values

Swept, flat, not fitted:

| parameter | range swept | spread in Z mean |
|---|---|---|
| `adapt_tau_sec` | 1.0 – 3.0 | 4.372 – 4.377 |
| `adapt_RS_mult` | 2.0 – 4.0 | 4.372 – 4.374 |
| pre-band noise floor | 0.06 – 0.20 | 4.374 – 4.374 |
| sigma band ratios | five shapes, 0.35/4 to 0.5/6 | 4.371 – 4.404 |

The noise floor is flat because the band now *refers* it through the band's own
white-noise variance gain rather than subtracting it raw, so the schedule barely
notices a factor of three in the pre-band value. The band ratios have no shape
that beats 0.5/4.0 on more than one channel.

### And the refinement schedule, which is phase noise, not a trend

| `mag_refine_start_sec` | 45 | 60 | **90** | 150 |
|---|---|---|---|---|
| Z mean | 4.367 | 4.367 | 4.369 | 4.388 |
| bias worst | 166.7 | 166.7 | 166.7 | 166.7 |

| `mag_refine_window_sec` | 15 | **30** | 60 |
|---|---|---|---|
| Z mean | 4.333 | 4.369 | 4.308 |
| 3D worst | 21.46 | 21.03 | 21.06 |
| yaw worst | 1.421 | 1.528 | 1.607 |
| bias worst | 134.0 | 166.7 | 137.1 |

The window sweep is non-monotonic on three of four channels. That is the phase
sensitivity of a single one-shot re-acquisition landing at a particular point in
a particular wave, not a trend in window length, so fitting it to these eight
records would be fitting noise. The values stay at OU-III's.

## 4. Robustness off the tuning set

The coefficients were fitted on the default seed. Repeating the full before/after
comparison under two other sensor-noise and bias realizations:

| | default | `W3D_SEED=20260813` | `W3D_SEED=777` |
|---|---|---|---|
| Z mean | 4.638 → 4.369 | 4.789 → 4.107 | 4.731 → 4.305 |
| Z worst | 5.209 → 4.778 | 6.682 → 4.634 | 5.190 → 4.679 |
| 3D worst | 25.78 → 21.03 | 28.34 → 21.24 | 25.99 → 21.24 |
| yaw mean | 2.149 → 1.077 | 5.023 → 3.217 | 3.401 → 2.073 |
| roll mean | 0.698 → 0.308 | 0.864 → 0.355 | 0.769 → 0.217 |
| pitch mean | 0.306 → 0.340 | 0.558 → 0.369 | 0.674 → 0.345 |
| bias worst | 398.2 → 166.7 | 389.1 → 121.9 | 368.2 → 105.8 |

Every channel improves on both held-out realizations, pitch included. The gates
in `kalman_tfg-sim.cpp` are cut against the default seed only, which is the
protocol the whole repository uses; the two extra seeds are evidence that the
retune is not fitted to one realization, not a second gate.

There is no held-out *record* set: `W3dSimulationRunner` scores only JONSWAP and
PM-Stokes, so the other wave families in the data are not available as one.
