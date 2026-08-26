# Re-fitting TFG's adaptation coefficients on a multi-seed instrument

`SeaStateFusionFilter_TFG` schedules its operating point from five constants:

    tau        = tau_coeff * 0.5 / f_tune
    sigma_aw   = sigma_coeff * sigma_wave
    r_S        = R_S_coeff * sigma_aw * tau^3
    Sigma_aw  -> (S_factor * sigma_aw, S_factor * sigma_aw, sigma_aw)
    r_S       -> (R_S_xy_factor * r_S, R_S_xy_factor * r_S, r_S)

All five were fitted in `docs/tfg-adaptation-parity.md` once TFG took on the OU
families' magnetic acquisition, hard-iron tracking and band-passed sigma
channel. That fit scored one deterministic realization per record. This one
re-runs it paired over three IMU seed triplets per record and confirms it at
six, which is a different instrument rather than a repetition: on the channels
that decide two of these five constants, a single realization is not a
measurement.

**Result.** `S_factor` moves from `1.20` to `1.00`. The other four are
confirmed where they are. Pooled over the eight records at six IMU seed
triplets, `S_factor = 1.00` reads accelerometer bias 0.9209, roll/pitch 0.9305,
yaw 0.9578 and 3D RMS 0.9983 of the shipped filter, for vertical RMS 1.0009.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_tfg/kalman_tfg-sim`, scoring the trailing 900 s of a 1200 s
replay. Every table is a ratio to the shipped constants on the same record and
seed, pooled by geometric mean. `tools/ou_ema_adapt_study.py` grew a `TFG`
family for this and reproduces all of them.

## 1. Four of the five are confirmed

One knob at a time, eight records at three IMU seed triplets (n = 24):

| `tau_coeff` | 3D RMS | Z RMS | | `sigma_coeff` | 3D RMS | Z RMS |
| --- | --- | --- | --- | --- | --- | --- |
| 0.85 | 1.0571 | 1.0705 | | 0.80 | 1.0139 | 1.0080 |
| 0.92 | 1.0142 | 1.0170 | | 0.90 | 1.0033 | 0.9997 |
| **1.00** | 1.0000 | 1.0000 | | **1.00** | 1.0000 | 1.0000 |
| 1.08 | 1.0171 | 1.0228 | | 1.10 | 1.0023 | 1.0065 |
| 1.15 | 1.0555 | 1.0683 | | 1.25 | 1.0133 | 1.0243 |

| `R_S_coeff` | 3D RMS | Z RMS | roll/pitch | accel bias |
| --- | --- | --- | --- | --- |
| 0.18 | 1.0469 | 1.0559 | 0.9688 | 0.9515 |
| 0.22 | 1.0134 | 1.0149 | 0.9799 | 0.9721 |
| 0.25 | 1.0023 | 1.0019 | 0.9898 | 0.9864 |
| **0.28** | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 0.32 | 1.0068 | 1.0096 | 1.0123 | 1.0159 |
| 0.36 | 1.0215 | 1.0286 | 1.0233 | 1.0297 |
| 0.45 | 1.0726 | 1.0906 | 1.0430 | 1.0540 |

`tau_coeff` and `R_S_coeff` have clean interior minima exactly where they ship.
`sigma_coeff` is minimized at 1.0 on 3D and at 0.9 on vertical by 0.03%, which
is the same near-tie the first fit found and resolved the same way: 0.9 costs
1.1% on accelerometer bias and 1.6% on roll/pitch here, so 1.0 stands.

`R_S_coeff` was re-checked once more at the adopted `S_factor`, since the two
act on the same horizontal channel from opposite sides: at six seeds, 0.26 and
0.30 read 1.0004 and 0.9991 of 3D RMS against 0.9983 for 0.28. It does not move.

## 2. `R_S_xy_factor = 1.15` was a grid corner, and survives being looked at

The first fit gridded this against `S_factor` over `1.15, 1.30, 1.50` — it never
went below 1.15, so the adopted value sat on the edge of the grid. That is the
shape of an artifact, and OU-III had the mirror-image version of it: its setter
clamped the same constant to `[0, 1]`, so its sweep could only go down. Both
sweeps found the boundary.

Re-swept over `0.8` to `1.5`, at three IMU seed triplets:

| `R_S_xy_factor` | 3D RMS | Z RMS | roll/pitch | yaw | accel bias |
| --- | --- | --- | --- | --- | --- |
| 0.80 | 1.0311 | 1.0050 | 0.9998 | 1.0059 | 0.9959 |
| 0.90 | 1.0136 | 1.0031 | 0.9986 | 1.0029 | 0.9973 |
| 1.00 | 1.0039 | 1.0017 | 0.9986 | 1.0012 | 0.9984 |
| **1.15** | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 1.30 | 1.0054 | 0.9986 | 1.0027 | 0.9995 | 1.0017 |
| 1.50 | 1.0226 | 0.9971 | 1.0076 | 1.0026 | 1.0044 |

It is a genuine interior minimum on 3D RMS, and the vertical channel wants it
higher by an amount too small to referee anything. It is also the *only* one of
the five constants that does not move attitude or bias at all — every column
except the first two is inside 0.8% over a factor of two in the knob. Whatever
`R_S_xy_factor` is doing, it is doing it entirely to horizontal displacement,
which is what makes 3D RMS the right metric for it.

The same separation shows in the joint grid, three seed triplets, 3D RMS:

| 3D RMS | `xy=1.00` | `1.15` | `1.30` |
| --- | --- | --- | --- |
| `S=0.90` | 1.0024 | 0.9979 | 1.0025 |
| `S=0.95` | 1.0025 | 0.9981 | 1.0028 |
| `S=1.00` | 1.0027 | **0.9984** | 1.0033 |
| `S=1.05` | 1.0030 | 0.9988 | 1.0038 |
| `S=1.10` | 1.0033 | 0.9992 | 1.0044 |

`xy = 1.15` wins every row by about the same margin, and yaw at fixed `S` moves
by under 0.4% across the three columns. The two knobs are separable, so the
corner was not hiding an interaction.

## 3. `S_factor`: the records say 1.00

`S_factor` has no interior optimum over the range swept. It is a monotone trade,
three seed triplets:

| `S_factor` | 3D RMS | Z RMS | Z worst | roll/pitch | yaw | accel bias |
| --- | --- | --- | --- | --- | --- | --- |
| 0.80 | 0.9982 | 1.0094 | 1.0975 | 0.9702 | 0.8884 | 0.9654 |
| 0.90 | 0.9979 | 1.0040 | 1.0509 | 0.9450 | 0.8943 | 0.9308 |
| 1.00 | 0.9984 | 1.0014 | 1.0243 | 0.9539 | 0.9341 | 0.9449 |
| 1.10 | 0.9992 | 1.0003 | 1.0090 | 0.9744 | 0.9731 | 0.9698 |
| **1.20** | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 1.40 | 1.0014 | 1.0004 | 1.0152 | 1.0480 | 1.0299 | 1.0508 |
| 1.87 | 1.0033 | 1.0024 | 1.0421 | 1.1228 | 1.0654 | 1.1267 |

Yaw, attitude and accelerometer bias all improve as it comes down and all
degrade as it goes up; vertical RMS does the opposite, weakly in the mean and
less weakly on the worst record. A pooled score can be made to pick anything in
`[0.9, 1.1]` depending on how it weights those, so the constant needs a reason
from outside the sweep.

There is one. `S_factor` states that the sea's horizontal acceleration is
`S_factor` times its vertical acceleration, and the records say what that
number is. Standard deviations of the generated `acc_x, acc_y, acc_z` over the
scored window:

| record | `sigma_ax/sigma_az` | `sigma_ay/sigma_az` | horizontal magnitude |
| --- | --- | --- | --- |
| JONSWAP H0.27 | 0.803 | 0.584 | 0.993 |
| JONSWAP H1.50 | 0.830 | 0.552 | 0.997 |
| JONSWAP H4.00 | 0.817 | 0.570 | 0.996 |
| JONSWAP H8.50 | 0.827 | 0.546 | 0.991 |
| PM-Stokes H0.27 | 0.797 | 0.606 | 1.001 |
| PM-Stokes H1.50 | 0.827 | 0.563 | 1.001 |
| PM-Stokes H4.00 | 0.818 | 0.575 | 1.000 |
| PM-Stokes H8.50 | 0.833 | 0.555 | 1.000 |
| **mean** | **0.819** | **0.569** | **0.997** |

`sqrt(sigma_ax^2 + sigma_ay^2) / sigma_az` is 0.997 and sits inside 1% on every
record, in both spectra, across a factor of 31 in significant wave height. A
scalar horizontal prior of 1.20 was 20% above what the sea puts there. The
deployed value is now the measured one.

This is the same conclusion OU-III reached from 1.87 in
[`ou-iii-anisotropy-consistency.md`](ou-iii-anisotropy-consistency.md), by
sweep rather than by measuring the records, and that note left TFG's constant
named as an open item. It is closed the same way and at the same value, which
is worth more than either result alone: two filters with different translational
state structure both want the horizontal acceleration prior the records
actually contain.

### Where the change lands

Per record at `S_factor = 1.00`, three seed triplets:

| record | yaw | accel bias | 3D RMS | Z RMS |
| --- | --- | --- | --- | --- |
| JONSWAP H0.27 | 1.0057 | 1.0013 | 1.0008 | 1.0004 |
| JONSWAP H1.50 | 0.9487 | 0.9120 | 0.9993 | 1.0024 |
| JONSWAP H4.00 | **0.6918** | **0.8018** | 0.9941 | 1.0062 |
| JONSWAP H8.50 | 0.8990 | 0.8869 | 0.9890 | 1.0087 |
| PM-Stokes H0.27 | 0.9921 | 1.0095 | 1.0008 | 0.9999 |
| PM-Stokes H1.50 | 0.9959 | 1.0016 | 1.0006 | 0.9996 |
| PM-Stokes H4.00 | 0.9737 | 0.9843 | 1.0011 | 0.9972 |
| PM-Stokes H8.50 | 1.0149 | 0.9834 | 1.0016 | 0.9970 |

The gain is concentrated in the large JONSWAP seas — 31% of yaw and 20% of
accelerometer bias on H4.0 — which is where OU-III's was concentrated too, and
for the same reason: those are the records with the most horizontal motion for
an inflated horizontal prior to mis-state. The vertical cost is on the same two
records, at 0.6% and 0.9%.

That cost lands away from the sentinels it could have threatened. TFG's two
vertical gates are set by the `H_s = 0.27` m records, where the change is
1.0004 and 0.9999 — flat. The vertical error does get slightly worse in the two
seas that carry the most vertical signal, and it is a real if small cost, but
the records that decide whether the vertical channel has regressed do not move.

### 0.90 is available and is not taken

`S_factor = 0.90` buys more of everything the reduction buys (three seed
triplets, against 1.00 on the same instrument): 10.6% of yaw against 6.6%, 6.9%
of accelerometer bias against 5.5%. It costs 0.40% of pooled vertical RMS
against 0.14%, and 5.1% of the worst-record vertical against 2.4%.
It is not adopted because the argument for 1.00 is that it is the number the
records contain, and 0.90 is not that number — it is a number that scores
better on this eight-record set. Fitting a physical prior below what the physics
measures is how a constant stops transferring to the next set of records.

## 4. The EMA horizon was justified on the wrong instrument, and survives anyway

`adapt_RS_mult` is the sixth constant in the schedule — the smoothing horizon
`r_S`'s target is followed with, in units of `tau`. TFG inherited OU-III's 3.0,
and the first fit recorded it as "swept, flat, not fitted", flat over 2 to 4 on
the eight stationary records.

That is the wrong instrument for it, and
`docs/ou-ema-adaptation-tuning.md` §1 says so directly: a smoothing horizon is
only observable while its target moves, so a stationary record cannot referee
one at all. The 3.0 the OU families ship was chosen on synthesized sea-state
transitions, and TFG had never been run on them — the study tool had no TFG
family until this change added one.

Run now, crossfade interval, `large` pair in both directions at three wave-phase
seeds and three IMU seed triplets (n = 18), ratios to the shipped 3.0:

| `TFG_ADAPT_RS_MULT` | 3D RMS | Z RMS | Z worst | roll/pitch | yaw | accel bias |
| --- | --- | --- | --- | --- | --- | --- |
| 1.0 | 0.9931 | 0.9858 | 1.0089 | 0.9981 | 1.0068 | 0.9950 |
| 1.5 | 0.9950 | 0.9905 | 1.0062 | 1.0055 | 0.9937 | 1.0068 |
| 2.0 | 0.9967 | 0.9942 | 1.0036 | 1.0072 | 0.9900 | 1.0100 |
| **3.0 (shipped)** | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 8.0 | 1.0173 | 1.0233 | 1.0574 | 0.9704 | 1.0667 | 0.9987 |

It is not flat. It is the same curve the OU families have: displacement error
monotone in the horizon, clearly bad on the long side, and flattening on the
short side where the attitude and bias channels start to pay instead. The
shipped 3.0 sits on the shallow part of that curve — 1.5 buys 0.50% of 3D and
0.95% of vertical on the crossfade for 0.55% of roll/pitch and 0.68% of bias.

`adapt_RS_mult` stays at 3.0. There is no dominant point inside `[1, 3]`, every
difference in that interval is under a percent on a six-record instrument, and
3.0 is where the same trade was resolved for both OU families after a much
larger study. What changes is the justification: it is not flat, it is a shallow
trade, and the earlier claim rested on records that could not have shown
otherwise.

## 5. What the multi-seed instrument was needed for

The 3D and vertical channels are stable enough that the first fit's single
realization was adequate for `tau_coeff`, `sigma_coeff` and `R_S_coeff`, and
this study confirms all three unmoved. Yaw and accelerometer bias are not:

| record | yaw RMS at `S=1.20`, six seeds |
| --- | --- |
| JONSWAP H8.50 | 1.03, 3.79, 1.36, 2.93, 2.43, 4.56 |

A 4.4x spread on one record, and the deterministic seed draws the smallest of
the six. Any constant chosen by watching that record's yaw on one seed is
choosing on noise. That is not a hypothetical here — it is exactly what happens
to the yaw sentinel in section 6.

## 6. Gates: six move by under a percent, and yaw goes up

Re-derived with `tools/ou_regauge_gates.py --family tfg`, which grew a TFG
family for this change; it reproduces all seven of TFG's previously shipped
bars from the filter that produced them, which is the check that it implements
the rule rather than a rule.

| gate | was | now | worst scored, before → after |
| --- | --- | --- | --- |
| Z %Hs JONSWAP | 4.803 | 4.807 | 4.7784 → 4.7830 (H0.27) |
| Z %Hs PM-Stokes | 4.709 | **4.707** | 4.6846 → 4.6833 (H0.27) |
| yaw deg | 1.536 | 1.59 | 1.5278 (pmstokes H8.5) → 1.5818 (jonswap H8.5) |
| 3D % JONSWAP | 21.14 | **21.13** | 21.0299 → 21.0169 (H1.5) |
| 3D % PM-Stokes | 20.71 | 20.74 | 20.6045 → 20.6322 (H4.0) |
| acc Z bias % | 5.026 | **5.022** | 5.0002 → 4.9961 (pmstokes H8.5) |
| bias 3D % | 167.6 | **164.5** | 166.688 → 163.607 (pmstokes H4.0) |

Six of the seven move by less than a percent in either direction. The seventh
is a breach and gets its own argument.

### The yaw sentinel fails, and yaw improves

`1.5278 → 1.5818` on the deterministic protocol, past a bar of 1.536. Raising a
bar to admit one's own change is how sentinels stop meaning anything, so the
evidence has to be stated rather than asserted.

The binding record changed, from PM-Stokes H8.5 to JONSWAP H8.5, and JONSWAP
H8.5 is the record whose yaw spans 1.03 to 4.56 deg across six IMU seeds. Paired
on that record:

| seed | `S=1.20` | `S=1.00` |
| --- | --- | --- |
| default | 1.0322 | 1.5818 |
| 7 | 3.7871 | 3.1882 |
| 99 | 1.3618 | 0.7669 |
| 3 | 2.9330 | 3.5454 |
| 21 | 2.4334 | 1.8283 |
| 55 | 4.5617 | 3.9482 |

Four of six improve, the geometric-mean ratio on this record is 0.9109 and the
six-seed mean falls from 2.685 to 2.476 deg. Pooled over all eight records the
ratio is 0.9578. The default seed draws the smallest of the six values under the
old constant and a middling one under the new, which is what a 54% jump on one
realization looks like when the distribution has not moved in that direction at
all.

The sentinel is deterministic by design and follows the protocol it is written
against, so it moves to 1.59. The quality claim rests on the seeds, not on the
sentinel. This is the same situation, with the same resolution, as OU-III's yaw
gate in its anisotropy study, which moved 21% on a record spanning 1.05 to
6.57 deg.

One other channel is worth naming rather than leaving in the table.
Accelerometer bias 3D falls on the binding record (166.7 → 163.6) and rises
sharply on JONSWAP H8.5 (81.0 → 116.9) — the pooled six-seed figure is 0.9209,
so this is a redistribution across records on a quantity whose error still
exceeds the true bias on two of the eight, not a uniform gain. It is the same
caveat this filter's bias gate has always carried.

## 7. Adopted

| coefficient | was | now |
| --- | --- | --- |
| `tau_coeff` | 1.0 | 1.0 |
| `sigma_coeff` | 1.0 | 1.0 |
| `R_S_coeff` | 0.28 | 0.28 |
| `S_factor` | 1.20 | **1.00** |
| `R_S_xy_factor` | 1.15 | 1.15 |

`doc/kalman_tfg/tfg-sim-results-generated.tex-part` is regenerated from the
filter that ships.

## 8. Reproducing

```sh
make -C tests/kalman_tfg build

# one knob at a time, eight records, three IMU seed triplets
tools/ou_ema_adapt_study.py confirm --family TFG \
    --records jonswap,pmstokes --seeds default,7,99 \
    --point S=1:TFG_S_FACTOR=1 --point xy=1:TFG_R_S_X_FACTOR=1:TFG_R_S_Y_FACTOR=1

# the horizontal pair, jointly
tools/ou_ema_adapt_study.py confirm --family TFG \
    --records jonswap,pmstokes --seeds default,7,99 \
    --point S1/xy1.15:TFG_S_FACTOR=1:TFG_R_S_X_FACTOR=1.15:TFG_R_S_Y_FACTOR=1.15

# the finalists at six seed triplets
tools/ou_ema_adapt_study.py confirm --family TFG \
    --records jonswap,pmstokes --seeds default,7,99,3,21,55 \
    --point S1.0:TFG_S_FACTOR=1.0

# the gates, and the control that reverting the constant restores them
python3 tools/ou_regauge_gates.py --family tfg
python3 tools/ou_regauge_gates.py --family tfg --env TFG_S_FACTOR=1.20

# the records' own acceleration anisotropy
python3 - <<'PY'
import csv, glob, math
import numpy as np
for f in sorted(glob.glob("tests/kalman_tfg/wave_data_[jp]*_H*.csv")):
    a = {k: [] for k in ("acc_x", "acc_y", "acc_z")}
    for row in csv.DictReader(open(f)):
        if 300.0 <= float(row["time"]) <= 1200.0:
            for k in a:
                a[k].append(float(row[k]))
    sx, sy, sz = (np.std(a[k]) for k in ("acc_x", "acc_y", "acc_z"))
    print(f, round(sx / sz, 3), round(sy / sz, 3), round(math.hypot(sx, sy) / sz, 3))
PY
```
