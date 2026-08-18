# Re-fitting OU-II's two pseudo-measurement variance coefficients

> **Superseded as the deployed law.** OU-II now schedules both pseudo channels
> from the joint physical-MSE derivation of
> `doc/kalman_ou_ii/ou2-dual-regularization-mse.tex`; see
> [`ou-ii-dual-mse-adaptation.md`](ou-ii-dual-mse-adaptation.md). The
> `(R_p0_coeff, R_v0_coeff)` pair fitted here is still the calibration of the
> `Empirical` law, which remains selectable, so this write-up stays the record
> of how that arm got its coefficients.

OU-II regularizes double-integration drift with two zero pseudo-measurements
fired on one periodic tick inside the MEKF, `p = 0` and `v = 0`, whose standard
deviations the wrapper schedules from the tuner's operating point:

    r_p0 = R_p0_coeff * sigma_aw * tau^2
    r_v0 = R_v0_coeff * sigma_aw * tau

OU-III's single analogue, `r_S = R_S_coeff * sigma_aw * tau^3`, was re-fitted
against the per-record optimum when its tuner moved into the wave band, and its
coefficient has been re-measured since. OU-II's pair had not been touched since
that same wave-band change, and one thing that entered the law afterwards was
never priced: the parity work moved the `sigma_aw` channel behind
`AdaptiveWaveBandPass`, and a band-passed variance channel reads lower than the
broadband one it replaced. `r = c * sigma_aw * tau^k` inherits that reduction
directly, so a coefficient fitted before the band leaves the estimator
over-regularized after it.

**Result.** `R_v0_coeff` moves from `1.1` to `1.3` and `R_p0_coeff` from `0.6`
to `0.65`. Pooled over the eight stationary records at six IMU seed triplets,
that is 3D RMS 0.9904, vertical RMS 0.9958, roll/pitch 0.9929 and yaw 0.9970 of
the shipped filter; on synthesized sea-state transitions the crossfade interval
goes to 0.9865 / 0.9811 / 0.9730. Nothing measured here is large. What makes it
worth doing is that it is close to free: all eight records improve in 3D, seven
of the eight improve in vertical RMS, and the eighth is 0.03% worse.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_ou_ii/kalman_ou_ii-sim`, scoring the trailing 900 s of a 1200 s
replay. Every table is a ratio to the shipped `(0.6, 1.1)` on the same record
and seed, pooled by geometric mean, so records spanning `H_s` 0.27 m to 8.5 m
can be averaged. `tools/ou_ema_adapt_study.py confirm` reproduces all of them —
its `confirm` stage takes arbitrary environment knobs, and these two are
`OU_R_P0_COEFF` and `OU_R_V0_COEFF`.

## 1. Both coefficients are sharply, and separately, identified

One knob at a time over more than a decade of each, eight records at three IMU
seed triplets (n = 24):

| `R_p0_coeff` | 3D RMS | Z RMS | roll/pitch |
| --- | --- | --- | --- |
| 0.15 | 1.4820 | 1.4842 | 1.2646 |
| 0.25 | 1.2399 | 1.2188 | 1.0799 |
| 0.40 | 1.0725 | 1.0568 | 1.0154 |
| 0.60 (shipped) | 1.0000 | 1.0000 | 1.0000 |
| 0.90 | 1.0181 | 1.0489 | 0.9997 |
| 1.30 | 1.1138 | 1.2030 | 1.0033 |
| 2.00 | 1.2938 | 1.5200 | 1.0075 |
| 3.00 | 1.4907 | 1.9277 | 1.0101 |

| `R_v0_coeff` | 3D RMS | Z RMS | roll/pitch |
| --- | --- | --- | --- |
| 0.30 | 1.4176 | 1.8326 | 1.5712 |
| 0.50 | 1.1762 | 1.2099 | 1.1911 |
| 0.75 | 1.0464 | 1.0359 | 1.0521 |
| 1.10 (shipped) | 1.0000 | 1.0000 | 1.0000 |
| 1.60 | 1.0044 | 0.9977 | 0.9865 |
| 2.40 | 1.0230 | 1.0006 | 0.9857 |
| 3.50 | 1.0362 | 1.0027 | 0.9873 |

Both curves are steep on the tight side and shallower on the loose side, which
is the expected shape: too small an `r` is a hard pull of the state toward zero
that attenuates the wave itself, while too large an `r` merely stops
suppressing drift. Neither coefficient can be moved by a factor of two in
either direction without paying tens of percent, so the shipped pair was in the
right neighbourhood — and there is no version of this study in which these
coefficients are unimportant.

The two channels do different jobs, and the scan shows it. `r_p0` owns the
displacement error and barely moves attitude; `r_v0` owns roll/pitch — the
velocity channel is what the tilt and bias estimates compete with — and its
displacement curve is nearly flat from 1.1 upward. That is the same split
`docs/ou-ema-adaptation-tuning.md` found when it swept the two smoothing
horizons.

## 2. Where the optimum actually sits

Refining one knob at a time near the shipped pair, same protocol:

| `R_p0_coeff` | 3D RMS | Z RMS |
| --- | --- | --- |
| 0.45 | 1.0432 | 1.0312 |
| 0.50 | 1.0225 | 1.0142 |
| 0.55 | 1.0085 | 1.0043 |
| 0.60 (shipped) | 1.0000 | 1.0000 |
| 0.65 | 0.9959 | 1.0004 |
| 0.70 | 0.9955 | 1.0046 |
| 0.80 | 1.0028 | 1.0222 |

| `R_v0_coeff` | 3D RMS | Z RMS | roll/pitch |
| --- | --- | --- | --- |
| 0.85 | 1.0232 | 1.0167 | 1.0286 |
| 0.95 | 1.0095 | 1.0067 | 1.0135 |
| 1.05 | 1.0021 | 1.0015 | 1.0036 |
| 1.10 (shipped) | 1.0000 | 1.0000 | 1.0000 |
| 1.20 | 0.9981 | 0.9982 | 0.9947 |
| 1.30 | 0.9983 | 0.9974 | 0.9912 |
| 1.45 | 1.0008 | 0.9973 | 0.9881 |

`R_v0_coeff = 1.1` is simply not the optimum of anything: 1.2 to 1.45 is better
on 3D, on vertical, on roll/pitch and on yaw at once. `R_p0_coeff = 0.6` is the
vertical optimum, and the 3D optimum is somewhere above it.

They are not independent — both scale the regularization of one linear block —
so the pair was then run as a grid, 3D RMS above and Z RMS below:

| 3D RMS | `v0=1.1` | `1.2` | `1.3` | `1.45` | `1.6` |
| --- | --- | --- | --- | --- | --- |
| `p0=0.55` | 1.0085 | 1.0081 | 1.0094 | 1.0128 | 1.0169 |
| `0.60` | 1.0000 | 0.9981 | 0.9983 | 1.0008 | 1.0044 |
| `0.65` | 0.9959 | 0.9926 | **0.9917** | 0.9932 | 0.9962 |
| `0.70` | 0.9955 | 0.9907 | 0.9888 | 0.9893 | 0.9918 |
| `0.75` | 0.9980 | 0.9919 | 0.9890 | 0.9885 | 0.9904 |

| Z RMS | `v0=1.1` | `1.2` | `1.3` | `1.45` | `1.6` |
| --- | --- | --- | --- | --- | --- |
| `p0=0.55` | 1.0043 | 1.0029 | 1.0024 | 1.0024 | 1.0027 |
| `0.60` | 1.0000 | 0.9982 | 0.9974 | 0.9973 | 0.9977 |
| `0.65` | 1.0004 | 0.9981 | **0.9971** | 0.9968 | 0.9972 |
| `0.70` | 1.0046 | 1.0019 | 1.0006 | 1.0002 | 1.0005 |
| `0.75` | 1.0121 | 1.0088 | 1.0073 | 1.0067 | 1.0071 |

The two surfaces disagree about `r_p0` and agree about `r_v0`. Raising `r_p0`
keeps buying 3D error out to 0.75 and beyond, and starts costing vertical error
at 0.70. That is not a contradiction: the 3D gain from `r_p0` is mostly
horizontal — at `(0.65, 1.3)` the horizontal RMS improves on all eight records,
by 0.34% to 1.42% — and the horizontal and vertical channels of this filter are
not scored the same way by the sea. There is no anisotropy fix hiding here
either; `R_p0_xy_factor` is 1.0 by measurement, and at 0.65 the vertical
improves on seven of the eight records as well.

## 3. Where it stops: 0.65

Confirming the finalists at six IMU seed triplets (n = 48 per candidate):

| candidate | 3D RMS | Z RMS | Z worst | roll/pitch | yaw | accel bias |
| --- | --- | --- | --- | --- | --- | --- |
| `(0.60, 1.20)` | 0.9978 | 0.9982 | 1.0004 | 0.9959 | 0.9986 | 1.0000 |
| `(0.60, 1.30)` | 0.9977 | 0.9975 | 1.0008 | 0.9932 | 0.9976 | 1.0001 |
| `(0.60, 1.40)` | 0.9989 | 0.9973 | 1.0018 | 0.9916 | 0.9970 | 1.0001 |
| `(0.65, 1.25)` | 0.9908 | 0.9962 | 1.0112 | 0.9941 | 0.9975 | 1.0005 |
| **`(0.65, 1.30)`** | **0.9904** | **0.9958** | 1.0116 | **0.9929** | **0.9970** | 1.0006 |
| `(0.65, 1.35)` | 0.9905 | 0.9956 | 1.0120 | 0.9920 | 0.9966 | 1.0006 |
| `(0.70, 1.30)` | 0.9869 | 0.9981 | 1.0263 | 0.9930 | 0.9966 | 1.0010 |

The `roll/pitch` column is the hypotenuse of the two, and the two do not move
together: at the adopted point pitch reads 0.9796 and roll 1.0003, so the
attitude gain is entirely pitch and roll is flat.

`(0.70, 1.30)` is better than the adopted point on 3D RMS by another 0.35%, and
it is not adopted. The reason is the per-record mean vertical error, which is
what the deployed sentinels and the article's primary endpoint are built on:

| record | Z RMS at `(0.65, 1.3)` | Z RMS at `(0.70, 1.3)` |
| --- | --- | --- |
| JONSWAP H0.27 | 0.9957 | 0.9999 |
| JONSWAP H1.50 | 0.9981 | 1.0018 |
| JONSWAP H4.00 | 0.9911 | 0.9917 |
| JONSWAP H8.50 | 0.9931 | 0.9924 |
| PM-Stokes H0.27 | 0.9997 | 1.0057 |
| PM-Stokes H1.50 | 1.0003 | 1.0048 |
| PM-Stokes H4.00 | 0.9951 | 0.9958 |
| PM-Stokes H8.50 | 0.9935 | 0.9928 |

At 0.65 seven records improve and the eighth is 0.03% worse; at 0.70 three are
worse, two of them by half a percent, and all three are small seas. All four
JONSWAP records — the declared primary endpoint of the validation study is the
mean normalized vertical error over exactly that ensemble — improve at 0.65 and
two of them lose at 0.70.

So 0.65 is the largest position coefficient that leaves the per-record vertical
error where it is, and that, rather than an optimum of the pooled score, is the
rule that sets it. Going further trades a channel this filter is measured on
for one it is not.

The `Z worst` column above is the worst single (record, seed) ratio rather than
a mean, and 1.0116 at the adopted point is one realization of PM-Stokes H1.5.
It is reported because a pooled win that hides a local loss is not a win, but a
1% excursion on one seed of a record whose mean is flat is realization noise,
not a finding: the same column reads 1.0004 to 1.0018 for the `(0.60, v)` row,
which differs from the adopted point by 0.05 in one coefficient.

## 4. Why the optimum moved: it is the sigma channel

The shipped pair was not badly chosen; it was chosen for a filter that no
longer exists. Rebuilding OU-II at the commit before the parity change
(`3795296`, the last commit before OU-II took OU-III's `sigma_aw` band, startup
policy and hard-iron correction) and running the same scan there, against that
build's own `(0.6, 1.1)` baseline:

| candidate | 3D RMS | Z RMS | roll/pitch |
| --- | --- | --- | --- |
| `(0.50, 1.10)` | 1.0068 | 0.9891 | 1.0048 |
| `(0.55, 1.10)` | 1.0007 | 0.9914 | 0.9979 |
| `(0.60, 1.10)` shipped | 1.0000 | 1.0000 | 1.0000 |
| `(0.65, 1.10)` | 1.0033 | 1.0121 | 1.0002 |
| `(0.70, 1.10)` | 1.0100 | 1.0276 | 1.0005 |
| `(0.60, 0.90)` | 1.0087 | 1.0032 | 1.0134 |
| `(0.60, 1.30)` | 1.0019 | 1.0007 | 0.9947 |
| `(0.60, 1.60)` | 1.0101 | 1.0032 | 0.9930 |
| `(0.65, 1.30)` | 1.0038 | 1.0128 | 0.9953 |

On the pre-parity filter the shipped pair *is* the optimum — the `r_p0` optimum
sits at 0.55-0.60, `r_v0` at or just above 1.1, and the pair adopted here costs
0.38% of 3D RMS and 1.28% of vertical RMS. On the filter that ships today the
same pair gains 0.96% and 0.42%. The optimum moved; the constants did not.

The mechanism is visible directly in the tuner output. Comparing the two builds
on the four JONSWAP records at the end of the replay:

| record | `tau` before → after | `sigma_aw` before → after | `sigma_aw` ratio |
| --- | --- | --- | --- |
| H0.27 | 1.2796 → 1.2803 | 0.3852 → 0.3391 | 0.880 |
| H1.50 | 2.1796 → 2.1807 | 0.8963 → 0.6889 | 0.769 |
| H4.00 | 3.5933 → 3.5906 | 1.4027 → 1.0612 | 0.757 |
| H8.50 | 4.2244 → 4.2218 | 1.8602 → 1.3418 | 0.721 |

`tau` is unchanged to four figures — the wave-band period estimator was already
in place before the parity change — so this is one variable moving. The
band-passed variance channel reads 12% to 28% lower than the broadband one,
because it stops counting energy outside `[0.5, 4] f_tune`, and
`r = c * sigma_aw * tau^k` hands that reduction straight to the regularizer.

The re-fit gives part of it back and not all of it. Multiplying the two
columns, the regularizers the filter now applies, relative to what the
pre-parity filter applied on the same record:

| record | `r_v0` now / before | `r_p0` now / before |
| --- | --- | --- |
| H0.27 | 1.04 | 0.95 |
| H1.50 | 0.91 | 0.83 |
| H4.00 | 0.89 | 0.82 |
| H8.50 | 0.85 | 0.78 |

The band change cut both schedules by about 22% in a developed sea; the records
want a cut of 9% to 15% on the velocity channel and 17% to 22% on the position
channel, and the coefficients supply the difference. That is the shape to
expect: the band-passed channel is a better estimate of the wave-band
acceleration scale, so the optimal `r` should follow it partway rather than not
at all or one for one.

TFG hit the same coupling from the other direction when it inherited OU-III's
`R_S_coeff`: there the coefficient had been calibrated on the *other* filter's
broadband channel, and `R_S_coeff` came down from 0.35 to 0.28. Both are the
same statement. These coefficients are not dimensionless with respect to how
`sigma_aw` is measured, so they have to be re-fitted whenever that measurement
changes — which is a thing to check on the next change to a variance channel,
in either family, rather than a fact about this one.

That is the standing lesson, and it is why `tests/validation/test_record_conventions.py`
now pins OU-II's two coefficients and four clamps to the header the way it
already pinned OU-III's.

## 5. It also holds when the sea is moving

The coefficients set the target the EMA smoothers follow, so a pair fitted only
on stationary records could in principle be a pair that tracks a moving sea
badly. `tools/ou_ema_adapt_study.py transition` builds the non-stationary
instrument — two endpoint pairs, both directions, three wave-phase seeds, two
IMU seeds, with the 900 s window split into the pre-crossfade, crossfade and
post-crossfade intervals.

Crossfade interval, ratios to the shipped pair (n = 24):

| candidate | 3D RMS | Z RMS | roll/pitch | yaw | accel bias |
| --- | --- | --- | --- | --- | --- |
| `(0.60, 1.30)` | 1.0011 | 0.9958 | 0.9776 | 0.9962 | 0.9994 |
| `(0.65, 1.25)` | 0.9863 | 0.9820 | 0.9773 | 0.9959 | 0.9994 |
| **`(0.65, 1.30)`** | **0.9865** | **0.9811** | **0.9730** | 0.9950 | 0.9993 |
| `(0.70, 1.30)` | 0.9755 | 0.9701 | 0.9695 | 0.9941 | 0.9992 |

The transition instrument likes the adopted point more than the stationary one
does — 1.35% of 3D RMS against 0.96%, and 1.89% of vertical RMS against 0.42% —
so nothing about this pair is a stationary-record artifact. The full 900 s
window, which is three quarters stationary, reads 0.9876 / 0.9890. The `(0.70, 1.30)` ordering is unchanged here as well:
it wins the crossfade and loses the post-crossfade vertical (1.0106 against
1.0035), which is the same small-sea vertical cost that bounded it in section 3.

## 6. Adopted values, and what they cost the sentinels

| constant | was | now |
| --- | --- | --- |
| `R_p0_coeff` | 0.6 | **0.65** |
| `R_v0_coeff` | 1.1 | **1.3** |

`tools/ou_validation.py` mirrors both for its fixed-tuning modes and moves with
them. The regularizer floors are unaffected in the direction that matters —
raising a coefficient moves the schedule further from its lower clamp, and
`regularizer_floor-test` reports the smallest calibrated demand rising to
6.5x the `r_p0` floor and 55.6x the `r_v0` floor.

All nine deterministic quality gates still passed on their previous values, so
this is not a change that needed the bars moved. They were re-derived anyway,
with `tools/ou_regauge_gates.py --family ou_ii`, because the rule is that a
gate sits about half a percent above what the filter currently produces, and
after the re-fit five of the nine sat above that and four below — pitch at
0.0001 deg of margin, which is below this family's own measured `-march`
rebuild drift on that channel, so it had become a bar a rebuild decides.

Reverting the pair through `OU_R_P0_COEFF=0.6 OU_R_V0_COEFF=1.1` puts all nine
back at the rule to the digit, which is the control that says this set moved
for the re-fit and for nothing else.

| gate | was | now | worst scored, before → after |
| --- | --- | --- | --- |
| Z %Hs JONSWAP | 6.899 | **6.865** | 6.8644 → 6.8300 (H0.27) |
| Z %Hs PM-Stokes | 6.841 | 6.848 | 6.8062 → 6.8139 (H0.27) |
| yaw deg | 1.095 | **1.089** | 1.0895 → 1.0833 (JONSWAP H1.5) |
| roll deg | 0.4778 | 0.4792 | 0.4753 → 0.4768 (JONSWAP H4.0) |
| pitch deg | 0.3639 | 0.3657 | 0.3620 → 0.3638 (JONSWAP H8.5) |
| 3D % JONSWAP | 21.1 | **20.92** | 20.9867 → 20.8140 (H1.5) |
| 3D % PM-Stokes | 21.3 | **21.03** | 21.1935 → 20.9203 (H8.5) |
| acc Z bias % | 5.435 | **5.324** | 5.4073 → 5.2969 (JONSWAP H8.5) |
| bias 3D % | 94.37 | 94.47 | 93.8979 → 93.9911 (JONSWAP H4.0) |

Both 3D displacement gates come down by about 1.3%, which is where a change
aimed at the translational regularizer should show up.

Three of the four that loosen are single-realization moves against an ensemble
that goes the other way. Pooled over the eight records at six IMU seed triplets
the re-fit reads pitch 0.9796, roll 1.0003 and accelerometer-bias 3D 1.0006, so
pitch improves 2% while its deterministic worst record rises by 0.0018 deg, and
roll is flat inside its own realization noise. The fourth is not noise: the
small-sea vertical is the trade section 3 bounded on purpose, and 6.8062 →
6.8139 %Hs on PM-Stokes H0.27 is that bound being spent.

## 7. Reproducing

```sh
make -C tests/kalman_ou_ii build

# one knob at a time, eight stationary records, three IMU seed triplets
tools/ou_ema_adapt_study.py confirm --family OU_II \
    --records jonswap,pmstokes --seeds default,7,99 \
    --point p0=0.4:OU_R_P0_COEFF=0.4 --point v0=1.6:OU_R_V0_COEFF=1.6

# the finalists, six seed triplets
tools/ou_ema_adapt_study.py confirm --family OU_II \
    --records jonswap,pmstokes --seeds default,7,99,3,21,55 \
    --point p0.65/v1.3:OU_R_P0_COEFF=0.65:OU_R_V0_COEFF=1.3

# the same candidates against synthesized sea-state transitions
tools/ou_ema_adapt_study.py transition --family OU_II \
    --transition-stage confirm --transition-pair large,small \
    --transition-seeds 11,12,13 --seeds default,7 \
    --point p0.65/v1.3:OU_R_P0_COEFF=0.65:OU_R_V0_COEFF=1.3

# the deterministic gates, and the control that reverting the pair restores them
python3 tools/ou_regauge_gates.py --family ou_ii
OU_R_P0_COEFF=0.6 OU_R_V0_COEFF=1.1 python3 tools/ou_regauge_gates.py --family ou_ii
```

Section 4's attribution needs the pre-parity build:

```sh
git worktree add /tmp/ou2-preparity 3795296
ln -s "$PWD"/tests/kalman_ou_ii/wave_data_*.csv /tmp/ou2-preparity/tests/kalman_ou_ii/
cp tools/ou_ema_adapt_study.py /tmp/ou2-preparity/tools/
make -C /tmp/ou2-preparity/tests/kalman_ou_ii build
/tmp/ou2-preparity/tools/ou_ema_adapt_study.py confirm --family OU_II \
    --records jonswap,pmstokes --seeds default,7,99 \
    --point p0.65/v1.3:OU_R_P0_COEFF=0.65:OU_R_V0_COEFF=1.3
```

Both coefficients are also reachable at runtime through `setR_p0_Coeff()` and
`setR_v0_Coeff()`, and are part of `tools/ou_tuning_sweep.py`'s search space,
whose windows are re-centred on the adopted pair.
