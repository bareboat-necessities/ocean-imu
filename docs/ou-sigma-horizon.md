# Taking the OU operating point off the acceleration-band frequency tracker, and re-gauging the `sigma_a` horizon

Two changes, one cause. The OU adaptation path had one remaining read of the
acceleration-band frequency tracker, and that read also fixed the units of the
`sigma_a` averaging horizon. Removing the read makes the horizon mean something
different, so it has to be re-measured.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_ou_ii/kalman_ou_ii-sim` and
`tests/kalman_ou_iii/kalman_ou_iii-sim`, scoring the trailing 900 s of a 1200 s
replay. `tools/ou_sigma_horizon_study.py` reproduces every table.

**Result.** The tracker is gone from the adaptation path and costs nothing:
pooled vertical RMS moves by 0.00%, worst record 0.11%. `K_periods` stays at
`2.0`. It is not the value anyone would pick from first principles now that it
means something different, but the sweep says the alternatives net to nothing:
re-measured against a crossfade fast enough to lag (section 5), the sign of the
effect reverses between a building and a decaying sea, so no constant captures
it, and `2.0` is where the ensemble's worst realization sits. What the sweep
does establish is that `2.0` is on the correct side of the optimum -- shortening
it costs on the stationary ensemble, 0.25-0.33% at `K = 1` and 0.4-1.0% at
`K = 0.5` (section 4) -- and the constant is now a settable, measured knob
rather than an inherited one.

## 1. What was still coupled

`SeaStateFusionFilter_OU_II` and `_OU_III` run two frequency estimators side by
side, and they measure different things:

| estimator | band | what it is for |
| --- | --- | --- |
| `TrackerPolicy` (KalmANF / Aranovskiy / PLL / Schmitt) | acceleration | carrier of the wave-direction demodulator |
| `WavePeriodEstimator` | wave (elevation) | `T_z`, from which the OU operating point is built |

They are not interchangeable. An ocean acceleration spectrum is the elevation
spectrum weighted by `(2*pi*f)^4`, so its apparent frequency sits well above the
spectral peak and barely moves as the sea grows: across the reference family the
elevation zero-crossing period spans 2.3-8.4 s while the acceleration band stays
near 0.42-0.55 Hz.

The operating point already came from `WavePeriodEstimator` — but only once that
estimator reported `isReady()`, and that gate does not clear early:

| record | wave period ready | Live |
| --- | --- | --- |
| JONSWAP H_s = 0.27 m | 62.2 s | 22.0 s |
| JONSWAP H_s = 1.50 m | 70.8 s | 52.1 s |
| JONSWAP H_s = 4.00 m | 82.7 s | 109.1 s |
| JONSWAP H_s = 8.50 m | 84.4 s | 107.1 s |
| PM-Stokes H_s = 0.27 m | 61.9 s | 22.0 s |
| PM-Stokes H_s = 8.50 m | 77.9 s | 150.0 s |

Until then `tuner_frequency_hz_()` returned the tracker's output, and three
things downstream consumed it: the corners of `AdaptiveWaveBandPass`, the
`sigma_a` averaging horizon inside `SeaStateAutoTuner`, and `tau` itself. On the
two smallest seas the filter is already Live 40 s before the period estimator
is, so it spent that time adapting on an acceleration-band frequency.

Two smaller couplings went with it. The tuning frequency was clamped to
`[MIN_TUNE_FREQ_HZ, MAX_FREQ_HZ]`, the upper bound being the *tracker's*, so
`setFreqBounds()` — a wave-direction knob — reached the OU operating point.
It is now clamped to `[MIN_TUNE_FREQ_HZ, MAX_TUNE_FREQ_HZ]`, neither of which
binds on any reference record.

The stillness detector still shares the tracker's low-pass input, and that is
not a tracker dependence: it reads the filtered acceleration, never the
tracker's output frequency.

## 2. The replacement

`TunerFrequencySource::WaveBand`, now the default in both families, changes two
things at once, so a third value exists to separate them:

| source | before the estimator has a value | when it has one |
| --- | --- | --- |
| `TrackerFallback` (legacy) | tracker output | at `isReady()` |
| `WaveBandGated` (ablation) | 0.2 Hz prior | at `isReady()` |
| `WaveBand` (deployed) | 0.2 Hz prior | as soon as it is finite |

The prior is a constant because that is what makes the schedule exogenous at
every instant of the run rather than only after the gate clears — a constant is
trivially a pure function of no measurement at all. 0.2 Hz is a 5 s
zero-crossing period, against the 2.3-8.4 s the estimator reports across the
reference family, and it is the same constant `SeaStateFusionFilter_TFG` has always used
in this position. This change brings the two OU families into line with TFG
rather than inventing a policy for them.

`WaveBand` also stops waiting for `isReady()`. That gate wants a settled
*statistic*; a value that has merely survived the integrators' settling
transient is already a far better wave-band estimate than a constant, and it
appears about 50 s into a run instead of 60-85 s.

### What it costs

Ratios to `tracker_fallback` on the identical realization, eight stationary
records plus one crossfade, default seeds:

| family | metric | `wave_band_gated` | `wave_band` | largest single-record change |
| --- | --- | --- | --- | --- |
| OU-II | vertical RMS | 1.0000 | 1.0000 | +0.02% |
| OU-II | 3D RMS | 1.0000 | 1.0000 | +0.03% |
| OU-III | vertical RMS | 1.0000 | 1.0000 | +0.01% |
| OU-III | 3D RMS | 1.0000 | 0.9999 | -0.11% |

Nothing, to three decimal places, and nothing on the crossfade intervals either
(largest segment ratio 0.9993, on OU-II's blend). That is the expected answer and the reason the
change is safe: the scoring window opens at 300 s, by which time both sources
have long since converged on the same wave-band frequency, and the filter's
memory of the first minute has washed out.

The value of the change is therefore not in the scored window. It is that the
adaptation path is now a function of one physically correct band at every
instant, that `setFreqBounds()` and the choice of tracker no longer reach it,
and that the exogeneity argument in the stability appendix holds from the first
sample rather than from `isReady()`.

`tests/kalman_ou_iii/tuner_coupling-test.cpp` pins this: two filters differing
in nothing but their `TrackerType` are driven with identical samples, and at
30 s — where the two trackers report 0.297 and 0.112 Hz — `tau`, `sigma_a`,
`r_S` and the sigma band's corners must be bit-identical. Under
`TrackerFallback` that check fails, on the settled comparison as well as the
early one: the divergence the tracker injects in the first minute survives in
the tuner's own EMAs long after the period estimator has taken over.

## 3. Why the `sigma_a` horizon had to be re-gauged

`SeaStateAutoTuner` averages the band-limited acceleration variance over
`tau_var = K_periods * T_eff`, where `T_eff` is the period of whatever tuning
frequency it is fed. It is a two-stage EWMA — mean and square, then the variance
— so the memory is about `2 * tau_var`.

`K_periods` has been 2.0 since the tuner was written, when `T_eff` was the
acceleration band's period. That band does not move with sea state, so `K = 2`
meant a horizon of about 4 s on every record. Once the operating point moved to
the wave band the same constant became `2 * T_z`:

| record | `T_z` | horizon at `K = 2` | memory | horizon at `K = 4` | memory |
| --- | --- | --- | --- | --- | --- |
| JONSWAP H_s = 0.27 m | 2.60 s | 5.2 s | 10.4 s | 10.4 s | 20.8 s |
| JONSWAP H_s = 1.50 m | 4.37 s | 8.7 s | 17.5 s | 17.5 s | 35.0 s |
| JONSWAP H_s = 4.00 m | 7.09 s | 14.2 s | 28.3 s | 28.3 s | 56.7 s |
| JONSWAP H_s = 8.50 m | 8.42 s | 16.8 s | 33.7 s | 33.7 s | 67.4 s |

Nobody chose those numbers; they are what the old constant became when the
frequency underneath it changed by a factor of two to four. That is what this
section measures.

## 4. The stationary records: shorter is clearly wrong, longer is nearly free

`tools/ou_sigma_horizon_study.py --axis k_periods --seeds 6`. Eight records at
six IMU seeds, paired on the realization, geometric mean of the per-run ratio to
`K = 2` with a normal-approximation 95% interval (n = 48):

| `K_periods` | OU-III vertical | OU-III 3D | OU-II vertical | OU-II 3D |
| --- | --- | --- | --- | --- |
| 0.5 | 1.0038 [1.0009, 1.0066] | 1.0056 [1.0041, 1.0071] | 1.0096 [1.0065, 1.0127] | 1.0078 [1.0059, 1.0096] |
| 1 | 1.0025 [1.0009, 1.0040] | 1.0023 [1.0016, 1.0030] | 1.0033 [1.0017, 1.0049] | 1.0031 [1.0021, 1.0042] |
| 2 (shipped) | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 3 | 0.9991 [0.9983, 0.9999] | 0.9996 [0.9992, 1.0000] | 0.9993 [0.9985, 1.0001] | 0.9991 [0.9986, 0.9996] |
| 4 | 0.9989 [0.9977, 1.0001] | 0.9996 [0.9990, 1.0003] | 0.9992 [0.9980, 1.0004] | 0.9988 [0.9981, 0.9995] |
| 6 | 0.9987 [0.9971, 1.0003] | 0.9995 [0.9987, 1.0004] | 0.9992 [0.9976, 1.0007] | 0.9985 [0.9975, 0.9994] |
| 8 | 0.9986 [0.9968, 1.0003] | 0.9994 [0.9984, 1.0004] | 0.9990 [0.9974, 1.0007] | 0.9982 [0.9972, 0.9992] |

The curve is asymmetric. Halving `K` costs 0.4-1.0%, and the interval excludes
1 comfortably in every column: at `K = 0.5` the horizon is 1.3-4.2 s and the
`sigma_a` estimate is passing wave-to-wave variance through into `r_S` rather
than estimating the sea state. Doubling it, by contrast, buys about 0.1% and
then saturates: `K = 3` and `K = 8` are indistinguishable. Whatever the
stationary records are measuring, it is already finished by three wave periods.

## 5. Sea-state transitions: the two directions disagree

A stationary record cannot price a lag. The instrument for that is the
crossfade: `tools/ou_ema_adapt_study.py transition` builds records that blend
two independently phase-randomized realizations, in both directions and on two
endpoint pairs an octave apart in wave period, and scores the intervals of the
900 s window separately (12 records: 2 pairs x 2 directions x 3 wave-phase
seeds, OU-III, here at 3 IMU seed triplets each, n = 36).

That instrument has since been sharpened twice. The crossfade now runs over
540-660 s instead of 420-780 s -- 120 s rather than 360 s, which is about three
times the longest averaging memory in the filter instead of ten times it -- and
the endpoint sea is scored as two intervals rather than one: a `run-on`
interval one crossfade long, where a schedule that averages too long is still
carrying the sea it just left, and the `settled` remainder. The numbers below
are that instrument, on the deployed filter (`ADAPT_RS_MULT = 1.5`).

Vertical RMS, geometric mean of the ratio to `K = 2`, and the worst single
realization:

| `K_periods` | blend | blend worst | run-on | settled | whole window |
| --- | --- | --- | --- | --- | --- |
| 1 | 1.0057 | 1.0378 | 0.9990 | 0.9975 | 0.9980 |
| 1.5 | 1.0018 | 1.0138 | 0.9999 | 0.9988 | 0.9990 |
| 2 (shipped) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 3 | 0.9968 | 1.0201 | 0.9989 | 1.0017 | 1.0010 |
| 4 | 0.9928 | 1.0318 | 0.9972 | 1.0027 | 1.0012 |
| 6 | 0.9848 | 1.0426 | 0.9943 | 1.0035 | 1.0008 |

The whole window still spans 0.998 to 1.001 across a six-fold change in `K`, so
the conclusion of the first pass survives: there is no fixed `K` worth moving
to. What the sharper instrument adds is *why*, and the reason is not the one
the first pass proposed. Splitting by transition direction:

| `K_periods` | blend, sea building | blend, sea decaying | run-on, building | run-on, decaying |
| --- | --- | --- | --- | --- |
| 1 | 1.0267 | 0.9947 | 0.9959 | 0.9988 |
| 1.5 | 1.0104 | 0.9969 | 0.9982 | 1.0010 |
| 3 | 0.9823 | 1.0056 | 1.0022 | 0.9941 |
| 4 | 0.9659 | 1.0091 | 1.0030 | 0.9874 |
| 6 | 0.9402 | 1.0125 | 1.0022 | 0.9776 |

(`large` endpoint pair, H_s = 1.5 m to 4.0 m and back; the `small` pair has the
same signs and a third of the amplitude.)

The effect does not cancel because two intervals disagree. It cancels because
the *same* interval wants opposite things depending on which way the sea is
going, and the ensemble weights both directions equally:

- A long `sigma_a` memory lags, and a lagging estimate is biased **low** while
  the sea builds and **high** while it decays.
- The filter is not symmetric in that bias. Under-estimating `sigma_a` shrinks
  `r_S ~ sigma_a tau^3`, which tightens the integral anchor, and during a fast
  change that is worth up to 6% of vertical RMS (`K = 6`, sea building).
  Over-estimating it loosens the anchor and costs about 1.2%.
- The run-on interval mirrors it with the sign flipped, because that is where
  the memory of the sea just left is still in the average: the decaying case
  gains there (0.9776 at `K = 6`) exactly where the building case pays.

So the pooled blend gain of a long horizon is not better tracking. It is a
favourable bias in one direction of travel, and it is paid back in the other
direction and in the settled sea. A constant cannot buy it: the sign of the
correction it would need is the sign of `d(sigma_a)/dt`.

The worst single realization is what closes the case. On the blend it is
minimized at `K = 1.5` and rises in both directions -- 1.0138 at 1.5, 1.0201 at
3, 1.0378 at 1, 1.0426 at 6 -- and the same holds on the settled interval
(1.0017 at 1.5 against 1.0087 at 3), on the start interval and over the whole
window (1.0022 at 1.5, 1.0057 at 3). The run-on interval is the one exception,
narrowest at `K = 3` with 1.0061 against 1.0108 at 1.5. So on every interval a
fixed `K` has to be safe on, the shipped value is at or one grid step from the
narrowest tail, which is more than "the means cancel here". `K` stays at 2.0.

Getting the blend gain *and* the settled recovery needs a horizon that
shortens when the sea state actually moves -- the same two-regime problem
`adaptiveSmoothingHorizonSec()` already solves for the `r_S` channel with its
`slew_log` term. That is not attempted here: the detector would have to compare
two variance timescales rather than a target against an applied value, and it
introduces a constant of its own that this evidence cannot fit. The `r_S`
channel is where the sharpened instrument did move a constant; see
docs/ou-ema-adaptation-tuning.md, section 7.

## 6. Two constants that turned out not to matter

**The wave-band prior.** `tools/ou_sigma_horizon_study.py --axis prior --seeds 3`
sweeps it over 0.1-0.4 Hz, a factor of four, i.e. a 2.5-10 s assumed
zero-crossing period. Every pooled ratio is `1.0000` with a 95% interval of
`[1.0000, 1.0000]`, in both families, on the stationary records and on every
transition interval. The prior is replaced by a measurement about 50 s into
the run and the window opens at 300 s, so this is the expected answer; it is
recorded because a constant that reaches the shipped filter should be shown not
to matter rather than assumed not to.

**The absolute horizon ceiling.** `--axis horizon_max` at `K = 8`, where the
ceiling is closest to binding:

| ceiling | OU-III vertical | OU-III blend |
| --- | --- | --- |
| 20 s | 1.0005 | 1.0185 |
| 30 s | 1.0002 | 1.0109 |
| 60 s (shipped) | 1.0000 | 1.0000 |
| 120 s | 1.0000 | 0.9999 |
| 240 s | 1.0000 | 0.9999 |

Raising it changes nothing, which is what confirms section 5 is measuring `K`
and not the clamp: at `K = 8` the requested horizon reaches 67 s on the largest
sea, so the ceiling trims it by 10% there and by nothing anywhere else. Lowering
it to 20 or 30 s does bind, and reproduces the short-horizon penalty from the
other direction.

## 7. What shipped

| | before | after |
| --- | --- | --- |
| tuning frequency source | tracker until `wavePeriodReady()` | `WavePeriodEstimator`, else a 0.2 Hz constant |
| tuning frequency ceiling | `MAX_FREQ_HZ` (tracker's) | `MAX_TUNE_FREQ_HZ` |
| `K_periods` | 2.0, inherited | 2.0, measured |

New API on both filters, and the environment overrides the simulators read:

| setter | environment variable |
| --- | --- |
| `setTunerFrequencySource()` | `W3D_TUNER_FREQ_SOURCE` |
| `setTuneFreqPriorHz()` | `OU_TUNE_FREQ_PRIOR_HZ` |
| `setTuneFreqBounds()` | — |
| `setSigmaVarianceKPeriods()` | `OU_SIGMA_VAR_K_PERIODS` |
| `setSigmaVarianceHorizonBounds()` | `OU_SIGMA_VAR_HORIZON_MIN_S`, `OU_SIGMA_VAR_HORIZON_MAX_S` |

`getSigmaVarianceHorizonSec()` reports the horizon currently in force, so the
number in the table of section 3 can be read out of a running filter instead of
recomputed.

## 8. One quality-gate sentinel moved

The OU-III simulator's pitch sentinel was cut at 0.2211 deg, half a percent
above the 0.2200 deg the filter produced on PM-Stokes H_s = 4.0 m at the
default seed. After the change that record reports 0.2218 deg, 0.83% up, so
the gate is re-cut to 0.223 by `tools/ou_regauge_gates.py` and the rule in
`docs/quality-gate-regauge.md`.

That is a re-draw rather than a pitch regression, and it is worth saying why
the distinction is not special pleading. Paired over five IMU seeds and all
eight records, the pitch ratio is 0.9999 with a 95% interval of
[0.9989, 1.0010] — no systematic effect, resolved ten times finer than the move
on this one record. Roll is 0.9999 [0.9997, 1.0001] and yaw 0.9998
[0.9992, 1.0004]. The one record that moves is the one the gate happens to be
written against, and it moves inside its own scatter.

None of the other eight limits moved, and none were re-cut.
