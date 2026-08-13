# Bringing OU-II to parity with OU-III

OU-II and OU-III differ in exactly one thing that matters: the translational
state structure. OU-III carries an extra integral-displacement state `S` and
regularizes it with a single `r_S`; OU-II regularizes `p` and `v` separately
with `r_p0` and `r_v0`. Everything else in the two files — the auto-tuner, the
wave-period estimator, the private attitude observer, the magnetic acquisition,
the startup sequencing — is shared machinery that does not know which estimator
it is driving.

That machinery had drifted. OU-III accumulated changes that OU-II never
received, and each of them is about the shared part rather than about the
translational states, so OU-III's reasoning carries over unaltered. This
document records what was ported, what was deliberately not, and what it
measures.

Everything below is measured on the versioned simulation records
(`bareboat-necessities/oceanography-waves-lib`, `v1.1.3`) with
`tests/kalman_ou_ii/kalman_ou_ii-sim`, scoring the trailing 900 s of a 1200 s
replay.

## What was ported

### 1. The variance channel runs in a JONSWAP-similar wave band

`sigma_aw` is now estimated from the same exogenous levelled acceleration the
wave-period estimator uses, after `AdaptiveWaveBandPass` — a band-pass whose two
corners are fixed multiples (0.5 and 4) of the tuner's own wave frequency. Away
from the absolute safety clamps the transfer shape is therefore fixed in
`f/f_tune`, which is the condition the `sigma_aw` similarity argument needs. The
bench noise floor is referred through that band's own time-varying coefficients
before subtraction, rather than being subtracted as a broadband constant.

Two things changed together here, exactly as they did in OU-III: the band
itself, and the tuner's input signal moving from the raw body-Z proxy to the
complementary-levelled vertical acceleration. `setWaveBandTuning(false)` /
`W3D_TUNING_BAND=acceleration` bypasses both, so the old broadband path stays
available as a coherent ablation.

### 2. The pseudo-measurement cadence is self-similar in `tau`

One pseudo update has covariance `r^2`; at one update per `T_S` seconds the
continuous-equivalent information rate goes as `1/(r^2 T_S)`. Scaling `T_S` with
`tau` while holding `r` fixed would therefore change the regularization strength
with sea state as a side effect of the cadence, which is not what either
schedule is supposed to say.

`T_S = (0.015/1.1) * tau`, clamped to `[5, 250]` ms, with **both** filter inputs
renormalized by `sqrt(T_0/T_S)`. OU-II's `p` and `v` pseudo-updates fire on one
periodic tick inside the MEKF, so one cadence and one renormalization factor
serve both. The renormalization is deliberately not re-clamped, so the
small-sea end is allowed below the base floor once `T_S > T_0`.

The nominal operating point is preserved exactly: the wrapper's initial applied
`tau` is 1.1 s and `(0.015/1.1) * 1.1 = 0.015`, the historical period.

### 3. Startup: the Mahony proxy owns tilt and magnetic acquisition

`StartupInitPolicy::MahonyProxy` is the new default on the wrapper. The
measurement-only front end — proxy, frequency tracker, wave-period estimator,
sigma band, auto-tuner, wave-direction stage — runs from the first sample
through `updateFrontEnd()` with the MEKF untouched; the private Mahony observer
supplies the tilt that gates the magnetometer and frames the world-reference
average; and `goLive()` seeds the MEKF with the finished attitude so it starts
live in one step and never occupies the staged warmup.

`StagedMekf` restores the previous path and is the matched ablation
(`W3D_STARTUP_INIT=staged_mekf`). The inner filter still defaults to
`StagedMekf`, because only something above it can perform the handoff and a
filter driven directly through `updateTime()` would otherwise park at
`TunerReady` forever.

Three parts of this were not optional:

- **The proxy needs its integral term on.** `VerticalAccelComplementary`
  defaulted to `two_ki = 0`, which leaves about `2b/two_kp` of standing tilt.
  Everything else the observer feeds is high-passed and never noticed. Nothing
  high-passes an attitude seed. `startup_init-test` measures it directly: a
  0.05 deg/s constant gyro bias settles at **0.711 deg** of tilt with
  `two_ki = 0` and at zero with the deployed `two_ki = 0.02`.
- **The handoff needs a timeout.** Past `proxy_startup_timeout_sec` (150 s, and
  internally raised so it can never cut the magnetic acquisition short) the
  handoff proceeds on proxy tilt alone. A filter that silently produces nothing
  is a worse failure than one that starts from a stale operating point.
- **The seeded attitude covariance is anisotropic.** Tilt has been integrated
  through the wave band; yaw is either gauged by the magnetometer or entirely
  unknown, and those two cases are an order of magnitude apart. This needed a
  new `Kalman3D_Wave_OU_II::initialize_from_attitude()`.

### 4. Two-stage magnetic acquisition

A provisional reference locks as soon as the gravity gate allows, giving heading
and a live filter on the old schedule, and a second pass at
`mag_refine_start_sec` (90 s, 30 s window) re-learns the reference and re-gauges
heading once the observer has converged.

Two properties of that second pass are load bearing, and are the same two
OU-III found:

- **It is framed on the observer, not the MEKF.** By then the MEKF looks like
  the better frame, but it has been steering to the provisional reference the
  pass exists to replace, so its tilt carries that reference's error and
  averaging the field in it returns the same reference.
- **Accelerometer-bias learning waits for it** (`setAccBiasHold`). Bias and tilt
  error are barely separable in waves and the bias state has a 5000 s
  correlation time, so a value fitted to the provisional reference outlives the
  record.

### 5. The accelerometer-bias projection

The mean-reverting residual accelerometer-bias model was **already at parity**:
both filters run the same first-order Gauss–Markov residual about a
temperature-calibrated mean, `tau_b = 5000 s`, the same exact discrete
covariance (`-0.5 tau_b expm1(-2 T_s/tau_b)`), the same `Q_bacc` default and the
same `set_acc_bias_time_constant()` / `set_acc_bias_ou_stationary_std()` API.
Nothing there needed changing.

What OU-II was missing is the other half of the same argument:
`set_accel_bias_limit()` and the projection onto a ball of that radius after
every state injection. Mean reversion bounds the bias in distribution, not
pathwise, and the bias is excluded from the certified performance coordinate and
enters the ISS bound only as an input — so that input has to be bounded for the
bound to say anything. Default 0.5 m/s^2, loose enough never to bind on a
healthy MEMS unit, so it is a guarantee rather than a tuning parameter.

### 6. The continuous magnetic hard-iron correction

The last of the shared-machinery gaps, and the one with the largest single
effect on a scored channel. OU-III runs a body-fixed magnetometer offset
estimator for the life of the filter, in the private observer's yaw-stripped
tilt frame, and moves the magnetic reference by the delta the applied offset
implies rather than replacing it. OU-II now runs the same
`ContinuousMagHardIronEstimator` on identical defaults, wired at the same three
points of `updateMag()`.

Nothing about it is OU-III-shaped. The estimator reads gravity-referenced tilt
and the raw magnetometer, neither of which is a filter state, and it writes
measurement-model parameters rather than any state — so the translational
structure that distinguishes the two families does not enter, and neither
family's ISS argument is touched.

What it removes is a gauge rather than a tracking error: the startup
acquisition averages the field in a tilt frame and declares that direction to
be north, so a body-fixed offset lands in the world reference itself and every
consistent estimator inside the filter then agrees with the wrong north.
Eight-record mean yaw RMS falls from 1.887 to 0.813 deg and the worst record
from 2.161 to 1.089, which puts OU-II within three hundredths of a degree of
OU-III's 1.063 — as it should be, since the error belongs to the magnetometer
and not to either translational model. Pitch improves 0.289 to 0.255 deg; roll,
vertical and 3D displacement do not move.

The price is also OU-III's: the correction walks the heading onto the corrected
field during the run and the horizontal accelerometer bias absorbs part of that
motion, so the worst accelerometer-bias figures rise from 91.66% to 93.90% of
the true bias and the acc-Z gate from 5.33 to 5.41. That quantity's error
already exceeds the quantity itself under every configuration this family has
shipped, which is what a figure near 100% means.

`SF_MAG_CONT_HI=0` is the matched ablation and reproduces the pre-correction
filter to within 2.6e-4 relative — the same order as rebuilding at a different
`-march`. Full mechanism, regularization argument and the failure mode it
cannot detect: `docs/continuous-mag-hard-iron.md`.

## What was deliberately not ported

**The `r_S` floor change has no OU-II counterpart.** OU-III's 0.4 floor was not
a safety limit, it was the binding constraint on every low-motion sea: the
schedule asks for 0.24 m·s at the calibrated `H_s = 0.27` m point, so the floor
clipped it, and a full sweep of the tuner multipliers left both low-motion
records constant to three decimals. Dropping it to 0.15 recovered 8–10% there.

OU-II's laws differ only in the exponent, so the same failure was available and
had to be checked rather than assumed. It is not present. At the smallest
calibrated operating point (PM-Stokes `H_s = 0.27` m: `tau = 1.159 s`,
`sigma_aw = 0.390 m/s^2`) the schedule asks for `r_p0 = 0.314` m against a
0.05 m floor and `r_v0 = 0.497` m/s against a 0.01 m/s floor — clear by factors
of 6.3 and 50. Scaled down to the near-still `H_s = 0.05` m stress case the
demands are 0.058 m and 0.092 m/s, still above both floors. Lowering them
further would weaken a guard that is currently costing nothing.

`regularizer_floor-test` pins this against the operating points the eight scored
records actually produce, so a future coefficient re-fit that walks the schedule
down into a clamp fails rather than passes quietly. `OU_II_R_P0_MIN` /
`OU_II_R_V0_MIN` (and the `_MAX` counterparts) are now exposed on the simulator
so the question can be re-measured instead of re-argued.

**The `RSAdaptationLaw` ablation is not ported.** OU-III carries three
regularizer laws — `Cubic`, `StrongRiccati`, `PosteriorRiccati` — so that the
amplitude exponent of `r_S ~ sigma_aw^p tau^(5/2)` could be resolved by
measurement. The ablation confirmed the deployed `Cubic` law (`p = 1`) and the
non-Cubic laws exist for that ablation rather than as candidate defaults. There
is no OU-II analogue of the pole-placement derivation the Riccati laws come
from, and porting the machinery would add a configuration surface with no
deployed effect. Recorded as a known gap rather than a disagreement.

**The congruent `a_w` covariance-synchronization ablation is not ported.**
`Kalman3D_Wave_OU_II` has no
`synchronize_aw_covariance_to_stationary_congruent()`. OU-III's study concluded
the congruent alternative is worse, so this is an ablation surface, not a
behaviour difference.

## Results

Paired ensemble: five IMU noise seeds across the eight scored records, `n = 40`
paired records per metric, 900 s scoring window. `*` marks a paired mean
difference exceeding two standard errors.

Decomposed in the order the changes stack, so each column is the marginal effect
of one change with the others held:

| metric | baseline | + sigma band | + tau cadence | + proxy startup (deployed) |
| --- | --- | --- | --- | --- |
| Z RMS, %Hs | 6.318 | 6.343 | 6.343 | 6.350 |
| 3D RMS, % of max abs disp | 18.997 | 19.074 | 19.073 | 19.042 |
| roll RMS, deg | 0.354 | 0.348 | 0.348 | 0.286 |
| pitch RMS, deg | 0.269 | 0.308 | 0.308 | 0.300 |
| yaw RMS, deg | 2.488 | 2.359 | 2.360 | 2.411 |
| accel-bias 3D RMS, m/s^2 | 0.0595 | 0.0616 | 0.0616 | 0.0520 |
| accel-bias Z, % of true | 4.141 | 4.049 | 4.049 | 4.057 |
| accel-bias 3D, % of true | 82.94 | 85.73 | 85.70 | 70.14 |

Net effect of the whole change, baseline against the deployed default:

| metric | delta | rel | significant |
| --- | --- | --- | --- |
| roll RMS | -0.068 +/- 0.022 deg | **-19.3%** | `*` |
| accel-bias 3D, % of true | -12.80 +/- 5.99 pp | **-15.4%** | `*` |
| accel-bias 3D RMS | -0.0076 +/- 0.0038 m/s^2 | -12.7% | |
| yaw RMS | -0.076 +/- 0.018 deg | **-3.1%** | `*` |
| accel-bias Z, % of true | -0.083 +/- 0.021 pp | -2.0% | `*` |
| 3D RMS % | +0.045 +/- 0.049 pp | +0.24% | |
| Z RMS %Hs | +0.032 +/- 0.018 pp | +0.50% | |
| pitch RMS | +0.032 +/- 0.007 deg | **+11.8%** | `*` |

Read together with the decomposition:

**Displacement accuracy does not move.** Vertical RMS goes 6.318 -> 6.350 %Hs
and 3D goes 18.997 -> 19.042 %, both inside two standard errors. This is a
change to the attitude front end and to how `sigma_aw` is measured, not to the
translational estimator, and it prices as one.

**Attitude and accelerometer bias improve, and it is the startup policy that
does it.** Roll -17.9% and the accelerometer-bias aggregate -18.2% come almost
entirely from the `+ proxy startup` column; the sigma band contributes roughly
nothing to either. That is the expected shape: the reference and the yaw gauge
are locked once, so removing the warming MEKF's tilt error from the frame they
are locked in is worth a standing bias for the whole run.

**The tau-scaled cadence is measurably a no-op**, and by construction. Every
metric moves by under 0.05% between the second and third columns. The
renormalization holds `r^2 T_S` fixed, so the continuous-equivalent
regularization is identical and only the discretization granularity changes —
and at `tau` of a few seconds, 15 ms and 30 ms are both far finer than anything
in the band. It is ported for consistency of policy across the family, not for
a gain it does not produce. Its value is that the two filters now answer the
same way to "what sets the regularizer strength?".

**Pitch is the one real cost, and it comes from the sigma band.** Pitch RMS
rises 11.8% overall, 14.8% of it in the `+ sigma band` column. Roll falls
further than pitch rises in absolute terms, so combined tilt RMS
(`sqrt(roll^2 + pitch^2)`) improves from 0.444 to 0.414 deg, a 6.7% reduction,
and yaw improves 3.1% on top. The redistribution between the two horizontal
axes is not symmetric because the records are not: all eight carry a fixed
`+/-30 deg` wave direction. Worth flagging rather than averaging away.

**Yaw improves overall, but not monotonically.** The sigma band takes yaw down
5.2%; the proxy startup gives 2.2% of that back, which is the same ~2% yaw cost
OU-III measured and attributed to the learned reference vector rather than to
the one-time gauge.

### Startup timing

First heading and time-to-live, deployed configuration, from the deterministic
run:

| | first heading | live | mag refined |
| --- | --- | --- | --- |
| best record | 22.0 s | 22.0 s | 120.1 s |
| worst record | 52.0 s | 109.1 s | 139.1 s |

The two records that reach live late are the large JONSWAP seas, where the
gravity-agreement gate takes longer to hold. That is the gate working, not
failing: the handoff waits for a tilt it can trust, and the timeout is there so
it can never wait forever.

### Quality gates

The `kalman_ou_ii-sim` regression sentinels were re-derived on the filter that
now ships, following the existing convention of worst scored value plus about
half a percent, rounded up in the last digit the channel is quoted in. The
column below marked *parity* is the state after items 1-5; *now* is after the
continuous hard-iron correction of item 6.

| gate | before | parity | now | worst scored |
| --- | --- | --- | --- | --- |
| Z %Hs, jonswap | 7.0 | 6.9 | 6.9 | 6.8638 (H0.27) |
| Z %Hs, pmstokes | 6.9 | 6.9 | 6.9 | 6.8061 (H0.27) |
| yaw deg | 2.2 | 2.18 | **1.10** | 1.0895 (jonswap H1.5) |
| 3D %, jonswap | 20.9 | 21.1 | 21.1 | 20.99 (H1.5) |
| 3D %, pmstokes | 20.7 | 21.2 | **21.3** | 21.19 (H8.5) |
| acc Z bias % | 5.9 | 5.4 | **5.5** | 5.41 (jonswap H8.5) |
| bias 3D % | 85.4 | 92.2 | **94.4** | 93.90 (jonswap H4.0) |

See the comment at `FAIL_LIMITS` for why the gates that loosen are realization
and gauge moves rather than quality regressions, and for the ensemble figures
that establish it.

## Knobs added

| setting | default | what it does |
| --- | --- | --- |
| `startup_init_policy` | `MahonyProxy` | selects the startup policy |
| `proxy_startup_min_sec` | 8 | earliest handoff |
| `proxy_startup_timeout_sec` | 150 | latest handoff; raised internally so it cannot cut mag acquisition short |
| `proxy_handoff_tilt_sigma_rad` | 0.035 | seeded tilt variance |
| `proxy_handoff_yaw_sigma_rad` | 0.087 | seeded yaw variance once north is gauged |
| `proxy_mag_settle_sec` | 0 | observer settling before the provisional lock |
| `mag_refine_enabled` | true | run the second-stage acquisition |
| `mag_refine_start_sec` | 90 | when the refinement begins |
| `mag_refine_window_sec` | 30 | its averaging window |
| `acc_bias_unlock_mag_updates` | 250 | magnetometer updates after going live before the bias gate opens |
| `sigma_band_low_ratio` / `_high_ratio` | 0.5 / 4.0 | sigma-band corners in units of `f_tune` |
| `sigma_band_min_hz` / `_max_hz` | 0.01 / 6.0 | absolute safety clamps on those corners |
| `setTauScaledPseudoUpdateCadence` | true | self-similar cadence; false restores fixed 15 ms |
| `setPseudoUpdateTauRatio` | 0.015/1.1 | `T_S/tau` |
| `setPseudoUpdatePeriodBounds` | 0.005 / 0.25 s | clamps on `T_S` |
| `set_accel_bias_limit` | 0.5 m/s^2 | radius of the accelerometer-bias projection ball |
| `mag_continuous_hard_iron` | true | run the continuous hard-iron estimator |
| `mag_hi_memory_sec` | 600 | exponential memory of its statistics |
| `mag_hi_model_ridge` / `_relative` | 4e-3 / 0.5 | absolute ridge floor, and the part that scales with excitation |
| `mag_hi_min_information` | 2.0 | excitation the window must reach |
| `mag_hi_min_effective_weight` | 500 | effective samples before any answer |
| `mag_hi_max_residual_rms_uT` | 3.0 | model-fit gate; what a turning hull trips |
| `mag_hi_max_bias_fraction` | 0.35 | plausibility bound against the field norm |
| `mag_hi_apply_fraction` | 1.0 | shrinkage on top of the ridge |
| `mag_hi_slew_tau_sec` | 45 | how fast the applied offset moves |

Simulator env overrides: `W3D_STARTUP_INIT`, `W3D_PSEUDO_CADENCE`,
`SF_PROXY_START_MIN_SEC`, `SF_PROXY_START_TIMEOUT_SEC`, `SF_PROXY_MAG_SETTLE_SEC`,
`SF_MAG_REFINE`, `SF_MAG_REFINE_START_SEC`, `SF_MAG_REFINE_WINDOW_SEC`,
`SF_PROXY_TILT_SIGMA`, `SF_PROXY_YAW_SIGMA`, `SF_ACC_BIAS_UNLOCK_MAG_UPDATES`,
`SF_MAG_CONT_HI`, `SF_MAG_HI_MEMORY_SEC`, `SF_MAG_HI_RIDGE`, `SF_MAG_HI_RIDGE_REL`,
`SF_MAG_HI_MIN_INFO`, `SF_MAG_HI_MIN_WEIGHT`, `SF_MAG_HI_MAX_RESID`,
`SF_MAG_HI_FRACTION`, `SF_MAG_HI_SLEW_TAU`,
`OU_II_PSEUDO_TAU_RATIO`, `OU_II_PSEUDO_PERIOD_MIN_S`, `OU_II_PSEUDO_PERIOD_MAX_S`,
`OU_II_R_P0_MIN`/`_MAX`, `OU_II_R_V0_MIN`/`_MAX`.

## Reproducing

```
make -C tests/kalman_ou_ii build

# deployed default, eight scored records
W3D_WRITE_TIMESERIES=0 W3D_COLLECT_ALL_GATES=1 ./kalman_ou_ii-sim

# the three ablations that decompose the change
W3D_STARTUP_INIT=staged_mekf W3D_PSEUDO_CADENCE=fixed ./kalman_ou_ii-sim   # sigma band only
W3D_STARTUP_INIT=staged_mekf ./kalman_ou_ii-sim                            # + tau cadence
W3D_TUNING_BAND=acceleration ./kalman_ou_ii-sim                            # legacy broadband path
```

The ensemble figures above come from replaying each configuration under
`W3D_IMU_SEED=1..5` and pairing on `(seed, record)`.
