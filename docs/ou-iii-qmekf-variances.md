# The MEKF variances OU-III never swept

**Outcome: the three sensor variances move, and every gated channel improves.**
Roll $-13.3\%$, pitch $-14.7\%$, yaw $-1.2\%$, vertical displacement $-2.4\%$,
3D displacement $-37.2\%$, accelerometer bias $-1.0\%$, gyroscope bias
$-36.5\%$, paired over eight records and five seeds. All nine deterministic
sentinels come down; none goes up. Two of the seven variances turn out to be
incapable of affecting anything at all.

`Kalman3D_Wave_OU_III` takes seven variances through its constructor:

```cpp
Kalman3D_Wave_OU_III(Vector3 const& sigma_a,
                     Vector3 const& gyro_noise_density_rad_sqrt_s,
                     Vector3 const& sigma_m,
                     T Pq0 = T(5e-4), T Pb0 = T(1e-6), T b0 = T(1e-11),
                     T R_S_noise_var = T(1.5),
                     T gravity_magnitude = T(STD_GRAVITY));
```

Every other tuning constant in this family has a study behind it — the
$r_S$ law, the anisotropy, the EMA horizons, the sigma horizon, the startup
policy. These seven had none. They were picked once and never revisited.

## Four of them could not be swept even in principle

`SeaStateFusion_OU_III::begin()` called the three-argument `initialize()`:

```cpp
impl_.initialize(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m);   // before
```

`initialize_ext()`, which takes the other four, existed and was called by
nothing. Every deployment of OU-III has therefore run on the header defaults
for `Pq0`, `Pb0`, `b0` and `R_S_noise_var`, and no caller could have changed
them without editing the header. They are `Config` fields now, carrying those
same defaults, and `begin()` routes through `initialize_ext()`.

The simulator exposes all seven. The three sensor sigmas arrive from the
shared harness already built as a fixed multiple of the white noise it
injects, so they are swept as scale factors on that multiple:

| variable | meaning | deployed |
| --- | --- | --- |
| `SF_SIGMA_A_SCALE` | scale on `sigma_a`, harness $2.8\times$ injected accel white | 1 |
| `SF_SIGMA_G_SCALE` | scale on `sigma_g`, harness $2.0\times$ injected gyro white | 1 |
| `SF_SIGMA_M_SCALE` | scale on `sigma_m`, harness $1.2\times$ injected mag white | 1 |
| `SF_PQ0` | initial attitude-error variance, rad² | 5e-4 |
| `SF_PB0` | initial gyro-bias variance, (rad/s)² | 1e-6 |
| `SF_GYRO_BIAS_RW_VAR` | gyro-bias random-walk variance density, (rad/s)²/s | 1e-11 |
| `SF_RS_NOISE_VAR` | initial integral pseudo-measurement variance, (m·s)² | 1.5 |

`tools/ou_qmekf_variance_study.py` runs the sweep. Every arm is paired against
the deployed point on the same `(record, seed)`, over the four JONSWAP and four
PM-Stokes records the gates score, on the trailing 900 s of each replay.

## Two of the seven are inert

`Pq0` moved over three decades and `R_S_noise_var` over four leave every
metric **bit-for-bit identical**. Not approximately: `+0.00 [+0.0, +0.0]` on
all six reported channels, all 40 pairs.

This is not a dead wire. `Pb0` and `b0` ride the same `initialize_ext()` call
and do move results. Both inert parameters are overwritten before they can be
read:

- `Pq0` initializes the attitude block of `Pext`. Every startup path —
  `initialize_from_acc()`, `initialize_from_acc_mag()`,
  `initialize_from_attitude()` — then calls
  `set_accel_only_attitude_covariance_()`, which rewrites exactly that block.
  It is not only the Mahony handoff doing this: `Pq0` stays inert under
  `W3D_STARTUP_INIT=staged_mekf` too.
- `R_S_noise_var` initializes `R_S`. `apply_RS_tune_()` rewrites it at
  `enterLive_()`, before the linear block issues its first $S=0$ update, and
  the tuner rewrites it on every adaptation tick thereafter.

Neither is a bug to fix — the values that overwrite them are the better ones,
derived rather than guessed. But they are documentation problems: both read
as tuning knobs, and one of them, `R_S_noise_var`, sits in the constructor
signature next to six parameters that do matter. Anyone tuning $r_S$ through
it is tuning nothing.

## `sigma_g` is a units error

The harness injects gyro white noise with per-sample standard deviation
$\sigma_{\rm g,sample} = 0.00157$ rad/s at 200 Hz, and passes
$2\sigma_{\rm g,sample} = 0.00314$ into the constructor argument named
`gyro_noise_density_rad_sqrt_s`. The filter uses it as a density:

```cpp
Q_AA = Qbase * Ts;      // Qbase is variance per second
```

A per-sample standard deviation and an angle-random-walk density are not the
same quantity. The density that reproduces the injected noise is

$$\sigma_{\rm ARW} = \sigma_{\rm g,sample}\sqrt{\Delta t}
  = 0.00157\sqrt{1/200} = 1.11\times10^{-4}\ \mathrm{rad}/\sqrt{\mathrm{s}},$$

so the deployed value overstates it by $\sqrt{200} = 14.1$ on top of its own
$2\times$ inflation — $28.3\times$ in standard deviation, $800\times$ in
variance. The filter was told its gyro was two orders of magnitude noisier in
power than it is, and it duly stopped trusting it.

`Kalman3D_Wave_OU_III` already ships the conversion, as a public static helper
that no caller used:

```cpp
static Vector3 gyro_noise_density_from_sample_std(const Vector3& sample_std_rad_s,
                                                  T sample_period_s);
```

The same pattern appears in `sensors/full_marine_ins/atomS3R_ins_kalman_ou3`,
which multiplies a reference `gyr_sigma_ref_rps` by a rate scale and passes it
straight in. **That sketch is unchanged by this note and carries the same
mis-specification.**

### What the sweep says

Paired change against the deployed point, percent, negative better. `*` marks
one sign on all 40 pairs.

| `SF_SIGMA_G_SCALE` | roll | pitch | yaw | 3D | gyro bias |
| --- | --- | --- | --- | --- | --- |
| 0.02 | −11.75 | −14.32 | −0.59 | −28.25* | −42.50 |
| 0.035 | −11.55 | −14.14 | −0.64 | −33.68* | −41.41 |
| 0.056 | −11.12 | −13.92 | −0.70 | −36.74* | −39.36 |
| 0.1 | −10.28 | −13.14 | −0.79 | −35.64* | −34.35 |
| 0.32 | −7.53* | −9.37 | −0.79 | −20.32* | −13.91 |
| 1 (deployed) | — | — | — | — | — |
| 2 | +10.12* | +17.39 | +1.71 | +16.25* | +0.49 |

Roll and pitch improve on all eight records — roll $-10.1$ to $-12.9\%$, pitch
$-9.6$ to $-22.9\%$ — and 3D displacement improves on every pair. The optimum
is broad and sits at 0.035 to 0.056; the unit-consistent value is 0.0354.
**Correcting the units is the measured optimum**, to within the width of the
basin, which is the strongest form this result could take: the sweep is not
fitting a constant to a simulator, it is recovering a physical one.

Yaw is flat across the whole axis. This is a tilt effect, not a heading one,
which is what a gyro-trust change should be: heading is anchored by the
magnetometer.

## The other five

**`sigma_m` wants to go up.** At $2\times$: roll $-6.5\%$, pitch $-2.4\%$, 3D
$-1.0\%$. At $4\times$ roll continues to improve but 3D turns around
($+6.4\%$), and at $8\times$ everything degrades. The harness injects hard-iron
offset, soft-iron misalignment and scale error alongside the white noise, none
of which the $1.2\times$ inflation accounts for, so a magnetometer variance
above the white-noise floor is the physically expected answer.

**`sigma_a` is flat for attitude and asymmetric for displacement.** Below 1 it
does almost nothing ($-0.1$ to $-0.3\%$ on roll); above 1 it costs vertical
displacement steeply ($3\times \Rightarrow$ z $+14.7\%$). OU-III carries the
wave acceleration as a state, so `Racc` is meant to be sensor noise, and
inflating it lets the accelerometer stop constraining $a_w$.

**`Pb0` and `b0` do not support a move.** `Pb0` at 1e-8 buys pitch $-4.4\%$ but
costs yaw $+2.2\%$ and doubles gyro-bias error ($+105\%$); `b0` at 1e-9 buys
pitch $-2.2\%$ and gyro bias $-35.9\%$ but costs roll and 3D. Both are noisy
and neither is one-signed. They are left alone. Note that `b0 = 1e-11` is
$10\times$ below the injected bias random walk of $10^{-10}$ (rad/s)²/s, so
the mis-specification is real — it is simply not worth anything.

## The adopted point

A joint refine round over the three live axes, 440 runs:

| arm | roll | pitch | yaw | z | 3D | acc bias | gyro bias |
| --- | --- | --- | --- | --- | --- | --- | --- |
| g0.05 | −11.25 | −13.99 | −0.69 | −1.90 | −36.27* | −0.28 | −39.98 |
| g0.05_m2 | −13.13 | −14.50 | −1.19 | −1.97 | −37.08* | −0.94 | −36.46 |
| g0.07_m2 | −12.80 | −14.22 | −1.22 | −1.96 | −39.14* | −1.06 | −34.63 |
| **g0.05_m2_a0.71** | **−13.33** | **−14.65** | **−1.19** | **−2.37** | **−37.16*** | **−1.04** | **−36.49** |
| g0.05_m2_a0.5 | −13.56 | −14.51 | −1.17 | −2.26 | −37.16* | −1.00 | −36.49 |

The adopted triple, in `FusionAdapter_OU_III`:

```cpp
static constexpr float SIGMA_A_RESCALE = 0.71f;  // 2.8x -> 2.0x injected accel white
static constexpr float SIGMA_G_RESCALE = 0.05f;  // 2.0x sample std -> sqrt(2)x density
static constexpr float SIGMA_M_RESCALE = 2.0f;   // 1.2x -> 2.4x injected mag white
```

`sigma_g` at 0.05 restores the density and keeps a $\sqrt{2}$ inflation over
it. The other two are empirical, measured against this harness's noise model
rather than derived, and both are small next to the gyro term.

`tau_applied` and `sigma_applied` are **bit-for-bit unchanged** across every
arm of both rounds. The OU schedule is untouched; this is a MEKF-side effect
only. That also means the result is orthogonal to the $r_S$, anisotropy and
sigma-horizon studies, none of which needs revisiting.

## The gates

`tools/ou_regauge_gates.py --family ou_iii`, default seeds, same rule as every
previous cut — worst across the eight records plus about half a percent:

| gate | was | now | worst was | worst now |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 4.735 | 4.72 | 4.7106 | 4.6961 |
| Z %Hs PM-Stokes | 4.682 | 4.666 | 4.6580 | 4.6426 |
| yaw deg | 1.297 | 1.27 | 1.2896 | 1.2630 |
| roll deg | 0.42 | 0.3637 | 0.4179 | 0.3618 |
| pitch deg | 0.223 | 0.195 | 0.2218 | 0.1940 |
| 3D % JONSWAP | 20.95 | 13.94 | 20.8367 | 13.8686 |
| 3D % PM-Stokes | 20.86 | 14.92 | 20.7483 | 14.8387 |
| acc Z bias % | 4.624 | 4.475 | 4.6004 | 4.4519 |
| acc 3D bias % | 81.84 | 78.61 | 81.4268 | 78.2145 |

All nine tighten. The two 3D displacement bars come down by about a third,
the largest single move any of these has made. Previous re-gauges in this
family have had to argue whether a moved sentinel was a regression or a
re-draw; that question does not arise here, because nothing got worse.

## What this does not cover

- **The AtomS3R sketch is unchanged.** It carries the same `sigma_g` units
  mis-specification. Correcting it is a change to deployed device firmware
  against a real IMU whose noise is not the harness's, and it should be
  measured on hardware rather than inherited from this sweep. The
  recommendation is the same conversion:
  `gyro_noise_density_from_sample_std(gyr_sigma_ref_rps * imu_rate_scale, 1/LOOP_HZ)`,
  times whatever inflation a hardware measurement supports.
- **OU-II and TFG are unchanged.** Both take the same harness sigmas through
  the same shared `W3dSimCommon.h`, so both are very likely carrying the same
  gyro units error. Neither was swept here, and the shared harness constants
  were deliberately left alone so that this change moves one family only.
- **The committed evidence bundles are now stale.** `reports/results/ou_validation`,
  `reports/results/ou_robustness` and the generated `doc/kalman_ou_iii/*.tex-part`
  macros were produced by the previous operating point. They are regenerated
  by the `ou-validation` workflow on `main`.

## Reproducing

```sh
make -C tests/kalman_ou_iii build

# Coordinate sweep, all seven axes (1440 runs).
python3 tools/ou_qmekf_variance_study.py --seeds 1,2,3,4,5

# Joint refine over the three live axes (440 runs).
python3 tools/ou_qmekf_variance_study.py --seeds 1,2,3,4,5 --arms \
    'g0.05:SF_SIGMA_G_SCALE=0.05' \
    'g0.05_m2:SF_SIGMA_G_SCALE=0.05,SF_SIGMA_M_SCALE=2' \
    'g0.05_m2_a0.71:SF_SIGMA_G_SCALE=0.05,SF_SIGMA_M_SCALE=2,SF_SIGMA_A_SCALE=0.71'

# The gates.
python3 tools/ou_regauge_gates.py --family ou_iii
```

The scale overrides multiply the adopted constants, so a scale of 1 now
reproduces the adopted point and a re-run re-centres on it. To recover the
previous operating point exactly:

```sh
SF_SIGMA_A_SCALE=1.4084507 SF_SIGMA_G_SCALE=20 SF_SIGMA_M_SCALE=0.5
```

Raw runs: `reports/results/ou_qmekf_variances/`.
