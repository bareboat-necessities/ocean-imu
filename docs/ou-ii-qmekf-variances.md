# The MEKF variances OU-II never swept

**Outcome: the three sensor variances move, and every gated channel improves.**
Roll $-12.5\%$, pitch $-11.6\%$, yaw $-1.2\%$, vertical displacement $-1.7\%$,
3D displacement $-23.8\%$, accelerometer bias $-0.7\%$, gyroscope bias
$-36.6\%$, paired over eight records and five seeds. All nine deterministic
sentinels come down; none goes up. Three of the seven variances turn out to be
incapable of affecting anything.

This is the OU-II half of [`ou-iii-qmekf-variances.md`](ou-iii-qmekf-variances.md),
run after that study found the same constants unswept in the other family. The
two share a harness, and — as it turns out — share its errors. Read the OU-III
note first; this one records where OU-II agrees, and the two places it does
not.

## The same structure, the same dead code

`Kalman3D_Wave_OU_II` takes seven variances, differing from OU-III only in the
regularizer: `R_p0_noise_var` and `R_v0_noise_var` where OU-III has a single
`R_S_noise_var`.

`SeaStateFusion_OU_II::begin()` called the three-argument `initialize()`, so
`initialize_ext()` was unreachable and five of the seven ran on header defaults
that no caller could change. They are `Config` fields now, carrying those same
defaults, and the simulator exposes all seven — `SF_RP0_NOISE_VAR` and
`SF_RV0_NOISE_VAR` in place of OU-III's `SF_RS_NOISE_VAR`, everything else
identical.

## `sigma_g`: the same units error, the same size

The shared `W3dSimCommon.h` hands both families
$2\sigma_{\rm g,sample} = 0.00314$ for an argument named
`gyro_noise_density_rad_sqrt_s`, and `Kalman3D_Wave_OU_II` integrates it as a
density exactly as OU-III does (`Q_AA = Qbase * Ts`). The injected density is
$\sigma_{\rm g,sample}\sqrt{\Delta t} = 1.11\times10^{-4}$ rad/$\sqrt{\rm s}$,
so the deployed value is $28.3\times$ too large in standard deviation and
$800\times$ in variance. **This is one error in two places, not two errors.**

| `SF_SIGMA_G_SCALE` | roll | pitch | yaw | 3D | gyro bias |
| --- | --- | --- | --- | --- | --- |
| 0.02 | −11.39 | −11.51 | −0.78 | −13.28 | — |
| 0.035 | −11.18 | −11.42 | −0.79 | −19.20* | — |
| 0.056 | −10.79 | −11.28 | −0.82 | **−22.55*** | — |
| 0.1 | −10.07 | −10.66 | −0.90 | −22.42* | — |
| 0.32 | −7.60 | −7.71 | −0.87 | −13.52* | — |
| 1 (deployed) | — | — | — | — | — |
| 2 | +10.54* | +16.42* | +1.82 | +10.86* | — |

Same broad optimum bracketing the unit-consistent 0.0354, same flat yaw, same
one-signed 3D improvement. The attitude gain matches OU-III almost exactly
(roll $-11$ vs $-11.5$); the displacement gain is smaller ($-22.6$ vs $-36.7$),
which is what one would expect from a family whose regularizer anchors position
and velocity directly rather than through a third integral.

## Three inert parameters, not two

`Pq0` is inert here for the same reason as in OU-III. So are **both**
regularizers: `R_p0_noise_var` over four decades and `R_v0_noise_var` over four
decades leave every metric bit-for-bit identical, because `apply_ou_tune_()`
overwrites both at `enterLive_()` before the linear block runs. OU-II therefore
carries three constructor parameters that cannot affect its output, out of
seven.

As in OU-III, this is not a dead wire — `Pb0` and `b0` ride the same
`initialize_ext()` call and do move results — and it is not a bug to fix, since
the values that overwrite them are the derived ones. It is a documentation
problem, and a slightly worse one than OU-III's: two of the five numbers a
caller sees in the constructor signature after the sensor sigmas do nothing.

## Where OU-II differs: the magnetometer, and why the axis sweep lies about it

In OU-III, `sigma_m` at $2\times$ improved every channel. In OU-II the
one-at-a-time sweep says the opposite — roll $-5.1\%$ but pitch $+1.7\%$ and 3D
$+1.9\%$, with $4\times$ and $8\times$ clearly worse.

That reading is an artifact of holding `sigma_g` at its mis-specified value.
Once the gyro is corrected, the sign flips:

| arm | roll | pitch | yaw | z | 3D |
| --- | --- | --- | --- | --- | --- |
| g0.05 | −10.91 | −11.33 | −0.81 | −0.91 | −22.01* |
| g0.05_m1.25 | −11.54 | −11.37 | −0.91 | −0.93 | −23.02* |
| g0.05_m1.5 | −11.96 | −11.46 | −1.03 | −0.94 | −23.43* |
| g0.05_m2 | −12.21 | −11.70 | −1.28 | −0.96 | −23.63* |

At the corrected gyro, $2\times$ mag is better than $1\times$ on *every*
channel, including the two it hurt before. A filter that cannot trust its gyro
leans on the magnetometer, so widening the magnetometer gate costs it; once the
gyro carries its proper weight, the same widening only rejects hard-iron and
misalignment error. **This is the case for running the joint round rather than
reading the optimum off the coordinate sweep**, and it is the one place in
either family where the two disagree.

## The adopted point

Confirmation round over the three best candidates:

| arm | roll | pitch | yaw | z | 3D | gyro bias |
| --- | --- | --- | --- | --- | --- | --- |
| **g0.05_m2_a0.5** | **−12.46** | **−11.57** | **−1.22** | **−1.66** | **−23.76*** | **−36.56** |
| g0.05_m2_a0.71 | −12.34 | −11.66 | −1.25 | −1.48 | −23.73* | −36.55 |
| g0.07_m2_a0.5 | −12.18 | −11.37 | −1.24 | −1.66 | −25.70* | −34.80 |

`g0.07_m2_a0.5` buys another 2 points of 3D displacement and gives back roll
and pitch. This study is about attitude, so the adopted point is the one that
wins roll and pitch:

```cpp
static constexpr float SIGMA_A_RESCALE = 0.5f;   // 2.8x -> 1.4x injected accel white
static constexpr float SIGMA_G_RESCALE = 0.05f;  // 2.0x sample std -> sqrt(2)x density
static constexpr float SIGMA_M_RESCALE = 2.0f;   // 1.2x -> 2.4x injected mag white
```

OU-III adopted the same gyro and mag values and $0.71$ on accel rather than
$0.5$; both are on the flat part of that axis in both families, and neither
family is sensitive enough there to justify forcing them equal.

`tau_applied` and `sigma_applied` are **bit-for-bit unchanged** across every
arm of all three rounds, so the OU schedule is untouched and this result is
orthogonal to the $(r_{p0}, r_{v0})$ tuning in
[`ou-ii-pseudo-variance-tuning.md`](ou-ii-pseudo-variance-tuning.md) and to the
EMA horizons in [`ou-ema-adaptation-tuning.md`](ou-ema-adaptation-tuning.md).

## The gates

`tools/ou_regauge_gates.py --family ou_ii`, default seeds, same rule as every
previous cut:

| gate | was | now | worst was | worst now |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 6.865 | 6.776 | 6.8300 | 6.7420 |
| Z %Hs PM-Stokes | 6.848 | 6.803 | 6.8139 | 6.7688 |
| yaw deg | 1.089 | 1.074 | 1.0833 | 1.0681 |
| roll deg | 0.4792 | 0.4352 | 0.4768 | 0.4330 |
| pitch deg | 0.3657 | 0.279 | 0.3638 | 0.2775 |
| 3D % JONSWAP | 20.92 | 16.54 | 20.8140 | 16.4569 |
| 3D % PM-Stokes | 21.03 | 17.67 | 20.9203 | 17.5728 |
| acc Z bias % | 5.324 | 4.802 | 5.2969 | 4.7776 |
| acc 3D bias % | 94.47 | 92.35 | 93.9911 | 91.8873 |

All nine tighten; pitch by 24 percent and the two 3D bars by about a fifth.
All eight scored records pass and the OU-II suite exits 0.

## What this does not cover

- **The AtomS3R OU-II sketch is unchanged**, as its OU-III counterpart is. Both
  build `sigma_g` the same way and both carry the same mis-specification, and
  correcting deployed firmware wants a hardware measurement rather than this
  sweep's number. The recommended form is
  `gyro_noise_density_from_sample_std(gyr_sigma_ref_rps * imu_rate_scale, 1/LOOP_HZ)`
  times whatever inflation hardware supports.
- **TFG is unswept.** It draws the same sigmas from the same shared harness and
  is very likely carrying the same gyro units error. It was not measured here.
- **The shared harness constants in `W3dSimCommon.h` are deliberately
  untouched.** Fixing the units there would move all families at once,
  including the one nobody has measured. Each family's correction lives in its
  own adapter.
- **The committed evidence bundles are stale for both OU families.** They are
  regenerated by the `ou-validation` workflow on `main`.

## Reproducing

```sh
make -C tests/kalman_ou_ii build

python3 tools/ou_qmekf_variance_study.py --family ou_ii --seeds 1,2,3,4,5
python3 tools/ou_qmekf_variance_study.py --family ou_ii --seeds 1,2,3,4,5 --arms \
    'g0.05_m2_a0.5:SF_SIGMA_G_SCALE=0.05,SF_SIGMA_M_SCALE=2,SF_SIGMA_A_SCALE=0.5'
python3 tools/ou_regauge_gates.py --family ou_ii
```

To recover the previous operating point exactly:

```sh
SF_SIGMA_A_SCALE=2 SF_SIGMA_G_SCALE=20 SF_SIGMA_M_SCALE=0.5
```

Raw runs: `reports/results/ou_qmekf_variances/ou_ii_*.csv`.
