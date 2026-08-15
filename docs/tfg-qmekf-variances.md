# The MEKF variances TFG never swept

**Outcome: nothing is adopted, and that is the result.** TFG carries the same
`sigma_g` units error as both OU families, and correcting it improves 3D
displacement by 20–43% one-signed on every pair. But at every setting tested
the correction either regresses yaw or destabilizes pitch, and a targeted
ablation shows why: the **accelerometer-bias estimator cannot cope with a
trusted gyro** in the largest seas. The variance is not the thing to change
first.

This is the third family in the series begun by
[`ou-iii-qmekf-variances.md`](ou-iii-qmekf-variances.md) and continued in
[`ou-ii-qmekf-variances.md`](ou-ii-qmekf-variances.md). Both of those adopted
the correction and improved every channel. TFG does not, and the difference is
worth more than the improvement would have been.

## What TFG has, and does not have

TFG is a right-invariant EKF on the two-frame group of `src/lie/TwoFrameGroup.h`,
so it is shaped differently from the OU pair:

- `gyro_noise_density` was already a plain `Config` field. There is no dead
  `initialize_ext()` here — the defect that hid four variances in both OU
  wrappers is absent.
- There is no constructor seed for the pseudo-measurement regularizer: the
  tuner owns `r_S` outright. So no `R_S_noise_var` axis, and no third inert
  parameter to report.
- `P` is seeded by one isotropic `initialize_identity()` value rather than a
  split `Pq0`/`Pb0`, so there is no separate initial gyro-bias variance.

Two variances *were* unreachable, hard-coded in the filter rather than orphaned
by a wrapper: the `1e-4` seed in `initialize_identity()` and `Q_bg_`. Both are
`Config` fields now, defaults unchanged.

One thing TFG already had right: `Q_bg_` defaults to `1e-10`, which **matches**
the harness's injected bias random walk of $10^{-10}$ (rad/s)²/s exactly. Both
OU families sat at `1e-11`, ten times low. It reports `+0.00` on every channel
at its default, which also confirms the new wiring is a true no-op.

## The same units error

`W3dSimCommon.h` hands TFG `sigma_g.x()` $= 2\sigma_{\rm g,sample} = 0.00314$,
and the filter integrates it as a density:

```cpp
Qd.noalias() += wq * (Lg * Q_gyro_ * Lg.transpose());   // wq quadrature-weighted over the step
```

The injected density is $\sigma_{\rm g,sample}\sqrt{\Delta t} =
1.11\times10^{-4}$ rad/$\sqrt{\rm s}$, so the deployed value is $28.3\times$ too
large in standard deviation and $800\times$ in variance — the same error, in a
third place.

Displacement responds exactly as it did in the OU families: 3D improves at every
`sigma_g` below 1, one-signed on all 40 pairs, bottoming near $-40\%$ at 0.056.
That part is not in doubt.

## Attitude: two constraints that do not overlap

TFG's paired spreads are an order of magnitude wider than either OU family's.
Going from five seeds to eight barely narrowed them, which settles that this is
the filter's behaviour and not thin sampling. **Pooled means are unsafe to read
here**, and the two constraints below are invisible in them.

### Yaw regresses in the aggressive region

Pooled, the aggressive settings improve yaw by 5–14%. Per record, they do not:
every leading arm degrades yaw on `jonswap_H1.500` by 3.9–8.8% at three to four
standard errors, and the pooled gain is carried entirely by the large seas where
yaw improves 20–35%.

| `SF_SIGMA_G_SCALE` | yaw on `jonswap_H1.5` | 3D pooled |
| --- | --- | --- |
| 0.035 | +10.81 ± 4.0 | −40.21 |
| 0.056 | +9.18 ± 3.2 | −40.57 |
| 0.1 | +5.43 ± 2.0 | −37.36 |
| **0.18** | **−0.54 ± 3.2** | −29.45 |
| 0.32 | −5.27 ± 4.8 | −19.54 |

The magnetometer term is the counter-test for the mechanism. `sigma_m` at
$2\times$ is what *both* OU families wanted, and it makes TFG's yaw worse:
`g0.056` is $-9.30$ pooled, `g0.056_m2` is $+1.24$. Heading is
magnetometer-observable only, so trusting the gyro more while trusting the
magnetometer less removes heading authority from both sides at once. This is
the one place the three families genuinely disagree.

### Pitch and accelerometer bias blow up in the safe region

At `sigma_g = 0.18`, where yaw is clean on all eight records, pitch and the
accelerometer-bias error blow up together on the two largest seas:

| record | pitch | accel bias |
| --- | --- | --- |
| jonswap_H8.500 | **+143.8 ± 41.3** | **+100.1 ± 32.7** |
| pmstokes_H8.500 | **+197.0 ± 32.9** | **+133.7 ± 28.3** |

Same two records, same sign, similar magnitude. That pairing is the signature
of tilt/bias confusion: a pitch error tips gravity into body X and reads back
as an X bias.

## The ablation

`TFG_ACC_BIAS_UNLOCK_SEC=100000` holds the accelerometer-bias estimator shut
for the whole replay. Pitch RMS (deg), trailing 900 s:

| record | seed | deployed | dep+frozen | g0.18 | g0.18+frozen |
| --- | --- | --- | --- | --- | --- |
| pmstokes_H8.5 | 1 | 0.6332 | 0.4651 | 1.5734 | 0.4703 |
| jonswap_H8.5 | 2 | 0.2166 | 0.2859 | 0.4371 | **0.2092** |
| jonswap_H8.5 | 3 | 0.1849 | 0.3183 | 0.7328 | **0.2197** |
| pmstokes_H8.5 | 2 | 0.2111 | 0.3019 | 0.8048 | **0.2110** |
| pmstokes_H8.5 | 3 | 0.3220 | 0.3395 | 0.9535 | **0.2198** |

Three readings:

1. **The gyro correction does no pitch harm on its own.** Frozen, `g0.18` and
   the deployed variances are the same filter to within a percent.
2. **Frozen, the corrected gyro is the best configuration tested**, and by a
   wide margin the steadiest: 0.2092–0.2198 across the last four rows, a 5%
   spread, against 74% for the deployed filter and 118% for `g0.18` with the
   estimator running.
3. **The estimator is what fails, and only under a trusted gyro.** It is not
   simply harmful — at the deployed variances it helps on three of those four
   rows. Trusting the gyro leaves less attitude innovation for the bias state
   to hide in, and in the largest seas it takes the wave-correlated specific
   force instead.

## Why nothing is adopted

There is no setting that is safe on all channels:

- 0.035–0.1 regresses yaw on `jonswap_H1.500`.
- 0.18–0.32 keeps yaw clean but degrades pitch 92–197% and accelerometer bias
  70–134% on both H8.5 records.

Displacement would gain 20–43%, but it would be bought by degrading attitude
and bias in the largest seas, which is where a marine attitude estimate matters
most. `SIGMA_A_RESCALE`, `SIGMA_G_RESCALE` and `SIGMA_M_RESCALE` therefore stay
at 1.0 and the deployed filter is bit-for-bit unchanged. The seven regression
bars are untouched, because nothing moved.

## What to do instead

The correction is worth having and the ablation says what blocks it. The work
is in the accelerometer-bias gating, not the variance:

- TFG holds the bias via `acc_bias_hold_`, released by `maybeUnlockAccBias_()`
  once the magnetic reference is refined. OU-III additionally requires 250
  magnetometer updates after going Live (`acc_bias_unlock_mag_updates`),
  explicitly so the bias cannot absorb a settling tilt error. TFG has no
  equivalent count.
- OU-III also bounds the state with `set_accel_bias_limit()`. Whether TFG's
  projection is as tight under a trusted gyro is untested.
- The failure is sea-state dependent — only the H8.5 records — which suggests
  the bias should be gated on the wave-band excitation the tuner already
  estimates, rather than on elapsed time and mag updates alone.

Once the estimator holds under a trusted gyro, this sweep should be re-run; the
frozen-bias rows suggest the result would then look like the OU families', with
pitch near 0.21 and steady.

## What this does not cover

- **No adoption, so no re-gauge.** `tools/ou_regauge_gates.py --family tfg` was
  not run, because the filter it would measure is the one already shipping.
- **The AtomS3R sketches are unchanged**, as in both OU notes.
- **The shared harness constants in `W3dSimCommon.h` are untouched.** All three
  families now demonstrably carry the same units error, so fixing it there is
  finally a defensible single change — but it would move TFG, which this note
  has just shown must not move yet.
- **The ISS question is separate.** TFG's stability can be argued along the same
  four steps as OU-III's `w3d-iss-stability.tex-part`, with the attitude step
  strictly easier — its Jacobians are the constant $[g_w]_\times$ and
  $-[B_w]_\times$ rather than estimate-dependent body vectors, and
  $\Phi_{\varphi\varphi} = I$. Three gaps need work: the `Φ[V, BG]` coupling
  OU-III lacks, the `[S]_\times` attitude term in the integral residual, and
  the absence of an ISS contract test. The group-affine/log-linear claim is
  *not* available for the complete estimator; see `tfg-math-hardening.md`.

## Reproducing

```sh
make -C tests/kalman_tfg build

python3 tools/ou_qmekf_variance_study.py --family tfg --seeds 1,2,3,4,5
python3 tools/ou_qmekf_variance_study.py --family tfg --seeds 1,2,3,4,5,6,7,8 --arms \
    'g0.18:SF_SIGMA_G_SCALE=0.18' 'g0.25:SF_SIGMA_G_SCALE=0.25'

# The ablation.
cd tests/kalman_tfg
SF_SIGMA_G_SCALE=0.18 TFG_ACC_BIAS_UNLOCK_SEC=100000 W3D_SEED=2 \
  W3D_VALIDATION_WINDOW_SEC=900 W3D_WRITE_TIMESERIES=0 \
  ./kalman_tfg-sim --input wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv
```

Raw runs: `reports/results/ou_qmekf_variances/tfg_*`.
