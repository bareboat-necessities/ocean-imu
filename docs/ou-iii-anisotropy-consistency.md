# Is the OU-III anisotropy self-consistent?

OU-III applies its stationary acceleration scale anisotropically,
$\sigma_{aw} \to (S_\sigma\sigma_{aw}, S_\sigma\sigma_{aw}, \sigma_{aw})$ with
`S_factor = 1.87`, and its integral regularizer isotropically,
$r_S \to (\rho_{xy}r_S, \rho_{xy}r_S, r_S)$ with `R_S_xy_factor = 1`. Both live
in `SeaStateFusionFilter_OU_III::apply_ou_tune_()` and `apply_RS_tune_()`, and
both are quoted side by side in the paper
(`doc/kalman_ou_iii/w3d-fus-methods.tex-part`, "Anisotropic parameterization").

The similarity theorem the $r_S$ schedule is derived from
(`w3d-rs-scale-theorem.tex-part`) gives the natural scale of the third integral
as $\sigma_S \sim \sigma_{aw}\tau^3$, per component. If the horizontal
acceleration prior is $S_\sigma$ times the vertical one, then so is the natural
horizontal $S$ scale, and the axis-consistent regularizer is

$$\rho_{xy} = S_\sigma = 1.87,$$

not 1. The deployed pair is therefore *not* the equal-similarity point. This
note asks whether that inconsistency is the source of the horizontal penalty
documented in [`ou-3d-error-attribution.md`](ou-3d-error-attribution.md), and
measures it.

## The inconsistency is real, and larger than the one that was removed

Three separate readings agree on its size and sign.

- **Similarity scales.** With $\rho_{xy}=1$ and $\sigma_{aw,H}=S_\sigma
  \sigma_{aw,Z}$, the dimensionless regularizer $r_S/(\sigma_{aw}\tau^3)$ is
  $S_\sigma = 1.87$ times *smaller* horizontally — a tighter horizontal anchor
  than the law asks for.
- **The Riccati group.** In $r = R_S/(\sigma_a^2\tau^6)$ of
  `w3d-rs-similarity-design-guidance.tex-part`, the same mismatch is
  $S_\sigma^2 = 3.5$.
- **The regularization pole.** The pole target inverts to $\kappa^3 =
  \sqrt{q_{\rm eff}}\,\tau^3/(r_S\sqrt{T_S})$, and the filter's own reading of
  the Cubic law takes $q_{\rm eff}\propto\sigma_{aw}^2$ (the sea-state-scaled
  branch argued for in `SeaStateFusionFilter_OU_III.h`). At the common $\tau$
  and $T_S$ the two channels share, $\kappa \propto
  (\sigma_{aw}/r_S)^{1/3}$, so the horizontal normalized corner sits
  $S_\sigma^{1/3} = 1.23$ times higher than the vertical one. The horizontal
  high-pass is stronger, in the channel the error attribution already
  identifies as drift-limited.

For scale: the explicit anisotropy that was removed in an earlier round,
`R_S_xy_factor = 0.36`, was a $2.8\times$ tightening in $r_S$ units. Setting it
to 1 left an implicit $1.87\times$ tightening of the same sign in place. The
claim in [`ou-iii-adaptation-and-travel-sense.md`](ou-iii-adaptation-and-travel-sense.md)
that "the horizontal-versus-vertical asymmetry is a tuning constant
(`R_S_xy_factor`), not a property of the OU model" is therefore too strong:
with $S_\sigma = 1.87$ the OU model itself is anisotropic and the regularizer
does not follow it.

What is *not* true is that $\rho_{xy}=1$ was picked on a vertical metric: the
$0.36 \to 1$ change was justified on 3D RMS, with vertical error noted as flat
across the sweep. The gap is that the sweep only ever went downward from 1,
because the setter would not accept anything above it.

## The ablation could not previously be run

`setRSXYFactor()` clamped its argument to $[0,1]$:

```cpp
R_S_xy_factor_ = std::min(std::max(k, 0.0f), 1.0f);   // before
```

The ceiling of 1 encoded the assumption that the horizontal anchor can only be
tighter than the vertical one, which is precisely the assumption in question.
`OU_III_R_S_XY_FACTOR=1.87` was accepted, silently reduced to 1, and returned
the deployed run bit for bit:

```
rho=1.0   disp_3d_rms_m=0.268828601
rho=1.87  disp_3d_rms_m=0.268828601      # same clamp, same result
```

Any earlier attempt at this experiment would have measured nothing and read as
a clean null. The bound is now 4; the default is unchanged, and
`iss_contract-test` pins both the default and the reachability of $\rho_{xy}>1$.

## Experiment

`tests/kalman_ou_iii/kalman_ou_iii-sim`, four JONSWAP records, five IMU/init
seeds each (`W3D_SEED=1..5`), trailing 900 s, adaptive tuning. Seven arms; each
is a paired per-seed comparison against the deployed point. `A` is deployed;
`B` is the axis-consistent point at the deployed $S_\sigma$; `C` and `F` are
the axis-consistent points at $S_\sigma=1$ and $0.8$; `D`, `E`, `G` separate the
two knobs.

| arm | $S_\sigma$ | $\rho_{xy}$ | consistent? |
|---|---|---|---|
| A | 1.87 | 1.00 | no (deployed) |
| B | 1.87 | 1.87 | yes |
| C | 1.00 | 1.00 | yes |
| D | 1.00 | 1.87 | no |
| E | 0.80 | 1.00 | no |
| F | 0.80 | 0.80 | yes |
| G | 1.87 | 0.53 | no (anti-law) |

Mean over five seeds of the per-seed percent change against arm A; negative is
better. `*` marks a delta with the same sign on all five seeds.

| record | arm | x | y | z | 3D |
|---|---|---|---|---|---|
| $H_s=0.27$ | B | -3.0\* | +38.3\* | -0.0 | +15.7\* |
| | C | +0.4\* | -0.3\* | +0.2\* | +0.1\* |
| | E | +0.7\* | -0.5\* | +0.4\* | +0.2\* |
| | F | +6.0\* | -9.4\* | +0.5\* | -0.2 |
| | G | +21.2\* | -18.4\* | +0.1 | +6.0\* |
| $H_s=1.50$ | B | -4.1\* | +20.1\* | -0.2 | +11.3\* |
| | C | +0.1 | -0.5\* | +0.4\* | -0.3 |
| | E | +0.8 | -0.7\* | +0.9\* | -0.2 |
| | F | +8.1\* | -3.3\* | +1.1\* | +1.0\* |
| | G | +27.4\* | -1.0 | +1.1\* | +9.7\* |
| $H_s=4.00$ | B | -7.0\* | +34.4\* | -0.3 | +7.5\* |
| | C | -3.2\* | -1.3\* | +1.1 | **-2.5\*** |
| | E | -1.8 | -2.1\* | +2.9\* | -1.8\* |
| | F | +5.7\* | -10.0\* | +3.3\* | +1.1 |
| | G | +24.3\* | -16.6\* | +1.7\* | +12.6\* |
| $H_s=8.50$ | B | +23.3\* | +17.8\* | -0.4 | +20.1\* |
| | C | -14.1\* | -2.2\* | -1.3 | **-7.7\*** |
| | E | -3.4 | -3.6\* | +0.3 | -3.8 |
| | F | +1.3 | -5.6\* | +0.5 | -2.5 |
| | G | +13.6\* | +0.9 | +1.7\* | +6.6\* |

Pooled over all four records:

| arm | x | y | z | 3D |
|---|---|---|---|---|
| B $(1.87,1.87)$ | +2.3 | +27.6 | -0.2 | **+13.7** |
| C $(1.00,1.00)$ | -4.2 | -1.1 | +0.1 | **-2.6** |
| D $(1.00,1.87)$ | -12.8 | +26.8 | -0.1 | +6.4 |
| E $(0.80,1.00)$ | -0.9 | -1.7 | +1.1 | -1.4 |
| F $(0.80,0.80)$ | +5.3 | -7.1 | +1.3 | -0.2 |
| G $(1.87,0.53)$ | +21.6 | -8.8 | +1.1 | +8.7 |

## Result: the asymmetry is real, the proposed repair is the wrong one

**Enforcing $\rho_{xy}=S_\sigma$ is worse, at both $\sigma$ anisotropies
tested.** At the deployed $S_\sigma=1.87$, closing the inconsistency costs
+13.7% of 3D RMS pooled and is worse on every record on all five seeds. At
$S_\sigma=0.8$ it costs 1.2 points against the same-$S_\sigma$ isotropic arm
(F vs E). The similarity law is not the right tuning rule for this parameter.

**$\rho_{xy}$ is an $x$-versus-$y$ trade, not a horizontal-versus-vertical
control.** Loosening it helps $x$ and hurts $y$ (B, D); tightening it does the
reverse (G); 3D loses either way. The records are generated at $\pm30^\circ$, so
world $x$ carries the larger share of horizontal motion and is signal-limited
while $y$ is drift-limited — one scalar cannot serve both, which is the same
conclusion the earlier $\rho_{xy}=0.36$ sweep reached from the other side.
Every arm moves $z$ by at most 1.7%, so the vertical metric cannot referee this
parameter at all — but the horizontal metrics, which can, also prefer 1.

**The mis-specification worth fixing is on the $\sigma_{aw}$ side.** Arm C
($S_\sigma=1$, still self-consistent, and matching the records' measured
$\sigma_{ax}/\sigma_{az}=0.81$, $\sigma_{ay}/\sigma_{az}=0.55$, combined 0.99)
takes 2.6% off pooled 3D RMS, 7.7% in the $H_s=8.5$ m sea and 14.1% off $x$
there, on all five seeds, at a vertical cost of +0.1%. Arm E ($S_\sigma=0.8$)
is directionally the same and smaller. This sharpens the earlier reading in
`ou-3d-error-attribution.md` that "the horizontal anisotropy is not the cause":
it is not the cause of the $H_s=1.5$ m result, where it is worth 0.3%, but it
is worth 7.7% of 3D in the largest sea.

**Why consistency loses.** The similarity law makes the *prior* self-similar
across axes; it does not make the *disturbance* self-similar. The horizontal
channels absorb tilt error and unobservable horizontal accelerometer bias whose
low-frequency content does not scale with the sea's horizontal acceleration
amplitude, so the drift they must reject is not $S_\sigma$ times the vertical
one. Holding $r_S$ isotropic while $\sigma_{aw,H}$ is inflated keeps the
horizontal corner 23% above the vertical one, and the measurement says the
horizontal channel wants at least that much extra high-pass.

## What was and was not changed

Changed: the `setRSXYFactor()` upper bound, 1 to 4, so the configuration exists;
and `iss_contract-test`, which now pins the deployed pair and the reachability
of $\rho_{xy}>1$.

Not changed: `S_factor` stays at 1.87 and `R_S_xy_factor` at 1. Taking
$S_\sigma$ to 1 is a deployed-behaviour change that invalidates both committed
evidence bundles, the eight-record quality gates, the ISS audit constant, and
the anisotropy figures in the paper. The measurement above is four records at
five seeds, not the ten-seed paired bundle, and it is a recommendation for the
next re-gauge rather than a re-tune. The natural version of that change is not
a smaller scalar but the per-axis stationary covariance already flagged as
Stage D in `ou-iii-adaptation-and-travel-sense.md`: $\sigma_{ax},\sigma_{ay}$
estimated separately and fed through `set_aw_stationary_cov_full`, which would
make $S_\sigma$ an estimate rather than a constant.

## OU-II has the same gap, unmeasured here

`SeaStateFusionFilter_OU_II` applies `P_factor = 1.5` to the horizontal
stationary acceleration and derives both its anchors from the same vertical
$\sigma$: $r_{p0}=c\,\sigma\tau^2$ and $r_{v0}=c\,\sigma\tau$. Axis
consistency would put both horizontal anchors at $1.5\times$ the vertical ones.
`R_p0_xy_factor` is 1 and is clamped to $[0,1]$ exactly as OU-III's was, and
$r_{v0}$ has no horizontal factor at all (`Vector3f::Constant`). Given the
OU-III result — consistency worse, the $\sigma$ side the one that matters — the
OU-II experiment worth running is `P_factor` $1.5 \to 1.0$, not the anchor
factors. Neither has been run; OU-II's setter was left clamped.

## Reproducing

```bash
make fetch-sim-data
make -C tests/kalman_ou_iii build
python3 tools/ou_anisotropy_ablation.py --seeds 1,2,3,4,5
```

That is the driver the tables above come from; it writes the per-run rows to
`reports/results/ou_anisotropy/ou_anisotropy_raw.csv` and prints the paired
summary. A single arm by hand:

```bash
cd tests/kalman_ou_iii
OU_III_S_FACTOR=1.87 OU_III_R_S_XY_FACTOR=1.87 W3D_SEED=1 \
  W3D_WRITE_TIMESERIES=0 W3D_VALIDATION_WINDOW_SEC=900 \
  ./kalman_ou_iii-sim --input wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv
```

Non-zero exit on non-default seeds only means a deployed quality gate tripped;
the `VALIDATION_METRICS` line is still emitted and is what the tables above
score. `gate` is carried in the raw CSV for that reason.

This is a four-record, five-seed deterministic sweep, not the ten-seed paired
bundle in `reports/results/ou_validation`, and it changes no default, so the
committed evidence is untouched.
