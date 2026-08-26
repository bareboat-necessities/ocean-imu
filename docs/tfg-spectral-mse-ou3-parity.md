# TFG SpectralMSE and OU-III front-end parity

This note records the deterministic comparison used to move TFG's default
integral-regularization schedule to the same reduced spectral-MSE law as
OU-III while retaining the previous TFG cubic schedule as the low-cost
embedded configuration.

## Laws

The new default is

\[
r_S=C_J q_{\rm eff}^{1/14}\sigma_{a,B}^{6/7}\tau^{24/7}/\sqrt{T_S},
\]

with `C_J = 0.0538`. With the deployed self-similar pseudo-measurement cadence
`T_S proportional to tau`, the effective tau exponent is `41/14` away from
cadence clamps.

`RSAdaptationLaw::LegacyCubic` preserves TFG's previous schedule exactly:

\[
r_{S,base}=0.28\,\sigma_{aw}\tau^3,
\qquad
r_{S,filter}=r_{S,base}\sqrt{T_{S,0}/T_S}.
\]

It is selectable with `setEmbeddedFriendlyLegacyRSLaw(true)` (or `setRSLaw`).

## Front-end parity audited against current OU-III

The TFG measurement-only adaptation/startup front end now uses the same
statistical timing and startup/magnetic constants as the deployed OU-III path:

- canonical `WavePeriodEstimator` state, including log-period smoothing;
- wave-period value used as soon as finite, with the same 0.2 Hz prior;
- `tau`/`sigma` EMA horizon `0.40 * T_sea`;
- `r_S` EMA horizon `1.5 * tau_target`;
- acceleration-variance horizon of four sea periods;
- period-scaled acceleration variance band, 0.5--4.0 times tuning frequency;
- startup Mahony gains `two_kp=0.2`, `two_ki=0.02`;
- world-frame gravity-trust average: 12 s horizon, 5 s warmup, 0.075 sine gate,
  2 s sustained hold;
- deployed wrapper tuner warmup of 10 s;
- magnetometer delay 7 s, acquisition minimum 128 samples / 15 s;
- magnetic refinement at 90 s for a 30 s window;
- continuous hard-iron memory 600 s, absolute ridge `5e-4`, relative ridge
  `0.5`, minimum information `2`, minimum effective weight `500`, maximum
  residual `3 uT`, maximum bias fraction `0.35`, full fitted fraction, 45 s
  slew.

TFG-specific estimator geometry and independently measured horizontal priors
remain TFG-specific; this change is parity of the shared front end, not an
attempt to make the TFG MEKF numerically identical to OU-III.

## Three-arm deterministic comparison

All arms use the same eight JONSWAP/PM-Stokes records, the same sensor-error
seeds (`1234/5678/9012`), and the same final 900 s scoring window.

- **main**: exact pre-PR `main` at `6d7fe815b1f218cce0007e481fb65c3f0d315f8a`.
- **parity + legacy**: the new front end with `LegacyCubic`; isolates front-end
  changes from the new law.
- **parity + SpectralMSE**: new front end and new default law.

| Record | Z %Hs main | Z %Hs legacy | Z %Hs MSE | 3D %refmax main | 3D legacy | 3D MSE |
|---|---:|---:|---:|---:|---:|---:|
| JONSWAP 0.27 | 4.7752 | 4.7392 | 4.6736 | 19.2179 | 19.1431 | 18.6272 |
| JONSWAP 1.50 | 4.3529 | 4.3716 | 4.3740 | 20.3271 | 20.1965 | 20.3211 |
| JONSWAP 4.00 | 4.2300 | 4.2855 | 4.2635 | 18.1782 | 18.0887 | 18.1110 |
| JONSWAP 8.50 | 3.9758 | 3.8723 | 3.8076 | 16.2425 | 16.1373 | 15.9877 |
| PM/Stokes 0.27 | 4.6862 | 4.6662 | 4.6077 | 19.5676 | 19.4935 | 18.9993 |
| PM/Stokes 1.50 | 4.0404 | 4.0669 | 4.1217 | 17.6517 | 17.7040 | 17.9667 |
| PM/Stokes 4.00 | 4.4139 | 4.3704 | 4.3963 | 19.2831 | 19.2930 | 19.5169 |
| PM/Stokes 8.50 | 4.4519 | 4.0300 | 3.9841 | 19.2593 | 19.9851 | 20.0491 |
| **mean** | **4.3658** | **4.3002** | **4.2786** | **18.7159** | **18.7552** | **18.6974** |

Relative to old `main`, SpectralMSE + parity improves mean vertical RMS by
about **2.00%** and mean 3-D RMS by about **0.10%**. The worst vertical value
falls from 4.7752 to 4.6736 %Hs. Mean accelerometer-bias 3-D error also falls
from 91.02% to 86.05% of the maximum true bias magnitude.

Comparing only the two parity arms isolates the law: SpectralMSE improves mean
vertical RMS by about **0.50%** and mean 3-D RMS by about **0.31%** relative to
LegacyCubic. Thus the new law gives a small but measurable displacement gain;
the larger per-record changes are predominantly caused by the requested
front-end/startup parity work.

The PM/Stokes Hs=8.5 m 3-D result is the important attribution example. Old
`main` gives 19.2593%, parity+legacy gives 19.9851%, and SpectralMSE gives
20.0491%. Almost all of that record's increase therefore comes from front-end
parity; changing the law adds only 0.064 percentage point. On the same record,
vertical RMS improves 4.4519 -> 4.0300 -> 3.9841 %Hs, yaw improves
1.2008 -> 0.5212 -> 0.5032 deg, and accelerometer-bias 3-D error improves
136.18 -> 63.53 -> 60.80%.

The comparison was produced by `tfg-validation` workflow run `32760913213`.
The workflow builds and runs all three arms and uploads their complete logs as
`tfg-rs-law-comparison`.

## Re-gauged deterministic gates

Both copies of TFG's eight-record quality gate — the `kRegressionBars` block in
`tests/kalman_tfg/kalman_tfg-sim.cpp` and the `Gate SpectralMSE default` step of
`.github/workflows/tfg-validation.yml` — are fitted by the rule in
`docs/quality-gate-regauge.md`: the worst value across the eight scored records
plus about half a percent, rounded up at the quantum the channel is quoted in.
The workflow copy was re-derived when this change landed and the in-binary copy
was not, so `tests/kalman_tfg/run_tests.sh` failed on the last record — PM/Stokes
3-D at 20.0469 against a bar of 19.64 — on every run after the merge. Both now
carry the same seven numbers, from `tools/ou_regauge_gates.py --family tfg`:

| gate | was | now | worst observed | margin |
| --- | --- | --- | --- | --- |
| Z %Hs JONSWAP | 4.812 | **4.698** | 4.6736 (jonswap H0.27) | 0.52% |
| Z %Hs PM-Stokes | 4.71 | **4.631** | 4.6077 (pmstokes H0.27) | 0.51% |
| yaw deg | 1.352 | **1.292** | 1.2851 (jonswap H1.50) | 0.54% |
| 3D % JONSWAP | 20.43 | **20.43** | 20.3211 (jonswap H1.50) | 0.54% |
| 3D % PM-Stokes | 19.64 | **20.15** | 20.0469 (pmstokes H8.50) | 0.51% |
| acc Z bias % | 5.12 | **4.532** | 4.5086 (jonswap H8.50) | 0.52% |
| bias 3D % | 164.9 | **155.6** | 154.786 (pmstokes H4.00, accel) | 0.53% |

Six bars come down and one is unchanged; the PM-Stokes 3-D bar is the only one
that rises, for the reason the table above this section gives — the binding
record moves from PM/Stokes 0.27 (19.5357 before, 18.9993 now) to PM/Stokes
8.50, whose rise is almost entirely front-end parity rather than the new law,
while the mean 3-D figure improves. `bias_3d_percent` gates the gyro channel
too; the accelerometer still sets it, with the gyro worst at 71.72% (PM/Stokes
8.50) against 155.6.

Rebuilding the simulator at `-march=x86-64` instead of the host's native
architecture moves the binding records by at most 4.1e-4 relative, so the
thinnest of these margins is about twelve times the spread it has to survive.
All seven pass under both builds.
