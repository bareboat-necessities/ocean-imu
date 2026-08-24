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