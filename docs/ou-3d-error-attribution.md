# Where OU-III's 3D displacement error comes from

The paired study in [`ou-validation.md`](ou-validation.md) reports that OU-III
lowers vertical RMS in every scenario while raising 3D RMS in four of the five
primary ones, and that `disp_y_rms_m` is higher in all nine scenarios. This note
attributes that result to a specific band and a specific mechanism. It is a
diagnosis of the committed filters, not a proposal to change them.

The companion question -- why the *vertical* error is a larger fraction of the
sea on the small records than on the large ones, in all three families -- is
measured in [`ou-low-sea-percent-error.md`](ou-low-sea-percent-error.md).

## What the committed bundle already says

From `reports/results/ou_validation/ou_validation_paired_effects.csv`
(`OU_III_minus_OU_II`, `Adaptive`, ten paired seed triplets):

- `disp_z_rms_m` is lower for OU-III in all nine scenarios;
- `disp_y_rms_m` is higher in all nine, with intervals excluding zero;
- `disp_x_rms_m` is mixed — higher in the two smallest and the $H_s=4.0$ m seas,
  lower in the $H_s=1.5$ m and $H_s=8.5$ m seas;
- `disp_3d_rms_m` is therefore higher everywhere except the $H_s=8.5$ m seas.

The channel that decides the 3D verdict is horizontal, and the axis that always
loses is the one carrying the smaller share of the sea's horizontal motion:
every record is generated at $\pm 30^\circ$, so world $y$ holds roughly 60% of
the horizontal displacement that world $x$ does.

## Method

The deterministic single-realization simulators reproduce the effect without the
surrogate machinery, which makes it cheap to instrument:

```bash
make fetch-sim-data
W3D_VALIDATION_WINDOW_SEC=900 tests/kalman_ou_ii/kalman_ou_ii-sim   --input wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv
W3D_VALIDATION_WINDOW_SEC=900 tests/kalman_ou_iii/kalman_ou_iii-sim --input wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv
```

Each run writes a `w3d_*_fusion_ou{2,3}.csv` time series holding `disp_ref_*`
and `disp_est_*`. Over the same trailing 900 s the study scores, the reference
and estimate were compared by Welch cross-spectrum (300 s Hann segments, 50%
overlap), giving a complex displacement transfer function $H(f)$ per axis. Error
power was then split at the bottom edge of the wave band, so that

- **wave band** collects the error the filter makes *on the wave*, which is
  $|H-1|$ times the reference amplitude plus whatever is incoherent with it;
- **sub-band** collects error at frequencies where the reference carries no
  energy at all. On these records the reference displacement spectrum is empty
  below the band to five decimal places, so everything scored there is drift the
  filter invented.

## Result: the vertical gain is real, and it is paid for below the band

JONSWAP, adaptive mode, trailing 900 s, displacement RMS in metres. `sub` is the
band below the wave band (cut at 0.20, 0.087 and 0.045 Hz for the three records),
`band` is the wave band.

| record | axis | band | OU-II | OU-III | delta |
|---|---|---|---|---|---|
| $H_s=0.27$ | x | sub | 0.0075 | 0.0267 | **+0.0192** |
| | x | band | 0.0359 | 0.0315 | -0.0044 |
| | y | sub | 0.0158 | 0.0421 | **+0.0263** |
| | y | band | 0.0231 | 0.0133 | -0.0097 |
| | z | sub | 0.0107 | 0.0127 | +0.0021 |
| | z | band | 0.0151 | 0.0064 | -0.0086 |
| $H_s=1.5$ | x | sub | 0.0356 | 0.0850 | **+0.0494** |
| | x | band | 0.1959 | 0.1372 | -0.0587 |
| | y | sub | 0.1011 | 0.1542 | **+0.0532** |
| | y | band | 0.0969 | 0.1237 | +0.0267 |
| | z | sub | 0.0610 | 0.0556 | -0.0054 |
| | z | band | 0.0790 | 0.0393 | -0.0397 |
| $H_s=8.5$ | x | sub | 0.1758 | 0.4212 | **+0.2454** |
| | x | band | 1.0367 | 0.7587 | -0.2779 |
| | y | sub | 0.5267 | 0.6509 | **+0.1241** |
| | y | band | 0.5105 | 0.6085 | +0.0980 |
| | z | sub | 0.2509 | 0.2338 | -0.0171 |
| | z | band | 0.4325 | 0.2166 | -0.2159 |

Two patterns hold across all three seas:

1. **In the wave band OU-III is the more faithful filter.** Its coherent error
   $|H-1|$ is smaller than OU-II's on every axis of every record: vertically
   0.09-0.10 against 0.20-0.23, horizontally 0.32-0.53 against 0.58-0.65 on $x$.
   OU-II under-reports horizontal displacement by about 30% with a $30^\circ$
   phase lag near the peak; OU-III's horizontal phase error is under $15^\circ$.
2. **Below the wave band OU-III invents two to three times more horizontal
   drift**, while its vertical drift is equal or slightly better.

The 3D metric sums the two. Vertically the first term dominates and OU-III wins
by 25-50%. Horizontally the second term is large enough to cancel the first on
$x$ and to overturn it on $y$, whose in-band signal is smallest.

The wave-band gain scales with the sea ($z$ band delta is -0.026 to -0.032 of
$H_s$ in all three records) while the horizontal drift penalty grows more slowly
than $H_s$ ($y$ sub-band delta is 0.097, 0.035 and 0.015 of $H_s$ as $H_s$ goes
0.27, 1.5, 8.5). That is why the sign of the 3D difference flips in the largest
sea and only there, and why the flip is not evidence of a different mechanism.

## Mechanism: one integral anchor instead of a position and a velocity anchor

The families differ in how the linear block is kept from drifting.
`Kalman3D_Wave_OU_II::time_update` periodically applies two zero
pseudo-measurements, on position and on velocity:

```cpp
measurement_update_position_pseudo(Vector3::Zero(), sigma_p0);
measurement_update_velocity_pseudo(Vector3::Zero(), sigma_v0);
```

`Kalman3D_Wave_OU_III::time_update` applies one, on the extra state $S=\int p\,dt$:

```cpp
applyIntegralZeroPseudoMeas();
```

Anchoring only the third integral is a softer, higher-order high-pass than
anchoring position and velocity directly. Both consequences are visible in the
measured $|H(f)|$ on the $H_s=1.5$ m record (coherent gain, $|S_{re}|/S_{rr}$):

| band (Hz) | 0.087-0.12 | 0.12-0.175 | 0.175-0.25 | 0.25-0.35 |
|---|---|---|---|---|
| OU-II x | 0.77 | 0.69 | 0.66 | 0.68 |
| OU-II y | 1.25 | 1.01 | 1.06 | 1.02 |
| OU-II z | 0.95 | 0.96 | 0.96 | 0.96 |
| OU-III x | 2.00 | 1.39 | 1.20 | 1.09 |
| OU-III y | 2.30 | 1.58 | 1.51 | 1.32 |
| OU-III z | 1.24 | 1.13 | 1.08 | 1.04 |

OU-II's response never exceeds unity by more than a few percent. OU-III's peaks
at 2.0-2.3 at the bottom edge of the wave band and decays through it, which is
the resonance of a lightly damped high-pass whose corner the wave-band operating
point places just under the sea's own energy. The same softness that lets real
low-frequency wave motion through unattenuated — the in-band gain — lets
accelerometer-derived low-frequency error through as drift.

Why the horizontal axes and not the vertical: the vertical acceleration channel
is well observed, because the filter knows where gravity points. The horizontal
channels absorb tilt error and horizontal accelerometer bias, neither of which is
observable in this geometry — on this record the estimated body-$x$ bias is off
by 0.09 m/s$^2$ (OU-II) and 0.12 m/s$^2$ (OU-III) against a true bias of
0.046 m/s$^2$. Those errors enter as a slowly varying horizontal acceleration
whose only removal path is the drift anchor.

## Ablations

All on the JONSWAP $H_s=1.5$ m record, trailing 900 s, displacement RMS in metres.

| run | x | y | z | 3D |
|---|---|---|---|---|
| OU-II | 0.2005 | 0.1443 | 0.1008 | 0.2668 |
| OU-III | 0.1682 | 0.2032 | 0.0684 | 0.2725 |
| OU-III `OU_III_S_FACTOR=1.5` | 0.1680 | 0.2030 | 0.0685 | 0.2723 |
| OU-III `OU_III_S_FACTOR=1.0` | 0.1677 | 0.2024 | 0.0687 | 0.2717 |
| OU-III `OU_III_R_S_XY_FACTOR=0.36` | 0.2432 | 0.1976 | 0.0693 | 0.3209 |
| OU-III `OU_III_R_S_COEFF=0.1` | 0.2726 | 0.2094 | 0.0834 | 0.3537 |
| OU-III `OU_III_R_S_COEFF=1.0` | 0.2034 | 0.3155 | 0.1100 | 0.3912 |
| OU-II `OU_II_R_V0_COEFF=20` | 0.2007 | 0.1683 | 0.1026 | 0.2813 |

- **The horizontal anisotropy is not the cause.** Taking `S_factor` from the
  deployed 1.87 down to 1.0 moves $y$ by 0.4% and leaves the measured transfer
  functions unchanged to two decimal places. `sigma_aw` enters both the process
  noise and, through `r_S = 0.35 * sigma_aw * tau^3`, the anchor, so the loop
  gain that sets the corner is nearly invariant to it. (It does matter in the
  largest sea, where `S_FACTOR=1.0` takes $x$ from 1.014 to 0.882 m.) A
  five-seed sweep of both anisotropy constants together, in
  [`ou-iii-anisotropy-consistency.md`](ou-iii-anisotropy-consistency.md), puts
  that largest-sea effect at 7.7% of 3D RMS on all five seeds, so "not the
  cause" holds for this record but not for the whole envelope.
- **Retuning the single anchor cannot recover the axis.** $y$ is already near
  its own optimum at the deployed `R_S_coeff = 0.35`: tightening to 0.1 or
  loosening to 1.0 makes both $y$ and 3D worse, because a tighter anchor buys
  drift suppression by cutting into the in-band gain that OU-III exists for.
  Restoring the old `R_S_xy_factor = 0.36` buys 3% on $y$ and costs 45% on $x$.
- **Removing OU-II's velocity anchor moves OU-II toward OU-III, on exactly the
  axes and in exactly the band predicted.** With `r_v0` pushed to its clamp,
  OU-II's $y$ rises 0.144 to 0.168 while $x$ is unchanged; spectrally, sub-band
  $x$ goes 0.036 to 0.050 and sub-band $y$ 0.101 to 0.119, while both vertical
  terms move by less than 0.003.

Two secondary observations, recorded because they bound the claim:

- With `--no-noise` the ranking reverses: OU-III scores 3D 0.186 against OU-II's
  0.233 on the same record, and its $y$ is still the worse channel (0.130 against
  0.096). The horizontal drift term needs the sensor error to feed it; the
  in-band advantage does not.
- OU-III's incoherent in-band horizontal error is also larger (0.083 against
  0.046 m on $x$), consistent with the same softer loop passing more
  accelerometer noise.

## What this does and does not establish

It establishes that the 3D result is a band trade and not a tuning slip: OU-III
is the better in-band estimator on all three axes and pays for it with
sub-wave-band horizontal drift that its single integral anchor does not remove.
It does not establish that the trade is necessary. The obvious untested question
is whether adding a velocity-zero pseudo-measurement to OU-III alongside the
integral one recovers the horizontal axes without giving back the vertical gain;
nothing here has tried it, and the ablations above only price the knobs that
already exist.

The measurements are single deterministic realizations of three of the nine
scenarios. They agree in sign and rough magnitude with the ten-seed paired
bundle, which remains the evidence for the population claim.
