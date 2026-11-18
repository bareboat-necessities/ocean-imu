# Transferring the horizontal TFG wave constraint to CARMA

Both TFG-CARMA(4,1) filters can express the OU horizontal wave constraint using their existing states. The two-axis correction improves pooled accelerometer-bias RMS but raises pooled displacement and yaw RMS. Large waves encountered bow-on account for the largest displacement regressions.

## Direct transfer

Each filter ran on eight JONSWAP and PM–Stokes records per sea direction at four sensor-noise and initialization draws. The rows pool the last 900 seconds of each record by taking the square root of the mean squared record RMS. Every comparison uses the same wave samples, draw and gyro profile in both arms.

| Filter | Matched records | 3D displacement RMS (m) | 3D accelerometer-bias RMS (m/s²) | Yaw RMS (°) |
| --- | ---: | ---: | ---: | ---: |
| TFG-CARMA, independent axes | 64 | 0.22307 → 0.22881 (+2.58%) | 0.03771 → 0.03593 (-4.71%) | 0.7419 → 0.7571 (+2.04%) |
| TFG-CARMA, shared modes | 64 | 0.20258 → 0.21069 (+4.00%) | 0.03744 → 0.03573 (-4.55%) | 0.7367 → 0.7513 (+1.99%) |

| Filter | Following-sea 3D RMS (m) | Head-sea 3D RMS (m) | Records with lower 3D RMS |
| --- | ---: | ---: | ---: |
| TFG-CARMA, independent axes | 0.22261 → 0.21864 (-1.78%) | 0.22353 → 0.23855 (+6.72%) | 38/64 |
| TFG-CARMA, shared modes | 0.19165 → 0.18909 (-1.33%) | 0.21296 → 0.23026 (+8.13%) | 35/64 |

## Direction and uncertainty

A separate experiment projected the correction across the propagation axis estimated by the onboard wave-direction tracker. It skipped unreliable axes and needed no propagation sign. The spectrum-weighted version also allowed for genuine broadband wave motion in the projected residual.

| Filter | Horizontal observation | Following-sea 3D change | Head-sea 3D change | Pooled 3D change | Pooled bias change | Pooled yaw change | Fixed-draw default records outside existing limits |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| TFG-CARMA, independent axes | Both axes | -1.78% | +6.72% | +2.58% | -4.71% | +2.04% | 9/16 |
| TFG-CARMA, independent axes | Across measured wave axis | -2.05% | +0.49% | -0.77% | -1.88% | -0.65% | 3/16 |
| TFG-CARMA, independent axes | Across measured wave axis, spectrum-weighted | -1.43% | +0.07% | -0.67% | -1.19% | -0.41% | 1/16 |
| TFG-CARMA, shared modes | Both axes | -1.33% | +8.13% | +4.00% | -4.55% | +1.99% | 10/16 |
| TFG-CARMA, shared modes | Across measured wave axis | -2.47% | +0.55% | -0.79% | -1.85% | -0.64% | 4/16 |
| TFG-CARMA, shared modes | Across measured wave axis, spectrum-weighted | -1.75% | +0.09% | -0.73% | -1.17% | -0.40% | 2/16 |

The archive includes fourfold larger observation variance, along-wave projection, spectrum weighting alone, motion-only gain, motion-plus-accelerometer-bias gain and an accepted-startup tilt uncertainty of 0.5°. The fixed-draw comparison also covers the interval-gyro regression profile. The original per-record limits were used throughout.

The existing CARMA defaults and replay limits were retained. The largest failure that persisted under the spectrum-weighted cross-wave correction occurred in the fixed-draw 8.5 m PM–Stokes head sea.

| Filter | Original Y RMS (m) | Direct-transfer Y RMS (m) | Across-wave, spectrum-weighted Y RMS (m) | Existing Y limit (m) |
| --- | ---: | ---: | ---: | ---: |
| TFG-CARMA, independent axes | 0.29204 | 0.31202 | 0.29773 | 0.29351 |
| TFG-CARMA, shared modes | 0.28691 | 0.30629 | 0.29219 | 0.28835 |

In the 8.5 m following PM–Stokes sea, the shared-mode candidate also exceeded the heave limit. Its Z RMS rose from 0.18692 m to 0.18790 m while its 3D RMS fell from 0.35002 m to 0.34391 m.

## Why CARMA can represent the constraint

TFG-OU already applies a pseudo-observation to the displacement integral $S$. Its horizontal change replaces the target $S=0$ with $S+\beta v=0$, where $\beta=(T_{\mathrm{enc}}/2\pi)^2$. The CARMA filters have stable displacement dynamics and add this as an extra observation. Their shipped CARMA(4,1) shape has a spectral zero at DC. If $p$ is displacement, $v$ velocity, $a$ acceleration and $\ell$ the fourth CARMA state, the model gives

$$\dot\ell=\omega_1^2p+2\zeta_1\omega_1v+a,\qquad S_*=(\ell-2\zeta_1\omega_1p-v)/\omega_1^2.$$

At fixed coefficients, $\dot S_*=p$. Shared-mode coloring preserves that state equation. The experimental residual is the horizontal part of $S_*+\beta v$. A frequency retune changes the coefficients defining $S_*$, so the trial uses a locally defined model quantity. Broadband motion can produce a real residual even when the filter state is correct.

The direct arm uses the filter’s tuned encounter period, its existing pseudo-update cadence, and 5/8 of its horizontal observation variance. It starts after live handover and accelerometer-bias release. Candidate measurements use the common right-invariant correction and covariance update. The directional projection accounts for how rotating the vessel changes the estimated world-frame axis.

For the spectrum-weighted arm, the transfer function from white process noise to $S_*+\beta v$ is $(1+\beta s^2)/(D_1D_2)$ between coefficient retunes. Its two-sided density at zero is $q/(\omega_1^4\omega_2^4)$. Dividing this density by the observation interval gives a low-frequency equivalent discrete variance. Each axis uses the model’s own horizontal noise intensity. The archive includes controls for this assumption.

## Reproduction

Download [study.tar.gz](study.tar.gz). It contains each tested source patch, the input and binary digests, all three observation windows, the late per-axis bias percentages and a single-record reproduction script. The manifest maps each arm to its patch. Apply each patch to a separate clean checkout. The unmodified control matched all 576 selected values from the two hosted fixed-draw gyro profiles. [Hosted run](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44299828783)

The example below repeats the large head-sea PM–Stokes direct-transfer run. Set the DATA environment variable to a directory containing the following and head v1.4.0 datasets. Use the GCC, Eigen and libm versions recorded in the manifest when comparing fixed-draw numbers.

```sh
mkdir -p /tmp/tfg-carma-study
tar -xzf study.tar.gz -C /tmp/tfg-carma-study
git worktree add --detach /tmp/tfg-carma-test 85e7284972e4dc046bdc34cbcb4f044fbbc2a4ab
git -C /tmp/tfg-carma-test apply /tmp/tfg-carma-study/base.patch
make -C /tmp/tfg-carma-test/tests/kalman_tfg_carma build
python3 /tmp/tfg-carma-study/reproduce.py --repo /tmp/tfg-carma-test --data "$DATA" \
  --family independent --arm match --case head \
  --record wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv
```

The archive contains 1,024 complete 20-minute replays. Its SHA-256 is `652381146e514bbe1a6ea5fb92c3076b1e8d50d31aee7c95363020d3a0932d4f`. The scores are simulation results for these moving-vessel wave records.

[Independent sensor-noise seeds, original logs and shipping CARMA limits](independent-seeds.md)
