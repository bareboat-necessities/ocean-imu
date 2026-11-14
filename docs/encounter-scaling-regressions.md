# Encounter-period scaling and RMS regressions

All comparisons use exact-clock oceanography-waves-lib v1.3.1 inputs: the zero-speed/free-drifting 28 ft response and the 5 kn response. The moving records are following or quartering seas. They change encounter period and acceleration without a comparable change in displacement amplitude. All eight input hashes and the two archive checksums are bound in the [paired study](../reports/results/encounter_scaling/pii_phase.json).

## Physical input change

| Spectrum | Nominal Hs (m) | Zero-speed Tz → 5 kn Tz (s) | Acceleration RMS after/before | Heave RMS after/before |
| --- | ---: | ---: | ---: | ---: |
| JONSWAP | 0.27 | 2.970 → 4.761 | 0.430 | 0.913 |
| JONSWAP | 1.50 | 4.926 → 6.963 | 0.427 | 0.915 |
| JONSWAP | 4.00 | 6.831 → 8.777 | 0.529 | 0.935 |
| JONSWAP | 8.50 | 8.813 → 11.126 | 0.574 | 1.027 |
| PM-Stokes | 0.27 | 2.956 → 4.458 | 0.503 | 0.936 |
| PM-Stokes | 1.50 | 4.633 → 6.689 | 0.412 | 0.919 |
| PM-Stokes | 4.00 | 6.329 → 8.398 | 0.513 | 0.953 |
| PM-Stokes | 8.50 | 8.120 → 10.535 | 0.547 | 1.018 |

Tz is the displacement spectral-moment period; it differs from the peak encounter period. For example, the largest JONSWAP record has peak encounter period 13.03 s and Tz 11.13 s. Slower acceleration is not evidence of still water. Sensor noise and physical bias drift retain their original units and are not time-scaled to manufacture an invariant RMS.

## PII phase correction

The existing amplitude schedule can drive acceleration-filter time constant tau_a to its 0.15 s lower bound as acceleration falls. Meanwhile the wave gets slower, and the restoring PII core has a larger phase lead relative to that wave. The old acceleration filter then supplies insufficient compensating lag. This is a dimensionless frequency-response error, rather than a reason to weaken the accuracy gates.

For a frozen continuous-time PII core without the separate bias channel, its response from true displacement to reconstructed displacement before acceleration smoothing is C = (i omega / (r + i omega))^3. Let A = Re(C), B = Im(C), and g = 2A - |C|². The positive minimizer of |C/(1+i x)-1|², x = omega tau_a, is x* = 2B/(g + sqrt(g²+4B²)).

The library computes omega from the already filtered measured scheduling frequency and r from its clamped command. On the positive-phase branch A > 0, B > 0 it raises the empirical tau_a command to at least x*/omega. Existing 0.15–0.90 s bounds and parameter smoothing still apply. Existing stronger noise smoothing is retained. The restoring pole, bias-channel law and Mahony attitude algorithm are unchanged. No boat speed, wave truth or simulator-only correction enters the filter. `match_accel_filter_phase=false` reproduces the previous schedule.

This scalar result does not guarantee improvement for every broadband or nonlinear execution. An unguarded phase replacement worsened the largest zero-speed case; retaining the empirical smoothing floor avoids that loss. The independent C++ tests check the complex frequency response, unchanged restoring poles, bounds, and actual sinusoidal reconstruction.

## Paired shipping-code measurements

The same rebuilt executable runs both settings, with default seeds and held-out W3D_SEED values 7 and 23. Each run has 240,000 samples at 200 Hz; RMS uses all 180,000 samples in the final 900 s. There are 96 runs and 48 before/after pairs. The original source snapshot is now published at c94389a3bcedb632abd025e7db8d7dd367c94b08; the subsequent config-field placement preserves existing aggregate member order.

| Dataset and seeds | Mean normalized heave RMS before → after (% Hs) | Worst before → after (% Hs) |
| --- | ---: | ---: |
| Zero-speed, default | 5.898 → 5.880 | 8.749 → 8.749 |
| Zero-speed, all three seeds | 5.901 → 5.883 | 9.148 → 9.148 |
| 5 kn, default | 12.166 → 8.758 | 15.387 → 10.487 |
| 5 kn, all three seeds | 12.128 → 8.734 | 15.387 → 10.487 |

All 24 moving-record pairs improve. The mean reduction is 28.0%, and the worst normalized error falls by 31.8%. The zero-speed worst case is unchanged. Four zero-speed pairs have small increases; the largest absolute cost is 0.316 mm RMS on a 0.27 m PM-Stokes record. Paired roll, pitch and yaw metrics are identical. These are controlled records and sensor draws, not population confidence intervals.

| Moving record, default seed | Previous heave (% Hs) | Phase-matched heave (% Hs) |
| --- | ---: | ---: |
| jonswap_H0.270_L14.047_A30.00_P60.00 | 8.886 | 7.008 |
| jonswap_H1.500_L50.710_A-30.00_P120.00 | 15.387 | 9.551 |
| jonswap_H4.000_L112.766_A30.00_P30.00 | 12.918 | 8.514 |
| jonswap_H8.500_L202.839_A-30.00_P72.00 | 13.766 | 10.487 |
| pmstokes_H0.270_L14.047_A30.00_P60.00 | 8.320 | 6.900 |
| pmstokes_H1.500_L50.710_A-30.00_P120.00 | 14.733 | 9.499 |
| pmstokes_H4.000_L112.766_A30.00_P30.00 | 11.549 | 8.481 |
| pmstokes_H8.500_L202.839_A-30.00_P72.00 | 11.767 | 9.623 |

Production moving-data heave limits tighten from 15.47/14.81 to 10.54/9.68% Hs for JONSWAP/PM-Stokes. The yaw limit remains 8.28°. No gate is relaxed. The PII Makefile now records header dependencies, so a library-header change rebuilds its simulator and regression tests.

## Cutoffs and clamps that were tested and rejected

The full [exploratory record](../reports/results/encounter_scaling/exploratory_ablation.json) contains input physics, period-only diagnostics, executable hashes, all ablation rows, and the temporary include-overlay patches. Those overlays are experiments and are not shipping-source evidence.

| Candidate | Observed result | Decision |
| --- | --- | --- |
| Increase acceleration-variance horizon cap, 35 → 60 s | The 35 s cap is active on the largest moving waves, but changing it alone barely changes OU-III RMS. | Not adopted as an RMS fix. |
| Lower period-estimator leak, 0.02 → 0.01 Hz | Small-wave OU-III heave rises from 4.25 to 8.35% Hs at zero speed and from 8.24 to 9.53% at 5 kn. Sensor/bias contamination enters the leaky displacement proxy. | Rejected. |
| Set leak to 0.1 times estimated frequency | Mixed results, including a 5.3% increase in largest-wave moving heave error. Small-wave live period drifts toward 7.66 s against truth Tz 4.76 s. | Rejected. |
| Remove OU-III S-noise floor | Moving small-wave heave rises from 8.24 to 26.01% Hs. | Rejected. |
| Raise fixed S-noise floor to 0.225 | Moving small-wave heave improves to 7.62%, but zero-speed error increases to 4.59%. | Not a general scaling repair. |
| Reduce proxy correction gains | Little heave gain; quarter gains increase accelerometer-bias RMS to about 0.09 m/s². | Rejected. |
| Increase NLO theta/frequency gain from 0.5575 to 1.0 | Small moving-wave heave improves from 21.90 to 9.48% Hs, but largest moving-wave error worsens from 7.72 to 13.39%. | Reject a global gain change; further work needs signal/noise-aware scheduling. |

The separate 0.2 Hz acceleration tracker floor feeds the OU direction carrier. It does not supply the OU tuning frequency, so lowering it alone would not fix these heave regressions. The wave tuner has its own 0.03 Hz guard and the period estimator. No broad claim is made that every remaining RMS increase is a software bug: lower encounter frequency also lowers acceleration signal-to-noise ratio at unchanged physical sensor noise.

## Reproduction

```sh
make -C tests/pii_observer build
python3 tools/encounter_scaling_study.py \
  --stationary-archive /path/to/sim-data-files-vessel-rao-28ft.zip \
  --stationary-dir /path/to/zero-speed-v1.3.1 \
  --moving-dir /path/to/moving-v1.3.1 \
  --seeds default,7,23 --jobs 3 \
  --output reports/results/encounter_scaling/pii_phase.json
```

The driver verifies the zero-speed archive, both sets of input hashes, exact clocks, complete scoring windows, finite metrics and unchanged executable identity across the pair set. Results include the source closure, compiler environment and all paired values. Source-bound hosted evidence and the full native regression suite remain separate required checks.
