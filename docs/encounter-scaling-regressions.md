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

## NLO low-acceleration bandwidth protection

The encounter-frequency schedule can reduce observer bandwidth while displacement remains substantial. The library adapter now blends in the configured base bandwidth when measured wave-band acceleration is weak. With band variance V, empirical scale s = 0.12 m/s² and weight w = 1/(1+(V/s²)²), the target is sqrt(theta_wave² + theta_base² w), inside the existing bandwidth limits and smoothing. Setting `theta_accel_scale_mps2=0` selects the previous schedule. The energy measurement uses the prior attitude to level measured specific force and the shared 0.5–4 times wave-frequency band; it does not use the observer's translational correction. The attitude and frequency estimates remain coupled to the observer. This is an empirical bandwidth safeguard, not a sensor-noise estimate or a still-water classification.

The three-seed [normal-noise comparison](../reports/results/encounter_scaling/nlo_band.json) contains 96 runs / 48 pairs. Separate [900 rpm engine](../reports/results/encounter_scaling/nlo_band_engine.json) and [noise-free](../reports/results/encounter_scaling/nlo_band_noise_free.json) controls each contain 32 runs / 16 pairs.

| Dataset and seeds | Mean heave before → after (% Hs) | Worst before → after (% Hs) |
| --- | ---: | ---: |
| Zero-speed, default | 6.953 → 6.845 | 7.322 → 6.984 |
| Zero-speed, all three seeds | 7.194 → 7.049 | 8.509 → 7.565 |
| 5 kn, default | 11.518 → 8.111 | 23.524 → 9.430 |
| 5 kn, all three seeds | 16.663 → 8.474 | 68.161 → 10.159 |

Sixteen of 24 moving pairs improve. The largest default-seed moving cost is on the largest JONSWAP record: 7.7168 → 7.7423% Hs (2.16 mm RMS). The all-seed roll/pitch maxima change by less than 0.000003°. The NLO moving heave limits tighten from 22.02/23.65 to 9.48/9.17% Hs for JONSWAP/low PM-Stokes; the existing 8.19% limit on the other PM-Stokes records remains.

Engine-noise moving mean/worst improve from 11.832/24.072 to 8.098/9.452% Hs. Zero-speed small JONSWAP at 900 rpm costs 0.407 mm RMS (6.961 → 7.112% Hs). A first broadband-energy implementation cost 2.84 percentage points in that case and was rejected; wave-band energy is the subsequent refinement. Noise-free moving mean increases from 5.989 to 6.995% Hs, including 4.106 → 7.775% on small JONSWAP (9.9 mm). The bandwidth safeguard therefore trades some ideal-sensor waveform accuracy for robustness to sensor bias and noise. It is not an across-the-board improvement.

The device-path synthetic test independently exercises sinusoidal motion plus slow additive acceleration drift: the short, small wave improves from 0.574 to 0.0448 m RMS and the larger, slower wave from 0.715 to 0.702 m, without detrending, fitting, or output suppression. Existing delayed/jittered device, local-gravity and magnetic-turn tests remain required.

## Shared-mode CARMA vertical band normalization

The amplitude observer measures acceleration after a band filter, while the
vertical CARMA prior expects total model acceleration RMS. The shared QMEKF,
TFG and tangent wrappers now account for that difference. The frozen default
CARMA(4,1) model passes 0.58598012624 of its acceleration variance through the
0.5–4 frequency band, giving a standard-deviation factor of 1.30634701277.
Two small covariance solves compute the fraction when shape or band ratios
change. Horizontal empirical factors retain their previous convention. This
uses measured scheduling parameters and the configured model, without vessel
speed or wave truth. It approximates a changing sampled band with a frozen
continuous model; it does not guarantee exact normalization of a real sea.

The [QMEKF](../reports/results/encounter_scaling/shared_band_csm.json),
[TFG](../reports/results/encounter_scaling/shared_band_tcsm.json) and
[tangent](../reports/results/encounter_scaling/shared_band_tangent.json)
qualifications each contain 96 full replays: eight seas, both speeds, default
seed plus seeds 7 and 23, and enabled/disabled arms of the same native executable.

| Shared-mode family | Moving mean, three seeds (% Hs) | Moving worst, three seeds (% Hs) | Zero-speed mean, three seeds (% Hs) |
| --- | ---: | ---: | ---: |
| QMEKF | 2.71437 → 2.69209 | 3.70165 → 3.58339 | 2.32357 → 2.28199 |
| TFG | 2.71450 → 2.69235 | 3.70215 → 3.58331 | 2.32342 → 2.28185 |
| Tangent | 2.71456 → 2.69242 | 3.70216 → 3.58334 | 2.32347 → 2.28192 |

For QMEKF, default-seed moving mean/worst change from 2.68488/3.70165 to
2.63526/3.53167% Hs. Across all seeds, 16 of 24 moving pairs improve;
the largest moving cost is 3.325 mm RMS and the largest zero-speed cost is
6.404 mm, both on the largest JONSWAP sea with seed 7. Mean moving 3-D
displacement improves from 0.154984 to 0.151292 m. Mean pitch slightly
increases from 0.121732 to 0.123154 degrees and yaw from 2.422058 to 2.426963
degrees. This is a modest improvement, not a universal one. Because the
zero-speed mean also improves, the new zero-to-moving gap need not shrink.

All eight default moving records pass their unchanged production gates for
each shared family. Held-out seeds are scored diagnostics; their harness does
not apply default-seed acceptance bars. Zero-speed runs retain any moving-bar
misses and are not presented as a separate gate qualification. The independent
axis wrappers retain their previous default: the earlier held-out amplitude
screen did not support enabling the correction there.

The separate 900 rpm controls also improve: QMEKF moving mean/worst
3.21590/4.07942 -> 3.12973/3.77766% Hs, and TFG
3.06692/4.20859 -> 2.93138/3.89966% Hs. Both zero-speed means improve.
The [QMEKF](../reports/results/encounter_scaling/shared_band_csm_engine.json)
and [TFG](../reports/results/encounter_scaling/shared_band_tcsm_engine.json)
controls each contain 32 replays. Normal-noise gates fail in both engine arms;
these are degradation comparisons, not claims that the engine cases satisfy
ordinary acceptance limits. No sea-transition qualification of this new
normalization is claimed.

The [additional localization record](../reports/results/encounter_scaling/residual_diagnostics.json)
contains the rejected OU candidates and oracle diagnostics, with source/input/
binary identities and overlay patches. Supplying the true period changes OU-III
screened moving mean from 6.04498 to 6.06986% Hs; removing cruise-turn forcing
changes small-wave heave only from 8.23605 to 8.23030%. Period error and that
forcing do not explain most of the residual. Removing accelerometer error has
a much larger effect than removing gyro or magnetometer error. This localizes
further work to motion-prior and acceleration/bias separation, without proving
an irreducible sensor floor. Fitting shape with wave truth worsens small-wave
CARMA heave from 3.550 to 5.084% Hs and is rejected. Vertical-only OU-III
stillness handling improves default moving mean but worsens the zero-speed
mean for every tested seed and fails yaw, 0.880571 > 0.875540 degrees. The one
refinement preserving acceleration weighting fails the same gate. Neither
OU-III candidate is deployed, and no acceptance limit is widened.

## Remaining OU, CARMA and TFG speed regressions

The [complete inventory](../reports/results/encounter_scaling/remaining_regressions.json), recorded before the shared-band normalization above, pairs all eight seas for every gated OU/CARMA/TFG replay family, using the same executable at both speeds. It includes displacement, roll, pitch, yaw and both bias-error RMS metrics. These values compare physical datasets after the earlier stillness/PII/NLO fixes; they are not before/after implementation improvements. The default sensor seed, production sampling profile and final 900 s are used. TFG legacy here means its regularization law, not the interval-rate sampling profile. Ungated CARMA spectral-shape research variants are not covered by these acceptance measurements.

| Family | Mean heave, zero-speed → 5 kn (% Hs) | Worst heave, zero-speed → 5 kn (% Hs) | Seas with higher heave error |
| --- | ---: | ---: | ---: |
| OU-II | 6.216 → 7.348 | 6.642 → 8.486 | 8/8 |
| OU-III | 3.996 → 5.601 | 4.249 → 8.236 | 8/8 |
| TFG spectral law | 4.058 → 5.520 | 4.383 → 7.819 | 8/8 |
| TFG legacy law | 4.196 → 5.926 | 4.322 → 8.247 | 8/8 |
| CARMA independent axes | 2.311 → 2.745 | 2.461 → 3.726 | 7/8 |
| CARMA shared modes | 2.273 → 2.685 | 2.455 → 3.702 | 6/8 |
| TFG-CARMA independent axes | 2.311 → 2.745 | 2.461 → 3.727 | 7/8 |
| TFG-CARMA shared modes | 2.273 → 2.685 | 2.455 → 3.702 | 6/8 |

The OU-II and independent-axis CARMA/TFG-CARMA stillness fixes reduce the worst moving errors from about 20.84 to 8.486% Hs and 4.67 to 3.73% Hs respectively. They do not remove the residual speed penalty above. The shared-mode correction also leaves a residual penalty; OU-III remains unresolved.

The [continuation experiments](../reports/results/encounter_scaling/continuation_trials.json) preserve candidate patches, compilation details, input/source/binary hashes and scored values. Their outcomes constrain the next changes:

* **Covariance timing:** a [48-run production screen](../reports/results/encounter_scaling/covariance_cadence.json) tests small, medium and large JONSWAP seas at both speeds. Disabling periodic OU-III covariance synchronization lowers worst screened heave 8.23605 → 8.16948% Hs but raises worst accelerometer-bias RMS 0.05680 → 0.09373 m/s². A 0.5 s tuning-commit interval also raises that bias metric to 0.08371 m/s². Neither is adopted. The TFG arm used the separately rejected cadence candidate and is labeled accordingly.
* **TFG fractional cadence:** preserving the fractional timer remainder slightly helps the spectral law, but the legacy-law replay fails three unchanged gates: heave 8.3259 > 8.28795% Hs, Z-bias error 5.52383 > 5.51543% of maximum true bias, and heave 8.10862 > 8.05605% Hs. The candidate is rejected and the original source restored. Both [spectral](../reports/results/encounter_scaling/tfg_cadence.json) and [legacy](../reports/results/encounter_scaling/tfg_legacy_cadence.json) paired results remain diagnostic evidence, separate from the accepted inventory.
* **TFG held-bias covariance:** a Schmidt-style startup overlay worsens small stationary pitch 0.0458 → 0.5440° and accelerometer-bias RMS 0.0506 → 0.11685 m/s²; large moving yaw rises 0.674 → 1.789°. The 12-run screen rejects it.
* **Vertical-only CARMA band compensation:** the trial uses a provisional acceleration-standard-deviation multiplier 1.26686. Its original exact continuous-band interpretation was incorrect: independent infinite-domain spectral quadrature and an augmented Lyapunov solve give variance fraction 0.58598012624 and multiplier 1.30634701277 for the stated frozen default CARMA(4,1) model and 0.5–4 frequency band. The recorded runs do not test that corrected factor. Applying the provisional multiplier vertically, with the empirically calibrated horizontal amplitudes retained, gives small default-seed gains across all four families and all eight-sea means. The 256-run screen also includes held-out moving seeds 7 and 23. For independent-axis CARMA those seeds worsen mean heave 2.79308 → 2.80640% Hs, with a maximum 6.46 mm cost. The tested multiplier does not reliably repair the regression and is not adopted; this does not rule out the corrected mathematical factor or a different band model. No held-out zero-speed or transition qualification is claimed.
* **Existing CARMA evidence mode:** in an eight-run small/large JONSWAP screen of both independent-axis cores at both speeds, replacing legacy excitation with the existing noise-aware mode changes small moving heave by only about 0.008 percentage point and slightly worsens 3-D RMS. Large-wave changes are negligible. Shared-mode filters already use this mode. It is not a material repair of these residual errors.

No acceptance bar is weakened. These finite diagnostics do not establish that a remaining regression is unavoidable: they reject the particular timing, covariance, startup and amplitude changes tested. The next useful mechanism must separate wave acceleration from attitude/bias error at the lower encounter acceleration, and preserve that separation across seeds and zero-speed controls.

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
| Increase NLO theta/frequency gain from 0.5575 to 1.0 | Small moving-wave heave improves from 21.90 to 9.48% Hs, but largest moving-wave error worsens from 7.72 to 13.39%. | Reject the global gain change; the wave-band safeguard above is evaluated separately. |

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

Use `--family nlo` for NLO, adding `--engine-rpm 900` or `--noise none` for its controls. Use `--family csm`, `tcsm`, or `tangent` for the shared-mode band-normalization ablation. `--resume` recovers a partial run only when sources, executable, inputs, seeds and settings still match.

The driver verifies the zero-speed archive, both sets of input hashes, exact clocks, complete scoring windows, finite metrics and unchanged executable identity across the pair set. Results include the source closure, compiler environment and all paired values. Source-bound hosted evidence and the full native regression suite remain separate required checks.
