# Horizontal wave observation in the TFG-CARMA filters

The standard TFG filter on PR #766 uses a period-matched horizontal velocity and displacement-integral observation after it releases accelerometer-bias learning. These replays test how the same physical relationship behaves in the independent-axis and shared-mode TFG-CARMA filters. All experimental observations ran in the live C++ filter through its normal startup and measurement path.

## Shipping comparison with main

The strict 300–1,200 s records cover eight wave records in each sea direction. All 20 CARMA global 3-D caps in the two gyro profiles have the same numerical values on main and the PR. The largest displacement gain in the parent study belongs to standard TFG. [Original shipping measurements and cap comparisons](shipping-main-comparison.json)

| CARMA filter | PR 3-D RMS | Main 3-D RMS | Change |
|---|---:|---:|---:|
| QMEKF-CARMA | 0.228162 m | 0.228164 m | -0.00074% |
| QMEKF shared-mode CARMA | 0.205282 m | 0.205276 m | +0.00283% |
| TFG-CARMA | 0.228237 m | 0.228222 m | +0.00635% |
| TFG shared-mode CARMA | 0.205371 m | 0.205369 m | +0.00094% |
| Tangent-group shared-mode CARMA | 0.205369 m | 0.205354 m | +0.00722% |

## Independent sensor-noise comparison

Each percentage below is the change in pooled RMS across 32 paired replays per filter, with reserved noise seeds 509 and 2719 on eight following and eight head-sea records. Lower is better. The last column reports the **separate 16 deterministic standard records** with their reviewed global and per-record gates. The shipping control passed all 16 for each filter. Sensor-noise sensitivity runs use the simulator’s ungated mode. [All matched RMS metrics, sample windows and pair counts](pooled-rms.csv.gz)

| Observation and filter | 3-D RMS | X RMS | Y RMS | Accelerometer-bias RMS | Yaw RMS | Standard gate passes |
|---|---:|---:|---:|---:|---:|---:|
| Direct two-axis · Independent-axis TFG-CARMA | +2.96% | +6.18% | +1.62% | −4.96% | −3.02% | 7/16 |
| Direct two-axis · Shared-mode TFG-CARMA | +4.63% | +10.59% | +1.41% | −4.96% | −3.02% | 6/16 |
| Weaker two-axis · Independent-axis TFG-CARMA | +0.04% | +0.02% | +0.10% | −1.53% | −0.99% | 11/16 |
| Weaker two-axis · Shared-mode TFG-CARMA | +0.47% | +1.01% | +0.00% | −1.53% | −0.99% | 10/16 |
| Motion gain only · Independent-axis TFG-CARMA | +7.91% | +20.77% | −1.10% | +0.29% | +0.15% | 4/16 |
| Motion gain only · Shared-mode TFG-CARMA | +7.61% | +22.89% | −2.30% | +0.23% | +0.09% | 0/16 |
| Motion and accelerometer-bias gain · Independent-axis TFG-CARMA | +7.92% | +20.76% | −1.05% | −2.16% | −0.91% | 5/16 |
| Motion and accelerometer-bias gain · Shared-mode TFG-CARMA | +7.19% | +22.44% | −2.36% | −2.30% | −1.03% | 5/16 |
| Spectrum-weighted two-axis · Independent-axis TFG-CARMA | +0.49% | +0.59% | +0.74% | −2.72% | −1.67% | 11/16 |
| Spectrum-weighted two-axis · Shared-mode TFG-CARMA | +1.25% | +2.19% | +0.73% | −2.69% | −1.64% | 9/16 |
| Along the estimated wave axis · Independent-axis TFG-CARMA | +4.09% | +5.63% | +5.30% | −2.99% | −1.63% | 6/16 |
| Along the estimated wave axis · Shared-mode TFG-CARMA | +5.78% | +8.74% | +5.87% | −3.00% | −1.61% | 6/16 |
| Across the estimated wave axis · Independent-axis TFG-CARMA | −0.52% | −0.57% | −0.82% | −1.92% | −1.16% | 13/16 |
| Across the estimated wave axis · Shared-mode TFG-CARMA | −0.40% | −0.23% | −0.79% | −1.91% | −1.15% | 12/16 |
| Spectrum-weighted across-wave · Independent-axis TFG-CARMA | −0.54% | −0.67% | −0.78% | −1.23% | −0.76% | 15/16 |
| Spectrum-weighted across-wave · Shared-mode TFG-CARMA | −0.52% | −0.58% | −0.81% | −1.22% | −0.76% | 14/16 |

| Observation and sea | Independent-axis 3-D RMS, m | Change | Lower / 16 | Shared-mode 3-D RMS, m | Change | Lower / 16 |
|---|---:|---:|---:|---:|---:|---:|
| Direct two-axis · following | 0.21514 → 0.21155 | −1.67% | 11 | 0.18284 → 0.18118 | −0.91% | 8 |
| Direct two-axis · head | 0.21992 → 0.23573 | +7.19% | 8 | 0.20900 → 0.22713 | +8.67% | 8 |
| Across the estimated wave axis · following | 0.21514 → 0.21106 | −1.89% | 14 | 0.18284 → 0.17880 | −2.21% | 13 |
| Across the estimated wave axis · head | 0.21992 → 0.22165 | +0.79% | 13 | 0.20900 → 0.21103 | +0.97% | 12 |
| Spectrum-weighted across-wave · following | 0.21514 → 0.21220 | −1.36% | 14 | 0.18284 → 0.17984 | −1.64% | 13 |
| Spectrum-weighted across-wave · head | 0.21992 → 0.22046 | +0.24% | 13 | 0.20900 → 0.20967 | +0.32% | 13 |

## Accelerometer bias in sensor axes

The late native transcripts also report bias RMS separately on X, Y and Z. These are the sensor-axis errors over 300–1,200 s for the 32 reserved pairs per filter. The raw axis values are printed to six significant digits. Both the vector norm and its higher-precision machine-readable score were checked for every replay. [Per-replay sensor-axis errors](late-bias-axes.csv.gz) · [Pooled accelerometer and gyro axes](pooled-bias-axes.csv.gz)

| Observation and filter | X bias RMS, m/s² | Y bias RMS, m/s² | Z bias RMS, m/s² |
|---|---:|---:|---:|
| Direct two-axis · Independent-axis TFG-CARMA | 0.01175 → 0.01158 (−1.52%) | 0.03428 → 0.03243 (−5.39%) | 0.00173 → 0.00173 (−0.06%) |
| Direct two-axis · Shared-mode TFG-CARMA | 0.01175 → 0.01158 (−1.44%) | 0.03432 → 0.03247 (−5.39%) | 0.00174 → 0.00174 (+0.04%) |
| Spectrum-weighted across-wave · Independent-axis TFG-CARMA | 0.01175 → 0.01168 (−0.63%) | 0.03428 → 0.03383 (−1.30%) | 0.00173 → 0.00173 (−0.01%) |
| Spectrum-weighted across-wave · Shared-mode TFG-CARMA | 0.01175 → 0.01167 (−0.66%) | 0.03432 → 0.03388 (−1.29%) | 0.00174 → 0.00174 (−0.01%) |

## When the difference begins

All comparisons use the same 200 Hz vessel and sensor inputs within each pair. The 40–60 s window contains 4,000 samples per case. The post-release window below contains 12,000 samples, and the late window contains 180,000. The discovery archive also measures 60–150 s. [Reserved records, binary hashes, activation and gate status](reserved-record-rms.csv.gz)

| Window and filter | Stock 3-D RMS | Direct observation | Across-wave observation | Spectrum-weighted across-wave |
|---|---:|---:|---:|---:|
| 40–60 s · Independent-axis TFG-CARMA | 1.23005 m | 1.23005 m | 1.23005 m | 1.23005 m |
| 40–60 s · Shared-mode TFG-CARMA | 1.26199 m | 1.26199 m | 1.26199 m | 1.26199 m |
| 60–120 s · Independent-axis TFG-CARMA | 2.81399 m | 2.81399 m | 2.81399 m | 2.81399 m |
| 60–120 s · Shared-mode TFG-CARMA | 3.09544 m | 3.09544 m | 3.09544 m | 3.09544 m |
| 300–1,200 s · Independent-axis TFG-CARMA | 0.21754 m | 0.22397 m | 0.21642 m | 0.21637 m |
| 300–1,200 s · Shared-mode TFG-CARMA | 0.19636 m | 0.20544 m | 0.19558 m | 0.19533 m |

## Per-record displacement regressions

Counts cover the same 32 reserved cases per filter. A small pooled reduction can coexist with a larger individual loss. Wave heights in the table are significant heights from the input filenames. [Per-case scores](reserved-record-rms.csv.gz)

| Observation and filter | Better by >0.5% | Worse by >0.5% | Largest relative increase | 3-D RMS | Change |
|---|---:|---:|---|---:|---:|
| Direct two-axis · Independent-axis TFG-CARMA | 18 | 13 | head, 8.500 m, seed 2719 | 0.43525 → 0.48083 m | +10.47% |
| Direct two-axis · Shared-mode TFG-CARMA | 16 | 14 | head, 8.500 m, seed 2719 | 0.41511 → 0.46519 m | +12.07% |
| Across the estimated wave axis · Independent-axis TFG-CARMA | 18 | 2 | head, 8.500 m, seed 509 | 0.43486 → 0.44455 m | +2.23% |
| Across the estimated wave axis · Shared-mode TFG-CARMA | 19 | 2 | head, 8.500 m, seed 509 | 0.41672 → 0.42690 m | +2.44% |
| Spectrum-weighted across-wave · Independent-axis TFG-CARMA | 16 | 2 | head, 8.500 m, seed 509 | 0.43486 → 0.43939 m | +1.04% |
| Spectrum-weighted across-wave · Shared-mode TFG-CARMA | 16 | 2 | head, 8.500 m, seed 509 | 0.41672 → 0.42149 m | +1.15% |

## Measurement and controls

For the deployed zero-at-DC CARMA(4,1) shape, the horizontal displacement integral at fixed model tuning can be expressed using three existing modal states. With position p, velocity v, fourth state l, first-mode angular frequency ω₁ and damping ζ₁, S = (l − 2ζ₁ω₁p − v)/ω₁². The CARMA transition gives dS/dt = p between retunes. Both CARMA implementations share that four-state transition. The adaptive trial uses the current model parameters and its estimated encounter frequency for β = 1/(2πf)², and observes the horizontal residual S + βv. The CARMA process already couples position and velocity, so a second observation changes how the filter apportions error among wave motion, attitude and sensor biases.

The direct two-axis arm uses 0.625 times the existing horizontal pseudo-observation covariance and the weaker arm uses 2.5 times that covariance. The motion-only trials constrain which Kalman gain rows can update attitude and sensor biases. The spectrum-weighted arm adds a model-based residual noise floor. The along-wave and across-wave arms project onto the estimated propagation axis or its perpendicular. The combined arm projects across the estimated wave axis and adds the residual noise floor. These projections use the onboard wave-direction estimate in the heading-aligned boat frame. Records where that estimate closes its confidence gate remain in the comparison. The source patches, variant settings and executable controls are in the original study archive. [Discovery scores, patches, manifests and replay driver](study.tar.gz)

The discovery controls reproduced 576 hosted score values from the branch strict matrix. Disabled-observation controls reproduced 7,579 of 7,776 recorded stock RMS values in their matched windows. All 197 differences occurred in the executable built for the restricted-gain trials. The largest relative difference was 0.0296%. Each reserved trial is compared with the disabled mode of the same binary. New tests used the same input hashes as the shipping replays, the pinned GCC 14.2 toolchain and the reference Ubuntu libm. The complete experiment contains 1,856 actual full replays, including 1,024 discovery runs and 832 additional runs on the two reserved noise seeds. [Source, input hashes, executable hashes and per-sea comparisons](study-summary.json)

To reproduce an original discovery record, unpack the study archive, apply the named patch in a clean checkout of 85e7284972e4dc046bdc34cbcb4f044fbbc2a4ab, build the corresponding native simulator, and invoke the included reproduce.py with its manifest. For the reserved comparison use the same binary and patch, the arm’s CARMA_TRANSFER value, W3D_SEED=509 or 2719 and W3D_REFERENCE_MOTION=underway. Set W3D_VALIDATION_SEGMENTS to early:40:60,release:60:120,late:300:1200 and W3D_VALIDATION_WINDOW_SEC to 900. Pass the pinned CSV with --input. The reserved CSV includes each input and binary digest, the exact original log path and SHA-256. The six compressed log files cover both implementations and all sensor seeds.

