# Heading regression correction

The shared tilt-only hard-iron fit could absorb residual magnetic scale, cross-axis and alignment errors into an additive offset. Changing mean heading exposed that error: the southbound low-wave OU-III replay had yaw RMS 2.542 deg and mean error -2.530 deg; disabling injected sensor errors reduced its RMS to 0.100 deg. The CSV gyro/quaternion consistency check passes.

The live common facade now uses causal gyro-relative magnetic calibration of field, offset, nuisance gyro bias and residual affine matrix. No truth attitude, speed, field direction or main-filter state enters the learner. Quiet/no-excitation data does not apply a calibration. Offset and affine correction slew together, with the reference recomputed from the same raw statistics. See `../../../docs/continuous-mag-hard-iron.md` for model and safeguards.

PII also had inconsistent body-axis conversions for magnetometer versus gyro/accelerometer and a fixed two-declination yaw workaround. All sensors and the output quaternion now use one NED/NWU basis conversion. The standalone AHRS replay used gyro-bias variance density 1000 times the firmware value; both now share the unchanged firmware default, 2e-7 rad^2/s^3.

## Complete replay matrix

Worst yaw RMS in degrees across each family's applicable default and interval profiles, all eight released records. The paired before/after columns use identical v1.4.0 inputs and default sensor draws. The release column is descriptive: v2.7.7 used v1.3.1 following seas with 120 s steering, whereas v1.4.0 adds head seas and 20 s steering.

| Family | v2.7.7 | Before fix, both cases | Fixed following | Fixed head |
| --- | ---: | ---: | ---: | ---: |
| ahrs | 7.080 | 8.810 | 1.426 | 1.654 |
| carma | 0.938 | 2.511 | 0.904 | 0.999 |
| csm | 0.977 | 2.501 | 0.909 | 0.934 |
| ou2 | 0.899 | 2.579 | 0.909 | 1.092 |
| ou3 | 0.871 | 2.542 | 0.807 | 1.093 |
| pii | 8.233 | 14.973 | 2.422 | 3.011 |
| tangent | not present | 2.494 | 0.894 | 0.930 |
| tcarma | 0.927 | 2.503 | 0.894 | 0.994 |
| tcsm | 0.947 | 2.494 | 0.896 | 0.930 |
| tfg | 1.069 | 2.607 | 0.932 | 1.100 |
| tfg_legacy | 3.063 | 3.848 | 3.738 | 1.262 |

Main wave-filter worst yaw is now 0.930–1.100 deg across both cases. The legacy cubic TFG comparator still reaches 3.739 deg in following seas (previously 3.848 deg); this is not claimed as resolved to the accuracy of the production SpectralMSE TFG.

## Tradeoffs and held-out draws

Individual high-sea yaw cells can worsen despite the lower family maxima. Relative to the steering-only implementation, the largest extra Y RMS is 8.63 mm, and PII head-sea heave rises from 9.570% to 10.327% of Hs in the 4 m JONSWAP record. The largest pitch/accelerometer-bias increases occur in the legacy cubic TFG comparator. Several head-sea gyro-bias errors also rise; the largest absolute increase is 3.995e-5 rad/s. These changes are preserved in the paired JSON and the new strict per-record caps.

Across 24 OU-III held-out sensor-initialization replays (four JONSWAP sea states, two headings, seeds 7/41/101), pooled yaw RMS falls 2.161 to 0.797 deg and worst RMS falls 4.165 to 1.513 deg; 19/24 improve. Pooled Y RMS changes 0.23055 to 0.23152 m, heave 0.180380 to 0.180367 m, accelerometer bias 0.05614 to 0.05590 m/s^2, and gyro bias 0.0001371 to 0.0001471 rad/s. These diagnostics supplement the controlling final-source CI matrix; they do not replace it.

## Gates and reproduction

`CurrentReplayLimits.h`, `RecordReplayLimits.inc` and their JSON sources are regenerated from all 544 CI measurement cells using the unchanged 0.5% margin and upward 1e-5 rounding. All 272 horizontal Y gates remain independent of the 3-D gate. The qualifying receipts retain their original failures against the previous caps. The fresh enforced run [43994230010](https://github.com/bareboat-necessities/ocean-imu/actions/runs/43994230010), at `9e73439e`, passes all 544 cells in both directions and reproduces the original numerical records exactly.

Run `python3 tools/sim_matrix.py --build --jobs 2` to reproduce both cases. Run `make all` for native regression coverage. The six-seed/headings synthetic calibration tests and Mahony basis/tracking regression are in the normal native suites. The full matrix source, artifacts and every paired metric are bound in `yaw-regression.json`; `yaw-heldout.json` includes input, binary and log hashes for the additional paired draws.
