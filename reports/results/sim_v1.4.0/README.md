# v1.4.0 bidirectional replay results

Both checksum-pinned 28 ft RAO datasets use 5 kn and ±3° steering with a 20 s period.
There are 272 filter/profile/sea cells per direction, 544 total. The measured limits
are fixed separately for each cell; normal runs do not fit new limits.

| Case | Full scores | Provenance | Heave | Attitude | Bias |
| --- | --- | --- | --- | --- | --- |
| Northbound, following seas | [CSV](following/metrics.csv) | [JSON](following/measurements.json) | [RMS](following/disp_z_pct_hs.svg) | [Roll](following/roll_rms_deg.svg), [pitch](following/pitch_rms_deg.svg), [yaw](following/yaw_rms_deg.svg) | [Accelerometer](following/accel_bias_3d_rms_mps2.svg), [gyro](following/gyro_bias_3d_rms_radps.svg) |
| Southbound, head/bow-quartering seas | [CSV](head/metrics.csv) | [JSON](head/measurements.json) | [RMS](head/disp_z_pct_hs.svg) | [Roll](head/roll_rms_deg.svg), [pitch](head/pitch_rms_deg.svg), [yaw](head/yaw_rms_deg.svg) | [Accelerometer](head/accel_bias_3d_rms_mps2.svg), [gyro](head/gyro_bias_3d_rms_radps.svg) |

[Worst-record summary](../../../docs/moving-vessel-results.md) ·
[Gate protocol](../../../docs/filter-simulator-gates.md) ·
[Inputs and reproduction](../../../docs/moving-vessel-replays.md)

Charts show production defaults. CSV/JSON also retain interval-rate profiles and
heave truth/proxy frontends. Per-case studies are kept under `CASE/studies/` and
published by the `simulation-studies` workflow; each receipt identifies its inputs,
executables, commands and completion status.

The full OU-III lever-arm sweep is also recorded for
[following seas](following/studies/ou3_lever_arm_study/lever_arm_report.md) and
[head seas](head/studies/ou3_lever_arm_study/lever_arm_report.md): 224 installation
runs plus 64 derivative-cutoff runs per direction, with five SVG charts each.

Detailed heave comparisons also pass all 136 full-record checks per direction:
[following seas](following/studies/heave_only_comparison/README.md) and
[head seas](head/studies/heave_only_comparison/README.md). Each generates 102
charts; the default-frontend overviews and all full-record metrics are retained
here, with the full trace/chart bundle available from its CI job.
