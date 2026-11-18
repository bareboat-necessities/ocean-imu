# Windowed TFG-CARMA replay results

The transverse-period implementation at `4b6d7fc` and its shipping predecessor at `00f1829` replayed 16 irregular-sea records using matched input, initialization, sensor noise and gyro timing. For each filter, 128 default-profile comparisons span eight draws. A further 16 comparisons use the legacy gyro timing on the fixed draw. Each record lasts 1,200 seconds. Every pooled value below is the square root of the mean of squared per-record RMS errors.

## Early and late displacement

| Filter | Draws | Early 3-D, 40–60 s (m) | Release 3-D, 60–150 s (m) | Late 3-D, 300–1,200 s (m) | Lower late 3-D | Largest increase |
|---|---|---:|---:|---:|---:|---:|
| Independent-axis | First six | 1.308632 → 1.308632 | 2.606991 → 2.606939 | 0.222231 → 0.218838 | 86/96 | +0.0688% |
| Independent-axis | 1063 and 3559 | 1.381913 → 1.381913 | 2.597780 → 2.597733 | 0.228851 → 0.225506 | 30/32 | +0.0433% |
| Shared-mode | First six | 1.335086 → 1.335086 | 2.851394 → 2.851347 | 0.201866 → 0.198516 | 87/96 | +0.0385% |
| Shared-mode | 1063 and 3559 | 1.396410 → 1.396410 | 2.803272 → 2.803228 | 0.207339 → 0.204104 | 30/32 | +0.0222% |

## Late displacement, bias and attitude by sea

Percentages below compare RMS on all eight draws. Each row covers 64 paired records. Lower is better.

| Filter | Sea | X | Y | Z | 3-D | Accel bias 3-D | Gyro bias 3-D | Roll | Pitch | Yaw |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Independent-axis | Following | -1.414% | -2.756% | -0.000% | -1.495% | -0.846% | -0.049% | -0.644% | -1.127% | -0.159% |
| Independent-axis | Head | -2.298% | -1.820% | +0.001% | -1.524% | -1.093% | +0.004% | -0.642% | -1.810% | -0.219% |
| Shared-mode | Following | -1.830% | -3.240% | -0.326% | -1.887% | -0.846% | -0.048% | -0.602% | -1.175% | -0.157% |
| Shared-mode | Head | -2.003% | -1.879% | +0.076% | -1.430% | -1.083% | +0.009% | -0.625% | -1.813% | -0.213% |

Across the eight draws, late 3-D improved on every head-sea record for both filters. On following seas, 12 of 64 independent-axis records and 11 of 64 shared-mode records had a higher late 3-D RMS. The largest increases were 0.069% and 0.039%, respectively. On the two post-freeze draws, following-sea yaw increased 0.112% and 0.085%. The largest individual yaw increase across all eight draws was 2.86% for independent-axis and 2.96% for shared-mode. Each filter also had two head-sea records where 3-D accelerometer-bias RMS rose by more than 0.5%.

All 64 fixed-draw combinations of filter, gyro profile, wave record and sea passed their pre-existing gates on both builds. The per-draw table reports every interval and count for nine primary channels, including changes greater than 0.5% in either direction. These comparisons use simulated inputs.

[Windowed pairs](windowed-pairs.csv.gz) · [Per-sea and per-draw loss counts](pooled-and-loss-counts.csv.gz) · [Original native logs](original-windowed-logs.tar.gz) · [Input, binary and log checksums](raw-input-index.json.gz) · [Verified inventory](manifest.json)

The earlier direct two-axis and spectrum-weighted variants used the vertical wave period. Their control replays and independent seed results are archived separately. [Earlier transfer experiments](../../tfg_carma_transfer/independent-seeds.md)
