# Startup learning: measured results

Base: `8c9f535709fe7ff4d87d75dc5c98f4d0a48daa7b`. Seeds: 101; power-on offsets: 17,107 s. Both pinned v1.4.0 datasets.

These are held-out sensor seeds and startup phases on the existing wave records, not independent newly generated seas. Values below pool one endpoint error per trial; they are not live-filter RMS over a navigation run.

## Equal-epoch handover comparison

Every candidate row in this table uses exactly the shipping handover timestamp. A shipping startup missing the 35 s deadline is counted separately below.

| Direction | Sensor profile | Method | N | Yaw RMS ° | Tilt RMS ° | Field direction RMS ° | Field magnitude RMS µT | Gyro bias RMS °/s |
|---|---|---|---:|---:|---:|---:|---:|---:|
| following | ideal | shipping | 16 | 1.502 | 0.671 | 0.033 | 0.011 | 0.0000 |
| following | ideal | uniform | 16 | 0.455 | 0.212 | 0.067 | 0.000 | 0.0000 |
| following | ideal | tapered | 16 | 0.071 | 0.037 | 0.019 | 0.000 | 0.0000 |
| following | ideal | magfit | 16 | 0.071 | 0.037 | 0.019 | 0.000 | 0.0000 |
| following | ideal | joint | 16 | 0.096 | 0.048 | 0.019 | 0.000 | 0.0048 |
| following | ideal | joint_cal | 16 | 0.096 | 0.048 | 0.020 | 0.002 | 0.0049 |
| following | interference | shipping | 16 | 3.038 | 0.728 | 4.029 | 2.044 | 0.0532 |
| following | interference | uniform | 16 | 2.497 | 1.113 | 4.620 | 2.312 | 0.0532 |
| following | interference | tapered | 16 | 2.746 | 1.227 | 4.606 | 2.312 | 0.0532 |
| following | interference | magfit | 16 | 3.107 | 7.425 | 4.459 | 2.446 | 0.4465 |
| following | interference | joint | 16 | 3.676 | 7.299 | 4.465 | 2.445 | 0.4366 |
| following | interference | joint_cal | 16 | 48.681 | 7.510 | 7.473 | 5.165 | 0.5217 |
| head | ideal | shipping | 16 | 2.207 | 0.837 | 0.244 | 0.019 | 0.0000 |
| head | ideal | uniform | 16 | 0.141 | 0.066 | 0.023 | 0.000 | 0.0000 |
| head | ideal | tapered | 16 | 0.020 | 0.019 | 0.017 | 0.000 | 0.0000 |
| head | ideal | magfit | 16 | 0.020 | 0.019 | 0.017 | 0.000 | 0.0000 |
| head | ideal | joint | 16 | 0.037 | 0.023 | 0.017 | 0.000 | 0.0017 |
| head | ideal | joint_cal | 16 | 0.038 | 0.023 | 0.017 | 0.001 | 0.0017 |
| head | interference | shipping | 16 | 4.497 | 1.006 | 4.342 | 1.661 | 0.0533 |
| head | interference | uniform | 16 | 4.776 | 1.196 | 4.921 | 1.834 | 0.0533 |
| head | interference | tapered | 16 | 4.774 | 1.196 | 4.920 | 1.834 | 0.0533 |
| head | interference | magfit | 16 | 7.345 | 8.632 | 5.047 | 1.666 | 0.5833 |
| head | interference | joint | 16 | 7.893 | 8.447 | 5.034 | 1.666 | 0.5647 |
| head | interference | joint_cal | 16 | 60.316 | 7.125 | 8.628 | 5.049 | 0.5992 |

## Fixed 35-second endpoint

| Direction | Sensor profile | Method | N | Yaw RMS ° | Tilt RMS ° | Gyro bias RMS °/s | Ready % | False ready % of ready |
|---|---|---|---:|---:|---:|---:|---:|---:|
| following | ideal | joint | 16 | 0.116 | 0.050 | 0.0076 | 100.0 | 0.0 |
| following | ideal | joint_cal | 16 | 0.117 | 0.050 | 0.0076 | 100.0 | 0.0 |
| following | interference | joint | 16 | 4.138 | 7.104 | 0.3882 | 0.0 | — |
| following | interference | joint_cal | 16 | 46.128 | 7.117 | 0.4798 | 0.0 | — |
| head | ideal | joint | 16 | 0.060 | 0.026 | 0.0029 | 100.0 | 0.0 |
| head | ideal | joint_cal | 16 | 0.060 | 0.026 | 0.0030 | 100.0 | 0.0 |
| head | interference | joint | 16 | 8.596 | 8.216 | 0.4909 | 0.0 | — |
| head | interference | joint_cal | 16 | 56.288 | 7.102 | 0.5221 | 0.0 | — |

The diagnostic readiness target is yaw ≤2° and tilt ≤1°. Readiness uses sensor statistics; target failures are evaluated using truth afterward. Neither the flag nor a small residual certifies accuracy.

## Deadline and latency

| Direction | Sensor profile | Method | Trials | Shipping misses at 35 s | Shipping median s | Accurate by 35 s | Median first ready s | False first readiness |
|---|---|---|---:|---:|---:|---:|---:|---:|
| following | ideal | joint | 16 | 0 | 31.84 | 16/16 | 20.0 | 0/16 |
| following | ideal | joint_cal | 16 | 0 | 31.84 | 16/16 | 20.0 | 0/16 |
| following | interference | joint | 16 | 0 | 31.84 | 0/16 | — | 0/0 |
| following | interference | joint_cal | 16 | 0 | 31.84 | 0/16 | — | 0/0 |
| head | ideal | joint | 16 | 0 | 31.84 | 16/16 | 20.0 | 0/16 |
| head | ideal | joint_cal | 16 | 0 | 31.84 | 16/16 | 20.0 | 0/16 |
| head | interference | joint | 16 | 0 | 31.84 | 0/16 | — | 0/0 |
| head | interference | joint_cal | 16 | 0 | 31.84 | 0/16 | — | 0/0 |

Host joint solve at 35 s: median 9.7 ms; maximum 13.5 ms in this run. This does not qualify ESP32 execution time or memory.

Accelerometer bias is retained as a calibration prior, so its error is unchanged. The gyroscope/field covariance is conditional and cannot be substituted for a full correlated live-filter handover covariance.

The reproducible raw metrics and actual handover components are in `metrics.csv.gz`; `coverage.json` retains every deadline miss and first-readiness result. `provenance.json` records input and source hashes.
