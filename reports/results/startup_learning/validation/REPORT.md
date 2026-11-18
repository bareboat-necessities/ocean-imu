# Startup learning: measured results

Base: `28b69057fcd49d741ab402dd064d53ff79ec583c`. Seeds: 11,12,13,14,15; power-on offsets: 67,307,607 s. Both pinned v1.4.0 datasets.

These are held-out sensor seeds and startup phases on the existing wave records, not independent newly generated seas. Values below pool one endpoint error per trial; they are not live-filter RMS over a navigation run.

## Equal-epoch handover comparison

Every candidate row in this table uses exactly the shipping handover timestamp. A shipping startup missing the 35 s deadline is counted separately below.

| Direction | Sensor profile | Method | N | Yaw RMS ° | Tilt RMS ° | Field direction RMS ° | Field magnitude RMS µT | Gyro bias RMS °/s |
|---|---|---|---:|---:|---:|---:|---:|---:|
| following | calibrated | shipping | 109 | 2.202 | 0.905 | 0.307 | 0.045 | 0.0537 |
| following | calibrated | uniform | 109 | 0.727 | 0.888 | 0.130 | 0.026 | 0.0537 |
| following | calibrated | tapered | 109 | 0.678 | 0.877 | 0.140 | 0.026 | 0.0537 |
| following | calibrated | magfit | 109 | 0.834 | 0.346 | 0.139 | 0.025 | 0.0298 |
| following | calibrated | joint | 109 | 0.802 | 0.333 | 0.139 | 0.025 | 0.0244 |
| following | calibrated | joint_cal | 109 | 1.214 | 0.333 | 0.519 | 0.314 | 0.0246 |
| following | residual | shipping | 109 | 2.494 | 0.905 | 1.102 | 0.490 | 0.0537 |
| following | residual | uniform | 109 | 2.589 | 0.888 | 1.092 | 0.509 | 0.0537 |
| following | residual | tapered | 109 | 2.601 | 0.877 | 1.094 | 0.509 | 0.0537 |
| following | residual | magfit | 109 | 2.398 | 0.391 | 1.094 | 0.509 | 0.0318 |
| following | residual | joint | 109 | 2.426 | 0.368 | 1.094 | 0.509 | 0.0262 |
| following | residual | joint_cal | 109 | 1.876 | 0.370 | 0.531 | 0.567 | 0.0247 |
| head | calibrated | shipping | 120 | 3.715 | 1.157 | 0.661 | 0.039 | 0.0535 |
| head | calibrated | uniform | 120 | 1.029 | 0.896 | 0.173 | 0.032 | 0.0535 |
| head | calibrated | tapered | 120 | 1.014 | 0.870 | 0.140 | 0.032 | 0.0535 |
| head | calibrated | magfit | 120 | 1.045 | 0.445 | 0.141 | 0.031 | 0.0333 |
| head | calibrated | joint | 120 | 0.966 | 0.408 | 0.141 | 0.031 | 0.0251 |
| head | calibrated | joint_cal | 120 | 1.161 | 0.409 | 0.624 | 0.334 | 0.0253 |
| head | residual | shipping | 120 | 4.680 | 1.157 | 1.034 | 1.088 | 0.0535 |
| head | residual | uniform | 120 | 1.774 | 0.896 | 0.819 | 1.089 | 0.0535 |
| head | residual | tapered | 120 | 1.756 | 0.870 | 0.821 | 1.089 | 0.0535 |
| head | residual | magfit | 120 | 1.628 | 0.467 | 0.821 | 1.089 | 0.0328 |
| head | residual | joint | 120 | 1.674 | 0.427 | 0.821 | 1.089 | 0.0253 |
| head | residual | joint_cal | 120 | 1.640 | 0.437 | 0.726 | 0.422 | 0.0219 |

## Fixed 35-second endpoint

| Direction | Sensor profile | Method | N | Yaw RMS / worst ° | Tilt RMS / worst ° | Gyro bias RMS °/s | Ready % | False ready % of ready |
|---|---|---|---:|---:|---:|---:|---:|---:|
| following | calibrated | joint | 120 | 0.797 / 1.387 | 0.337 / 0.574 | 0.0197 | 100.0 | 0.0 |
| following | calibrated | joint_cal | 120 | 1.158 / 2.389 | 0.337 / 0.570 | 0.0198 | 100.0 | 5.0 |
| following | residual | joint | 120 | 2.397 / 3.782 | 0.376 / 0.743 | 0.0213 | 100.0 | 60.0 |
| following | residual | joint_cal | 120 | 1.857 / 4.194 | 0.375 / 0.743 | 0.0205 | 100.0 | 35.0 |
| head | calibrated | joint | 120 | 0.926 / 1.545 | 0.398 / 0.610 | 0.0228 | 100.0 | 0.0 |
| head | calibrated | joint_cal | 120 | 1.125 / 3.233 | 0.399 / 0.597 | 0.0230 | 100.0 | 8.3 |
| head | residual | joint | 120 | 1.641 / 2.830 | 0.408 / 0.610 | 0.0221 | 100.0 | 22.5 |
| head | residual | joint_cal | 120 | 1.644 / 3.867 | 0.416 / 0.605 | 0.0186 | 100.0 | 23.3 |

The diagnostic readiness target is yaw ≤2° and tilt ≤1°. Readiness uses sensor statistics; target failures are evaluated using truth afterward. Neither the flag nor a small residual certifies accuracy.

## Deadline and latency

| Direction | Sensor profile | Method | Trials | Shipping misses at 35 s | Shipping median s | Accurate by 35 s | Median first ready s | False first readiness |
|---|---|---|---:|---:|---:|---:|---:|---:|
| following | calibrated | joint | 120 | 11 | 31.84 | 120/120 | 20.0 | 2/120 |
| following | calibrated | joint_cal | 120 | 11 | 31.84 | 114/120 | 20.0 | 18/120 |
| following | residual | joint | 120 | 11 | 31.84 | 48/120 | 20.0 | 72/120 |
| following | residual | joint_cal | 120 | 11 | 31.84 | 78/120 | 20.0 | 44/120 |
| head | calibrated | joint | 120 | 0 | 31.84 | 120/120 | 20.0 | 0/120 |
| head | calibrated | joint_cal | 120 | 0 | 31.84 | 110/120 | 20.0 | 10/120 |
| head | residual | joint | 120 | 0 | 31.84 | 93/120 | 20.0 | 29/120 |
| head | residual | joint_cal | 120 | 0 | 31.84 | 92/120 | 20.0 | 33/120 |

Host joint solve at 35 s: median 11.5 ms; maximum 38.7 ms in this run. This does not qualify ESP32 execution time or memory.

Accelerometer bias is retained as a calibration prior, so its error is unchanged. The gyroscope/field covariance is conditional and cannot be substituted for a full correlated live-filter handover covariance.

The reproducible raw metrics and actual handover components are in `metrics.csv.gz`; `coverage.json` retains every deadline miss and first-readiness result. `provenance.json` records input and source hashes.
