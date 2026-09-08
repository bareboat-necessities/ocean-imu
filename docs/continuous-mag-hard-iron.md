# Continuous magnetic hard-iron correction

The OU-II, OU-III and TFG facades estimate an additive body-fixed magnetometer offset continuously. The OU simulators expose the matched ablation as `SF_MAG_CONT_HI=0`; TFG exposes `TFG_MAG_HARD_IRON=0`. The correction starts after magnetic startup acquisition, using statistics accumulated from the first sample.

## Estimator and observability

`src/tuner/ContinuousMagHardIronEstimator.h` accumulates exponentially weighted magnetometer and yaw-stripped private-observer tilt statistics. With `Abar = mean(R_i)`, it solves the regularized system

```
(I - Abar^T Abar) b = mbar - Abar^T wbar.
```

Eliminating the world field leaves a three-dimensional offset solve. At a fixed attitude, the offset and reference field cannot be identified separately. Roll and pitch excitation provide information; small vessel tilts can make the system poorly conditioned. Multiplicative soft-iron and alignment errors can alias into an additive solution, so a small fitting residual alone does not establish calibration accuracy.

The relative ridge scales with the normal matrix's mean eigenvalue; an absolute ridge remains for weak excitation. The selected vessel configuration uses relative ridge 0.25 and minimum information 0.1 in all three facades. The minimum-information setting controls calibration activation; it is a filter parameter, not a simulation acceptance threshold. Other residual, offset and conditioning checks remain in force.

## Applying the correction

The estimated offset is subtracted from the magnetic measurement. The magnetic reference retains its canonical horizontal direction; only magnitude and dip are adjusted by a delta from the acquired reference. Rotating that reference with the inferred offset would preserve the heading ambiguity the correction is intended to resolve. The estimator reads the private observer and raw magnetometer, not the MEKF state. Nevertheless, its output changes the measurement model: stability admission must account for the actual time-varying reference and covariance and cannot be inferred merely from this separation.

TFG additionally supports `TFG_MAG_REFINE=0` and a staged startup ablation. Its selected refinement start is 30 seconds. The OU startup schedule retains its existing timing.

## Vessel-response evidence

`tools/startup_ablation.py` replays all eight pinned v1.2.1 28 ft vessel records. Five paired calibration draws for each OU family, plus four TFG arms, produce 192 measurements. The sensor draw is held fixed in this startup study. Results and input/binary hashes are in `reports/results/startup_ablation/`; `tools/startup_publication_sync.py` generates the paper tables.

The separate parameter-selection study varies both sensor and initialization draws, with held-out seeds 101, 1009, 2027 and 3037. Its full measurements are in `reports/results/rao_parameter_tuning/`. These finite-record comparisons do not guarantee improvement on every sea, sensor realization or hull. Simulation quality thresholds are unchanged.
