# Continuous magnetic hard-iron correction

The OU-II, OU-III and TFG facades continuously estimate residual magnetic calibration, shared by their CARMA, shared-mode and tangent variants. The OU simulators expose the matched ablation as `SF_MAG_CONT_HI=0`; TFG exposes `TFG_MAG_HARD_IRON=0`. The correction starts after magnetic startup acquisition, using statistics accumulated from the first sample.

The AtomS3R OU-II, OU-III and TFG sketches (`atomS3R_ins_kalman_ou2.ino`, `atomS3R_ins_kalman_ou3.ino`, `atomS3R_ins_tfg.ino`) set `mag_continuous_hard_iron = false`; the facade defaults are unchanged. Under heel excitation the additive-offset fit absorbs calibration residuals (axis gain, soft iron, magnetometer latency) as an offset with a horizontal component, and applying it rotates magnetic north until the memory window clears.

## Estimator and observability

The legacy additive mode in `src/tuner/ContinuousMagHardIronEstimator.h` accumulates exponentially weighted magnetometer and yaw-stripped private-observer tilt statistics. With `Abar = mean(R_i)`, it solves the regularized system

```
(I - Abar^T Abar) b = mbar - Abar^T wbar.
```

Eliminating the world field leaves a three-dimensional offset solve. At a fixed attitude, the offset and reference field cannot be identified separately. Roll and pitch excitation provide information; small vessel tilts can make the system poorly conditioned. Multiplicative soft-iron and alignment errors can alias into an additive solution, so a small fitting residual alone does not establish calibration accuracy.

The relative ridge scales with the normal matrix's mean eigenvalue; an absolute ridge remains for weak excitation. The selected vessel configuration uses relative ridge 0.25 and minimum information 0.1 in all three facades. The minimum-information setting controls calibration activation; it is a filter parameter, not a simulation acceptance threshold. Other residual, offset and conditioning checks remain in force.

## Applying the correction

The applied calibration transforms the magnetic measurement; in legacy mode this is just subtraction of the estimated offset. The magnetic reference retains its canonical horizontal direction; only magnitude and dip are adjusted by a delta from the acquired reference. Rotating that reference with the inferred offset would preserve the heading ambiguity the correction is intended to resolve. The estimator reads the private observer and raw magnetometer, not the MEKF state. Nevertheless, its output changes the measurement model: stability admission must account for the actual time-varying reference and covariance and cannot be inferred merely from this separation.

TFG additionally supports `TFG_MAG_REFINE=0` and a staged startup ablation. Its selected refinement start is 30 seconds. The OU startup schedule retains its existing timing.

## Vessel-response evidence

The archived legacy-mode study `tools/startup_ablation.py` replays all eight pinned v1.3.1 28 ft / 5 kn moving-vessel records (±3° heading over 120 s). Five paired calibration draws for each OU family, plus four TFG arms, produce 192 measurements. The sensor draw is held fixed in this startup study. Results and input/binary hashes are in `reports/results/startup_ablation/`; `tools/startup_publication_sync.py` generates the paper tables.

In the default draw, mean yaw RMS falls from 2.241° to 0.708° (OU-II), 2.041° to 0.595° (OU-III), and 1.546° to 0.876° (TFG). Across the five calibration draws, 4 of 40 OU-II pairs and 8 of 40 OU-III pairs worsen.

The separate archived stationary-vessel parameter-selection study varies both sensor and initialization draws, with held-out seeds 101, 1009, 2027 and 3037. Its full measurements are in `reports/results/rao_parameter_tuning/`. These finite-record comparisons do not guarantee improvement on every sea, sensor realization or hull. The current moving-vessel quality gates are documented in [filter-simulator-gates.md](filter-simulator-gates.md).

## Rotation-aware residual calibration

With `mag_continuous_hard_iron` and `mag_gyro_calibration` enabled (the facade
defaults), `GyroMagCalibration` replaces the fixed-heading additive solve.
It uses raw gyro and magnetometer samples. Neither true heading nor navigation
attitude, navigation bias estimates, vessel speed or the declared replay track
enters it. The old additive fit remains available with
`mag_gyro_calibration=false` (`SF_MAG_GYRO_CAL=0` in OU replays,
`TFG_MAG_GYRO_CAL=0` in TFG).

The local magnetic field `f` is transported between magnetic packets by the
relative gyro rotation. The measurement is `m = A f + b + n`, where `A` is a
small residual scale/axis/alignment matrix and `b` is the magnetic offset.
The 18 scalar states are `f` (3), `b` (3), residual gyro bias (3), and `A-I`
(9). Gyro bias is a nuisance parameter of this sensor-only learner; it is not
written into the navigation filter. Magnetic offset and gyro bias have slow
random-walk priors; the calibration matrix is constant. These are estimation
priors, not assertions that physical sensor errors follow those laws.

The first magnetic reading seeds `f`. Its offset/matrix uncertainty is carried
as correlated initial covariance, rather than treating that first reading as
an independently known field. Prediction preserves the full covariance and
the bias sensitivity of the relative rotation. Correction uses an expanded
Joseph update with only three measurement rows, limiting embedded stack use.
Invalid IMU intervals clear the rotation history; duplicate magnetic calls
without new IMU time cannot add information. Innovation, finite-value and
calibration-magnitude checks reject unsuitable updates.

The matrix prior describes small residual errors after ordinary sensor
calibration: diagonal standard deviation 0.00866, symmetric cross-axis
standard deviation 0.00577 and alignment standard deviation 0.01008 rad.
Magnetic white noise is 0.8 uT, initial offset standard deviation 1.155 uT,
and offset random-walk density 0.01 uT/sqrt(s). Gyro white noise is
0.00157 rad/s per 5 ms sample, initial bias standard deviation
0.000504 rad/s and bias random-walk density 1e-5 rad/s/sqrt(s).
Unexcited calibration directions remain prior-dependent: this is not an
unrestricted replacement for a full installation calibration.

Application still waits for live operation, completed magnetic refinement,
finite residuals and private-observer tilt excitation. The offset and inverse
matrix slew together over the existing 45 s time constant. Raw sufficient
statistics evaluate `mean(R C (m-b))` for both the applied and anchored
calibration. Their difference updates only the acquired reference's horizontal
magnitude and vertical component. Its canonical horizontal direction is kept.
The startup acquisition and private gravity observer remain unchanged.

The public hard-iron diagnostic still exposes the legacy tilt-only fit;
`magContinuousCalibration()` exposes the active gyro-based posterior. Disabling
continuous calibration leaves the measurement matrix at identity. The current
AtomS3R sketch opt-outs remain opt-outs. Calibration uncertainty is not jointly
augmented into the navigation covariance; finite replay improvement does not
establish statistical consistency or a stability theorem for this coupling.

`gyro_mag_calibration-test` exercises opposite headings, three noise seeds,
scale/cross-axis/alignment errors, covariance semidefiniteness, missing IMU
intervals, duplicate magnetic packets and absent tilt excitation. Its heading
score uses truth only in the test evaluator. `mahony-frame-test` independently
checks the full NWU/NED basis identity and sustained three-axis motion at four
headings; the PII and heave replay adapters share that conversion.
