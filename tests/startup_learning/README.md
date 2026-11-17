# Isolated wave-startup experiment

This host-only experiment reads the current v1.4.0 following-sea and head-sea
records. `StartupLearner` has no access to truth, filenames, wave height, vessel
speed, inclination, declination or field magnitude. It consumes calibrated
body-axis accelerometer, instantaneous gyroscope and timestamp-matched
magnetometer samples. No production filter includes this code.

`ShippingBaseline` runs the unchanged shared `ProxyStartupFusion` through the
OU-III wrapper and stops at its **first live handover, before a live magnetic
correction**. The test-only private-access shim reads the acquired magnetic
reference; it does not modify any estimator equation or API. Bias and
quaternion seeds are read from the actual live state.

```sh
make -C tests/startup_learning
tests/startup_learning/startup-test
python3 tools/startup_learning.py --seeds 1,2,3 --offsets 0,7
```

The runner verifies both archive/CSV pins using `tools/sim_dataset.py`, makes
one immutable noisy power-on history per trial, and feeds all methods that
same prefix. Only sensor generation and scoring access the CSV truth. The
baseline is scored at its own handover epoch; candidate `matched` rows use
that exact epoch. Separate `fixed` rows use 10, 15, 20, 25, 30 and 35 seconds.
No estimator receives samples beyond the horizon being scored.

## What is being estimated

| Quantity | Startup treatment |
|---|---|
| Gravity magnitude | Known, 9.80665 m/s² for these simulations; never fitted |
| Gravity direction | Unit vector, two physical orientation degrees of freedom |
| Attitude | Current-epoch body-to-magnetic-NED quaternion, all three axes |
| Local magnetic field | Unconstrained 3-D vector in the initial gyro frame; transformed to magnetic NED |
| Gyroscope bias | Regularized finite-window fit, retaining a prior in weak directions |
| Accelerometer bias | Calibration-prior seed; **not claimed learned** in a short moving window |
| Angular rate | Latest measured body rate minus gyro-bias seed, in rad/s |
| Next prediction history | Last raw gyro sample and its acquisition timestamp |
| Position, velocity, wave state | Outside this experiment; cannot be initialized from these attitude statistics alone |
| True north | Not learned; external declination is needed to convert magnetic north to geographic north |

Existing hard/soft-iron and axis calibration remain prerequisites. A coherent
magnetic offset can look like a different local field on a short, weakly
excited record. This experiment explicitly measures that degradation.

## Methods

All candidates propagate a relative body-to-initial-frame quaternion from the
first sample, with no fictitious prediction before that sample. Instantaneous
gyro endpoints use trapezoidal integration. The terminal attitude is therefore
at the **last acquisition epoch**, not the mean heading in the window.

* `uniform`: average gyro-transported force and magnetic observations using
  the gyro calibration prior.
* `tapered`: replace the force mean by the article's normalized
  `30 x² (1-x)²` quadrature. Its endpoint value and derivative vanish, which
  suppresses bounded-displacement contamination. Only the direction of this
  intermediate force statistic is used; gravity magnitude is fixed.
* `magfit`: estimate gyro bias by minimizing the scatter of the transported
  magnetic samples, eliminating the unknown constant field analytically.
  Every bias trial repropagates the complete gyro history.
* `joint`: additionally penalize the tapered linear slope of transported
  force. This informs drift about the magnetic axis. Its 35-second model
  standard deviation is 0.005 m/s³, enlarged as `(35/T)^3` for shorter windows.
  It is a finite-window wave-contamination assumption, not thousands of
  independent gravity measurements or a proof that wave acceleration averages
  exactly to zero.
* `joint_cal`: the same joint fit also estimates a constant residual hard-iron
  vector, with a 1.2 µT calibration prior. The field is still eliminated
  analytically. This is a separate ablation, because fitting an unnecessary
  calibration parameter costs accuracy when the magnetometer is already clean.

The gravity direction and horizontal magnetic projection define a right-handed
magnetic-NED frame. That frame is composed with the **terminal** relative
rotation. Magnetic magnitude and inclination are outputs, not inputs.

## Measurement contract and limits

Sensor axes are body NED; quaternion maps body NED to magnetic world NED.
Specific force is negative down at rest. Magnetic values are microtesla, rates
rad/s and acceleration m/s². Magnetometer samples run at 25 Hz, IMU at 200 Hz.
An arriving magnetic sample must first be associated with its acquisition
epoch; arbitrary delivery-time reuse is not supported. Invalid samples,
duplicate/backward timestamps, gaps above 0.1 s and windows beyond 35 s are
rejected without mutating the history.

Noise profiles separate ideal sensors, noisy/biased IMU with calibrated
magnetometer, the repository's residual magnetic calibration model, and a
10 µT body-frame field step at 15 s. All three sensor streams have separate
seed streams; residual calibration and measurement noise are also separate.

The development set uses seeds 1–3 at offsets 0 and 7 s. The recorded validation
set uses seeds 11–15 at 67, 307 and 607 s, after freezing both fits. The focused
CI workflow repeats both directions and checks paired yaw improvement, tilt,
gyro-bias error and complete fixed-horizon outputs. Run it locally with:

```sh
python3 tools/startup_learning.py --seeds 11,12,13,14,15 \
  --offsets 67,307,607 --profiles calibrated,residual --check \
  --output reports/results/startup_learning/validation
python3 tools/report_startup_learning.py reports/results/startup_learning/validation
```

The optional report renderer needs Matplotlib; the experiment and checks use
only a C++20 compiler, Eigen and the Python standard library. `metrics.csv.gz`
contains six checkpoints plus exact matched-time comparisons, including the
actual quaternion, field, bias/rate seeds and conditional gyro-bias covariance.

Readiness is a diagnostic heuristic based on elapsed time, horizontal field,
magnetic scatter, force norm, weighting sensitivity, fresh magnetic data and
four-sigma consistency with the supplied calibration priors. The report counts
false readiness against an explicit target of yaw error ≤2° and tilt ≤1°.
A small residual is not an accuracy certificate. In particular, an unavailable
heading, a deadline fallback, and a verified accurate seed must remain distinct
when this is eventually integrated.

The field-step control revealed that a constant-field optimizer can absorb
interference into a very large gyro/hard-iron correction. Such fits remain in
the evidence but cannot be marked ready. The fit equations and accuracy
settings remain fixed; the prior-consistency guard is an additional model
rejection, not a truth-based choice between candidate attitudes.

The current gyro-fit covariance is conditional on the calibration and wave
model. It is not a complete live-filter covariance: attitude, reference,
calibration and wave uncertainties are correlated. A future integration must
transfer those correlations, qualify coverage and continue the gyro clock
without replaying the last interval. This PR deliberately does not integrate
the candidate into any live filter.
