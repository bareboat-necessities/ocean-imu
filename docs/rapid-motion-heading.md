# Rapid-motion attitude protection

All ten `.ino` files are accounted for. Seven contain estimators: compass qMEKF,
compass Mahony, OU-II, OU-III, TFG, NLO and PII. They use the same motion and
recovery policy in `src/util/AttitudeProtection.h`. The basic IMU exposes
acquisition data. The magnetic diagnostic also integrates relative gyro rotation
for its axis checks; it does not estimate fused heading. The dummy example has
no estimator.

## Failure mechanism

The [equation audit and causal controls](rapid-motion-math.md) distinguish the
actual qMEKF covariance errors from the acceleration-model and antipodal
correction failures. The corrected qMEKF propagates its right attitude error
with the full SO(3) rotation and exact bias integral, then transports covariance
after finite quaternion injection. These algebraic repairs do not by themselves
remove the acceleration-driven flip.

The common host replay reproduces large attitude errors with exact gyro and
magnetic measurements. Rotation about a displaced sensor produces tangential
and centripetal acceleration. Treating its direction as gravity can corrupt
tilt and the gyro-bias estimate. qMEKF's magnetic-dip adaptation then uses the
corrupted attitude. Local magnetic corrections also have poor capture near a
180-degree heading error. These are filter failure mechanisms; a hardware fault
is not needed to reproduce them.

In OU-III, navigation/attitude cross-covariance can carry the disturbance beyond
the turn: subsequent integral-displacement corrections can move attitude even
after the accelerometer has become quiet. The protection therefore covers
pseudo measurements and recovery transitions as well as the accelerometer.

## Runtime behavior

The guard suspends attitude and bias aiding when angular rate exceeds 1.5 rad/s,
filtered angular acceleration exceeds 2 rad/s², specific-force magnitude differs
from local gravity by more than 20%, or an update gap exceeds 20 ms. Angular
rate threshold uses the raw gyro norm, so a short rate spike also triggers a
hold. Two 30 ms low-pass stages estimate angular acceleration. Release requires
500 ms below the lower rate, angular-acceleration and force thresholds.

Gyro propagation continues with the learned bias. Mahony retains its integral
without adapting it, NLO freezes attitude/bias feedback, and the Kalman filters
zero attitude, gyro-bias and accelerometer-bias gain rows. Kalman navigation
updates and the continuous translational observers still receive the measured
specific force. The private startup/tuner Mahony observers are protected too;
dynamic force cannot seed their initial tilt.

A large discrepancy (35 degrees) starts independent attitude capture; the
15-degree lower threshold provides hysteresis. Capture uses gravity and the
horizontal magnetic direction directly, so the antipodal heading has no
zero-innovation trap. It requires a stable window of at least 750 ms and twelve
new magnetic observations, with no magnetic gap over 120 ms and no acceleration
gap over 20 ms. Gyro-transported held acceleration is not new gravity evidence;
an acceleration dropout holds attitude/bias aiding and clears the capture window.
Acceleration and magnetic observations need not arrive on the same tick.
Angular rate must be below
0.12 rad/s and acceleration within 4% of local gravity. Field strength, horizontal
observability, direction consistency, and the previously established field
strength/dip fingerprint must also pass.

After a dynamic interval the filter also attempts this settled capture when
the discrepancy is small. Extra settling hold is limited to two seconds once
the dynamic latch clears, unless a large error remains pending. Thus ordinary
vessel motion can regain aiding without coming to rest. The three navigation
Kalman filters discard cross-covariance between held attitude/bias and driven
navigation blocks when aiding resumes. This preserves means and marginal
covariances and keeps covariance positive semidefinite.

Capture replaces attitude and its local covariance without clearing navigation
history, learned bias means, or the magnetic world reference. TFG also transports
its bias-error covariance coordinates so the physical body bias covariance is
preserved. Debug output includes `[ATT] hold=… dynamic=… recoveries=…`.

## Acquisition and time

The shared reader now delivers each fresh gyro sample after the first
accelerometer sample, including gyro-only updates. A pending accelerometer
sample is consumed once by discrete Kalman updates. Continuous observers
transport the held force vector into the new body axes between fresh samples.
Calibration capture still counts only fresh acceleration observations.

Compensated BMM150 XYZ and RHALL values, and their freshness decision, come from
the same eight-byte AUX read. A separate M5 cache change no longer gates that
read. Failed reads do not advance the source sequence. Existing stale-data
rejection remains in effect. Compass, NLO and PII preserve the measured elapsed
time across scheduling gaps instead of replacing it with a shorter interval.
The diagnostic's relative gyro integrator also preserves elapsed time, rejects
invalid gyro data, and discards axis-check pairs that would span a gap. Its
magnetic freshness uses acquisition sequence numbers; an invalid first field
cannot poison later comparisons. Acceleration statistics count fresh samples.

## Regression coverage

`tests/common/RapidMotionReplay.h` runs the same cases through seven thin
adapters of the production estimator implementations. The suite is part of
`make all` and covers:

- 180, 360 and 720 deg/s turns; 0, 5 and 15 cm pivot offsets; yaw and combined
  rotation; 200 Hz gyro, 100 Hz acceleration and 25 Hz magnetic delivery.
- A 75 ms scheduling gap, 40 ms unmodelled magnetic latency, and seeded gyro,
  acceleration and magnetic noise with nonzero physical gyro bias.
- 90, 150, 179.99, 180, 180.01 and -179.99 degree heading errors at two absolute
  headings and nonzero roll/pitch, plus noisy 180-degree capture at both headings;
  stale, invalid and strength-disturbed fields.
- Continued moderate roll and wave acceleration after a turn, bias preservation,
  covariance validity, navigation-state preservation and protected startup.
- Accelerometer dropout after a contaminated sample, isolated acceleration
  frames, and a missing initial acceleration stream, including clock wrap and
  recovery after asynchronous fresh acceleration resumes, across all seven backends.
- Exactly antipodal, horizontal-field controls: ordinary local updates remain
  at 180 degrees, while finite-angle capture escapes the equilibrium in all seven.
- Actual acquisition function bodies with controlled driver masks and AUX reads,
  and the diagnostic's integration and freshness code with gaps and invalid data.

Protected rapid cases must remain below 12 degrees of total attitude error,
end below 5 degrees after seven seconds, and keep gyro-bias magnitude below
0.06 rad/s. Every large-angle case must capture exactly once and finish below
3 degrees within three seconds. The matched unprotected control is printed
alongside protected results. These are engineering regression bounds, not a
stability theorem or a hardware accuracy specification.

The [source-bound replay manifest](../reports/results/rapid_motion/manifest.json)
records all 22 protected rapid cases and all fourteen large-angle cases per backend.
The maxima below
include ideal and noisy cases in both categories. Capture error is measured three
seconds after introducing the heading error; rapid tail error is measured
seven seconds after the turn ends. Values are total attitude error in degrees,
rounded upward to two decimals.

| Backend | Peak during rapid replay | Rapid tail | Large-angle capture tail |
| --- | ---: | ---: | ---: |
| qMEKF | 3.57 | 0.53 | 0.56 |
| Mahony | 3.78 | 0.91 | 0.75 |
| OU-II | 4.16 | 1.46 | 0.95 |
| OU-III | 6.60 | 2.66 | 1.61 |
| TFG | 4.08 | 0.60 | 1.09 |
| NLO | 7.06 | 2.80 | 2.11 |
| PII | 3.53 | 2.09 | 1.39 |

## Qualification limits

Hardware turn tests remain necessary. The BMM150 conversion epoch is still
unknown; a host read timestamp does not identify it. A gyro sample cannot
reconstruct motion that occurred during an unobserved gap. The guard cannot
distinguish every translational acceleration from tilt, and a persistent rotated
magnetic field with the same strength and dip is not distinguishable from a
heading error using these sensors alone. Existing OU-III proof results retain
their original ordinary-aiding scope; they do not certify the new guard,
covariance decorrelation or global capture transitions.
