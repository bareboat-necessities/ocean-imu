# Calibration accuracy checks

The reference is release **v2.3.2**, commit
`2ba31e4e255a2b03f014cfbb9c9352d58d49f791`. The comparison executes its
actual calibration headers, rather than treating the older six-pose method
labelled `OLD` in the accelerometer test as the released procedure.

The accelerometer's released full-matrix fitter, five useful seconds per pose,
thermal qualification, additional poses and later rechecks remain intact.
The current progress display counts accepted data and diagnostics do not wait
for serial output. The 30-seed accelerometer campaign has identical rounded
acceptance and error metrics in all ten scenarios. Typical-session results:

| Accelerometer metric | v2.3.2 | Current |
| --- | ---: | ---: |
| Acceptance | 30/30 | 30/30 |
| Median bias error, m/s² | 0.00231 | 0.00231 |
| 90th-percentile bias error, m/s² | 0.00556 | 0.00556 |
| Median held-out vector RMS, m/s² | 0.00406 | 0.00406 |

Gyro capture averages fresh readings into 25 ms observations while applying
stillness gates to the raw stream. Two-second magnetic means establish and
check a fixed yaw reference with the existing displacement threshold. Sixty
simulated low-power BMM150 noise histories complete an uninterrupted hold
within eight seconds. Thermal fitting remains separately qualified.

Magnetic capture retains the mean and population covariance of at least twelve
distinct readings over at least 400 ms, and keeps 160 such observations.
The covariance preserves the corrected field energy through motion inside each
window; hand turns require no angular-rate cutoff. Slow sampling extends the
window rather than discarding partial data. Independent verification uses
140 new observations with frozen coefficients. At 20–30 Hz, guided collection
takes about 70–100 seconds, followed by about 60–85 seconds of verification;
the display continues to show data and coverage progress throughout.
This is an accuracy/usability tradeoff: the extra time measures new data.
Hardware precision settings are optional; software averaging is also used
with driver defaults. Noise, coverage, field-change and information thresholds
are unchanged.

The comparison uses 30 seeds, known symmetric soft-iron distortion and hard-iron
offset, smooth three-dimensional turns, repeated 200 Hz polls of 20/30 Hz
magnetic readings, 0.3 µT quantization, and an additional correlated-noise
component. The noisy magnetic cases use per-axis white noise of 1/1/1.4 µT;
the quiet case uses 0.3/0.3/0.3 µT. Angular error is measured on 1,000 unseen
clean field directions, including appreciable pitch and roll. The driver
sequences its random draws explicitly, so GCC and Clang replay identical
histories. The 90/90 counts belong to those histories: three alternative draw
orders of the same noise model each leave one noisy 30 µT session failing
independent verification (89/90).

| Metric | v2.3.2 | Current |
| --- | ---: | ---: |
| Gyro bias RMS error, rad/s | 0.000284 | 0.000105 |
| Magnetic direction RMS, 30 µT noisy field | 1.067° | 0.292° |
| Magnetic direction RMS, 50 µT noisy field | 0.545° | 0.167° |
| Magnetic direction RMS, 50 µT quiet field | 0.150° | 0.062° |
| Successful magnetic fits | 90/90 | 90/90 |
| Successful independent magnetic checks | Not present | 90/90 |

These are simulations, not device accuracy measurements. They do not establish
absolute north accuracy or identify physical sensor mounting rotation. The
full wizard still needs validation on real recordings and hardware. In
particular, a failed independent check must never be bypassed to obtain a pass.

`calibration_accuracy-test` enforces improvement against the release baseline
for all three magnetic noise/field groups and gyro bias. It also verifies
that a later field step is rejected and that verification never changes fitted
coefficients. `test_release_accuracy.py` checks the deployed accelerometer
against the committed release reference. `calibration_workflow-test` covers
repeated and frozen registers, rapid turns, display pauses, realistic magnetic
noise during gyro holds, and optional-preset rollback. An algebraic regression
checks that window moments give the same corrected field energy as the raw
readings, including an off-diagonal calibration matrix and nonzero bias.

`mag_hand_motion-test` exercises 120 noisy/quiet sessions with continuous and
stop/start hand motion, 75 ms display pauses every 250 ms, and 10/20/30/100 Hz
magnetic sampling. All complete fitting and frozen independent verification,
and all reject a subsequent field step. Held-out direction RMS is 0.322° in the
noisy 30 µT group and 0.057° in the quiet 50 µT group. The slowest 10 Hz case
takes 195 seconds for capture and 170 seconds for verification, within the
220/180-second limits. These worst-case durations are not the normal 20–30 Hz
guidance times.

Run `make -C tests/imu_calibrate test` for the current checks. To reproduce the
reference from the Git tag, run `bash tests/imu_calibrate/compare_release.sh`
from the repository root. This extracts the release into a temporary directory
and builds the same gyro/magnetometer history driver plus the release's own
accelerometer campaign. `make all` is the complete repository check.
