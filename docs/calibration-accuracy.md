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

Magnetic capture averages distinct readings over 400 ms, rejects excessive
rotation within a window, and retains 160 averaged observations. Independent
verification uses 140 new means with frozen coefficients. Normal guided
collection takes about 70 seconds, followed by about one minute of verification;
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
clean field directions, including appreciable pitch and roll.

| Metric | v2.3.2 | Current |
| --- | ---: | ---: |
| Gyro bias RMS error, rad/s | 0.000284 | 0.000105 |
| Magnetic direction RMS, 30 µT noisy field | 1.067° | 0.293° |
| Magnetic direction RMS, 50 µT noisy field | 0.545° | 0.172° |
| Magnetic direction RMS, 50 µT quiet field | 0.150° | 0.095° |
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
repeated and frozen registers, rapid turns, lost gyro input, realistic magnetic
noise during gyro holds, and optional-preset rollback.

Run `make -C tests/imu_calibrate test` for the current checks. To reproduce the
reference from the Git tag, run `bash tests/imu_calibrate/compare_release.sh`
from the repository root. This extracts the release into a temporary directory
and builds the same gyro/magnetometer history driver plus the release's own
accelerometer campaign. `make all` is the complete repository check.
