# Where the gravity estimate changes after handover

The first live accelerometer correction can move the attitude several degrees even when the handover supplies the exact direction of gravity. On the high-wave TFG example below, the change takes one 5 ms sample. Each filter uses a fixed world gravity magnitude. What changes is the estimated attitude that projects gravity into the moving body frame. Its effect on the accelerometer can resemble both wave acceleration and sensor bias.

This experiment traces the guarded TFG configuration at `abfd578`, which assigned roughly 2° uncertainty to an accepted fit. Production subsequently adopted 0.5° for an accepted TFG wave fit. Its [separate matched comparison](startup-wave-learning.md#tilt-uncertainty-after-an-accepted-tfg-fit) reports the current results.

## One update at a time

The table follows the default sensor draw on the 8.5 m following JONSWAP record. The exact-tilt arm changes only body down at the handover sample and keeps the estimated yaw, gyro bias and field. The first live sample arrives 5 ms later. All angles are body-down errors in degrees.

| Filter and handover | Installed | After gyro prediction | After accelerometer update | After 1 s, accelerometer |
|---|---:|---:|---:|---:|
| OU-III, deployed | 0.2896 | 0.2893 | 1.0361 | 0.6431 |
| OU-III, exact tilt | 0.0000 | 0.0004 | 0.8114 | 0.4911 |
| OU-III, measured wave seed | 0.2896 | 0.2893 | 0.2689 | 0.5375 |
| TFG, deployed | 0.2898 | 0.2898 | 4.6633 | 4.8001 |
| TFG, exact tilt | 0.0000 | 0.0003 | 4.6424 | 4.7984 |
| TFG, measured wave seed | 0.2898 | 0.2898 | 0.2661 | 1.6652 |

Both deployed handovers start the wave-acceleration mean at zero. The gyro prediction carries the supplied tilt forward on this sample. The subsequent correction must distribute the measured specific force between the allowed error states under the current covariance. The gyro and accelerometer stages are logged separately, before the later magnetic update.

Across all records and three sensor draws, first-sample tilt RMS starting from exact tilt changed between prediction and acceleration correction as follows. OU-III 0.009° → 0.520° (0 of 48 above 1°). TFG 0.018° → 2.562° (32 of 48 above 1°).

## First correction with the production 0.5° setting

On the default sensor draw across all 16 records, the TFG first accelerometer correction gave 2.6514° pooled body-down RMS with the 2° prior and 2.1753° with the 0.5° prior. Starting both from exact tilt gave 2.6463° and 2.1530°. One second after handover, the deployed arms measured 2.4645° and 2.5386°, respectively. Giving the fit more weight reduces the initial correction, while an appreciable attitude transient remains.

On the 8.5 m following JONSWAP record, the new production handover starts at 0.2898°, changes to 3.8427° after the first accelerometer update, and yields 0.07857° late pitch RMS. Starting the same filter from exact tilt gives 3.7762° on the first correction and 0.07760° late pitch RMS. All nine early and late window metrics in the fitted instrumented arm matched the ordinary production executable on every record. [Fitted trace, control comparison and run inventory](../reports/results/startup_learning/tfg_tilt_uncertainty/fitted_trace)

## What persists

This table pools the 16 records on the default sensor draw. Each cell is the sample-weighted RMS for the previous startup, deployed startup and deployed startup with exact tilt, in that order. The final window contains 180,000 samples per record.

| Filter and time since power-on | X displacement (m) | Y displacement (m) | 3-D displacement (m) | Accel-bias vector (m/s²) |
|---|---:|---:|---:|---:|
| OU-III, 40–60 s | 0.87617 / 0.84417 / 0.84314 | 0.49929 / 0.45737 / 0.45731 | 2.01414 / 1.99516 / 1.99482 | 0.05465 / 0.05465 / 0.05465 |
| OU-III, 300–1,200 s | 0.29893 / 0.29896 / 0.29896 | 0.25029 / 0.25026 / 0.25026 | 0.43372 / 0.43372 / 0.43371 | 0.03579 / 0.03578 / 0.03578 |
| TFG, 40–60 s | 5.01026 / 5.02038 / 5.01835 | 2.20682 / 2.11475 / 2.11575 | 5.94345 / 5.91572 / 5.91430 | 0.05465 / 0.05465 / 0.05465 |
| TFG, 300–1,200 s | 0.30743 / 0.30736 / 0.30736 | 0.24885 / 0.24890 / 0.24887 | 0.41826 / 0.41825 / 0.41823 | 0.05394 / 0.05398 / 0.05396 |

On the 8.5 m following-sea example, late TFG pitch RMS changed from 0.08604° to 0.09180° with the production handover. Giving the previous startup exact tilt produced 0.08557°, and giving the deployed startup exact tilt produced 0.09224°. Supplying the exact magnetic field at the deployed handover produced 0.11743°, and supplying both exact tilt and field produced 0.11794°. The worse late pitch persists with perfect initial tilt. The production accelerometer-bias mean is held during the early window. The reference-refinement and bias-release schedules govern when the live filter can adjust it. OU-III opens that gate at roughly 120 s on these records.

## Carrying wave acceleration into the live update

At accepted handovers we also tried a sensor-only seed. It subtracts calibrated acceleration bias and the estimated turn acceleration from the endpoint specific force, rotates the result into NED using the learned attitude, and adds world gravity. The live wave acceleration takes that value once. Velocity, displacement and bias means retain the production initialization. This tests what the initial zero mean costs.

Each pair below compares deployed TFG with the measured wave seed on the same 16 records. Each draw varies the sensors on the same wave histories.

| Sensor draw | Early X (m) | Early Y (m) | Late X (m) | Late Y (m) | Late accel-bias vector (m/s²) |
|---|---:|---:|---:|---:|---:|
| default | 5.02038 → 1.46809 | 2.11475 → 0.51307 | 0.30736 → 0.31127 | 0.24890 → 0.24598 | 0.05398 → 0.06755 |
| 7 | 4.89119 → 1.33595 | 2.18868 → 0.64637 | 0.30024 → 0.30392 | 0.24482 → 0.24141 | 0.05963 → 0.05606 |
| 41 | 5.00818 → 1.39065 | 2.08547 → 0.62676 | 0.30584 → 0.30664 | 0.24452 → 0.23430 | 0.08838 → 0.06411 |

On the 8.5 m following-sea default draw, early TFG pitch changes from 5.5257° to 1.0693°. Late pitch changes from 0.0918° to 0.7883°. The early benefit comes with a repeatable late X penalty and draw-dependent bias changes. The measured wave seed stays in the replay experiment.

A separate 48-run test at `c3481ce` retained the zero wave-acceleration mean and increased its initial covariance. On the default 8.5 m following JONSWAP record, the production reference had 5.5257° early and 0.0918° late pitch RMS. Using the measured specific-force covariance gave 1.9102° early and 1.0519° late. An isotropic covariance based on the same total variance gave 1.9249° early and 1.0535° late. Both variants were rejected. The production 0.5° change adjusts the attitude uncertainty alone. [All 48 runs and the covariance definitions](../reports/results/startup_learning/tfg_tilt_uncertainty/rejected_wave_covariance)

## Scope and evidence

Completed 672 full 1,200 s replays on production source `abfd57856fc008556bf64e004d82f4b99f5e90f7` with g++ (Debian 14.2.0-19) 14.2.0. Both v1.4.0 archive directions contribute eight records. Sensor draws were default, 7, 41. Physical wave histories are reused across draws. Observed handover times ranged from 31.83 to 41.85 s. Windows are measured from power-on. The seed-7 0.27 m following JONSWAP run hands over at about 41.85 s, so its 40–60 s window includes approximately 1.85 s of startup. All default-draw early windows and all late windows occur after handover.

For both families and both sea directions, the deployed instrumented control matched all nine RMS values from a separately compiled, unmodified production executable in all 160 default-draw windows at the precision emitted by the simulator. [Production-control comparison](../reports/results/startup_learning/gravity_oracle/production-comparison.json)

**Structures preserved.** Each matched run uses the same 200 Hz IMU stream, measurement ordering, startup deadline and fallback policy, live attitude and bias equations, magnetic-reference refinement, wave tuner, and covariance updates. The underlying physical history continues through handover. The production flag selects the previous and deployed controls.

**Relaxations introduced.** The exact arms replace only the named handover quantities with simulation truth. The measured-wave arm changes one initial state mean using only sensor data. Trace instrumentation samples truth after handover solely for scoring. These finite replay comparisons establish behavior on the listed records. They supply no stability counterexample or uniform observability bound. Keep release gates tied to the ordinary production executable.

[Reproduction command](../tests/startup_learning/README.md#tracing-gravity-through-the-live-handover) · [Handover, updates, window scores and provenance](../reports/results/startup_learning/gravity_oracle)
