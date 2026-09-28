# Stationary heave and magnetic heading

## AtomS3R virtual-constraint cadence

The deployed full marine INS sketches do not all use the same
virtual-constraint (pseudo-measurement) cadence:

| Sketch | Call | Deployed cadence |
| --- | --- | --- |
| `atomS3R_ins_kalman_ou2` | `ff.setTauScaledPseudoUpdateCadence(false)` | fixed 15 ms |
| `atomS3R_ins_kalman_ou3` | `ff.setTauScaledPseudoUpdateCadence(false)` | fixed 15 ms |
| `atomS3R_ins_tfg` | `fusion_.setTauScaledPseudoCadence(true)` | tau-scaled |

The library default of all three filters remains tau-scaled for simulations
and applications that do not select a policy. For TFG the fixed 15 ms cadence
is a comparison configuration in the native regression, not the deployed one.
The selected cadence is passed to the existing pseudo-noise laws; fitted
coefficients, measurement gates and bias-learning gates are unchanged. Bounds
and adaptation transients mean that equal steady information-rate formulas do
not imply identical individual Kalman gains.

A sparse virtual correction changes the raw displacement estimate at the
correction instant. With a large residual acceleration and a long quiet-water
time scale, periodic corrections can create a sawtooth even before output
detrending. More frequent corrections reduce those individual state jumps.
The sketches do not clip displacement, lock it to zero, or introduce an output
smoother. TFG retains fractional scheduler credit using the common OU scheduler.

The nominal 15 ms service period of the OU sketches is not a guarantee of
hardware throughput. A loop that runs more slowly cannot deliver that rate.
Check actual sample timestamps and processing rate on the target; frequent
virtual corrections increase compute demand compared with a long tau-scaled
interval.

## TFG magnetic acquisition and handoff

The gravity-only startup observer cannot correct yaw drift from axial gyro
bias. Its magnetic acquisition gauge must therefore not remain fixed while
its unobservable yaw continues integrating until handoff.

For the deployed acquisition mode, each magnetic sample is rotated using the
independent proxy tilt and its own magnetic north alignment. The window learns
horizontal field strength and dip without attenuation from an actual turn or
from gyro-only yaw drift. Before handoff the yaw alignment is refreshed from
valid magnetic packets. Refinement uses the current magnetic heading while
retaining the core's tilt. The horizontal-field validity threshold and existing
acquisition, refinement and accelerometer-bias gates remain in force.

The optional joint startup hard-iron solve still uses its independent
world-frame identification model: an attitude inferred from the same magnetic
sample must not be treated as independent information for that bias fit.
That non-default mode does not receive the same drift-invariant window model.
The separate continuous hard-iron estimator and its information gates are
unchanged. OU-II/OU-III heading logic and the yaw-free startup compass output
are unchanged.

## Native regression

`stationary_device-test` in each of `tests/kalman_ou_ii`,
`tests/kalman_ou_iii` and `tests/kalman_tfg` exercises the full startup and live
orchestrator at 200 Hz IMU and 25 Hz magnetic sampling. It is run by each
existing native suite. Its additive makefile is included only for the device
regression; the original simulation build recipes remain byte-for-byte intact.
Sensors are generated from truth; no truth state is injected into the
estimator, and raw position is checked before any detrender. Every replay
runs from power-on through handoff without resetting the estimator.

The regression has three tiers:

1. **Stationary heading** at each family's deployed cadence, with axial gyro
   residual 0.025 rad/s. It distinguishes handoff, the transient
   bias-learning interval, and the final stationary minute.
2. **Supported operating envelope (acceptance).** Vertical accelerometer
   residual of +0.01 and -0.01 m/s^2, axial gyro residual 0.002 rad/s and
   0.0148 m/s^2 accelerometer noise, at the deployed cadence. The residuals
   are about twice the simulated post-calibration accelerometer p90 bias
   error and twenty times the gyro RMS error in
   [calibration-accuracy.md](calibration-accuracy.md). The stationary replay
   requires peak raw |p_z| below 0.3 m, below 5 cm after 240 s, correction
   steps below 2 cm and a settled heading. The rest-to-wave-to-rest replay
   (0.3 m, 5 s wave with a twice-continuously differentiable envelope)
   requires wave error below 0.1 m RMS and, after return to rest, raw |p_z|
   below 2 cm peak and 1 cm RMS with small correction steps.
3. **Extreme stress**, outside the supported envelope: vertical accelerometer
   residual 0.2 m/s^2 (beyond the commissioned 0.13 m/s^2 per-axis residual
   envelope), axial gyro residual 0.01 rad/s and 0.12 m/s^2 noise, with both
   cadence policies. It checks finiteness, that the fixed cadence reduces
   correction teeth, and that raw heave recovers below 0.3 m once the
   accelerometer bias is learned. The rest-wave-rest cadence comparison with
   0.03 m/s^2 residual also belongs here.

Existing turning-acquisition and vertical-field rejection tests remain separate.

## Scope and limits

These are deterministic software replays, not a hardware validation. They
establish a reproducible scheduler-induced ripple mechanism and a TFG startup
heading defect; they cannot identify every component of a particular device
trace without that trace.

The stationary startup excursion is a bounded transient, not drift. After
handoff the accelerometer-bias state is held (`freeze_acc_bias_until_live`,
the hold during magnetic refinement, then `acc_bias_unlock_sec`). An
unlearned vertical residual therefore integrates into raw displacement until
the bias is released, then is learned within about ten seconds and raw heave
returns to the millimetre to centimetre level. The peak scales with the
residual: in the TFG supported-envelope replay it is about 0.12 m at
0.01 m/s^2. The extreme stress fixture reaches about 12 m (tau-scaled) and
11.5 m (fixed 15 ms) raw TFG displacement for the same reason, and recovers
to under 0.1 m afterwards. That excursion is a stress response to a residual
outside the supported envelope, not a normal device result, and is reported
rather than clipped or detrended. No regional or global stability claim
follows from these tests.
