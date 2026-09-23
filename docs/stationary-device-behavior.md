# Stationary heave and magnetic heading

## AtomS3R virtual-constraint cadence

The OU-II, OU-III and TFG full marine INS sketches select the existing fixed
15 ms virtual-constraint cadence. The library default remains tau-scaled for
existing simulations and applications that select that policy. The selected
cadence is passed to the existing pseudo-noise laws; fitted coefficients,
measurement gates and bias-learning gates are unchanged. Bounds and adaptation
transients mean that equal steady information-rate formulas do not imply
identical individual Kalman gains.

A sparse virtual correction changes the raw displacement estimate at the
correction instant. With a large residual acceleration and a long quiet-water
time scale, periodic corrections can create a sawtooth even before output
detrending. More frequent corrections reduce those individual state jumps.
The sketches do not clip displacement, lock it to zero, or introduce an output
smoother. TFG retains fractional scheduler credit using the common OU scheduler.

The nominal 15 ms service period is not a guarantee of hardware throughput.
A loop that runs more slowly cannot deliver that rate. Check actual sample
timestamps and processing rate on the target; frequent virtual corrections
increase compute demand compared with a long tau-scaled interval.

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
Sensors are generated from truth; no truth state is
injected into the estimator.

The tests cover stationary heading with axial gyro bias, large biased/noisy
stationary raw-position corrections with both cadence policies, and an
uninterrupted rest-to-wave-to-rest history. The latter has a twice-continuously
differentiable envelope and acceleration obtained from the displacement's
second derivative. Wave error and return-to-rest behavior are checked without
resetting the estimator or applying a detrender. The heading test distinguishes
handoff, the transient bias-learning interval, and the final stationary minute.
Existing turning-acquisition and vertical-field rejection tests remain separate.

## Scope and limits

These are deterministic software replays, not a hardware validation. They
establish a reproducible scheduler-induced ripple mechanism and a TFG startup
heading defect; they cannot identify every component of a particular device
trace without that trace.

The deliberately large residual-acceleration fixture still exposes a slow raw
TFG startup displacement excursion of about 11.5 m despite much smaller
per-correction teeth. Reducing the teeth does not resolve that separate slow
bias/handoff error. Its magnitude is reported, not clipped, detrended or
labelled a passing absolute-heave result. No regional or global stability
claim follows from these tests.
