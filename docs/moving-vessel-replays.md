# Moving-vessel simulation inputs

Every external oceanography-waves-lib replay uses `v1.3.0`, archive
`sim-data-files-vessel-rao-28ft-5kn.zip`, SHA-256
`324e28cf77cbef9dc1c4bed6e0d8ce54c5c805769384e5997d97a4824882baf3`.
`tools/sim_dataset.py` owns the pin, validates each cached member, and records
individual hashes for the eight scored JONSWAP/PM-Stokes inputs. Filenames,
incident heights and spectra remain compatible with existing harnesses.

The 28 ft hull moves at 5 kn (2.572222 m/s). In north/west/up coordinates,
heading is `psi(t) = 3 degrees * sin(2*pi*t/120 s)`, with velocity
`U*(cos(psi), sin(psi), 0)`. Position is its analytic primitive from time zero;
acceleration is its derivative. The upstream model samples the incident wave
field along that path, applies estimated RAOs with signed encounter frequency
and changing heading, and exports total CG travel, velocity and acceleration.
It includes rotating response terms and full yaw/pitch/roll body rates. It
uses quasi-steady estimated RAOs, without hull-specific forward-speed
hydrodynamics or a propulsion/rudder model.

## Scoring convention

The estimators reconstruct bounded wave motion about a local equilibrium.
They have no external velocity/position aid for absolute cruise navigation.
`UnderwayReference.h` subtracts the declared track from **reference position
and velocity only**, after running the filter. The 3D percentage denominator
uses that residual displacement, not kilometres of forward travel. This is
an analytic truth-side coordinate definition, never an output-dependent fit
or a sensor correction. Horizontal travel is parsed in double precision
before subtraction to retain small wave motion late in a record.

Recorded specific force, angular rate, attitude, magnetic truth and total
acceleration remain intact. Steering acceleration is therefore an actual
low-frequency disturbance to the wave estimator. `acc_ref_*` continues to
mean total acceleration; it includes the track acceleration in addition to
the derivative of `vel_ref_*`. No known track state is supplied to a filter.
Standalone synthetic diagnostics can set `W3D_REFERENCE_MOTION=stationary`;
the low-level runner also defaults to total reference coordinates unless a
replay entry point explicitly selects the underway convention.

Phase surrogates remove the declared p/v/a/yaw track, randomize the residual
wave chain, restore the same track, and rebuild body IMU fields. Wave-height
scaling likewise scales only residual motion. Both sea-state crossfade
endpoints retain the identical track, so the existing exact derivative cross
terms preserve continuous total motion through a transition.

## Evidence and regression limits

The eight full records retain 200 Hz sampling, 1200 s duration, default noise
and seeds, and trailing 900 s W3d scoring. AHRS uses its trailing 60 s window.
Production and interval-rate profiles are measured separately. Limits are
rebaselined from this changed physical experiment, with a 0.5% relative margin
and upward rounding; a change of dataset is not an estimator improvement.
Nonfinite, partial, missing or duplicate records remain failures.

Older evidence naming v1.2.3 describes the stationary-vessel experiment. Its
hashes and numerical conclusions must not be relabelled as underway evidence.
The full OU evidence workflow regenerates source-bound statistical bundles
on the branch, then synchronizes the articles; it rejects stale bundles.
