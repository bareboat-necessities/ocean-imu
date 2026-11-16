# Heading oscillations, horizontal wave motion and bias

The v1.4.0 replay exposes physical steering acceleration in every tested horizontal filter. The scoring target is wave displacement after removing the entire curved cruise track, while the filters receive the total specific force and have no separate cruise-motion state. The characteristic 20 s displacement error and 10 s body-bias ripple are predominantly responses to this unmodelled motion. They are not explained by rotating a displacement axis alone.

## What the generator and scorer mean

The released generator uses fixed geographic axes: north, west, up. At constant cruise speed U and heading ψ, the cruise velocity and acceleration are

```math
v_c = U(\cos\psi,\sin\psi,0),\qquad
a_c = U\dot\psi(-\sin\psi,\cos\psi,0).
```

For U = 5 kn and ψ = ±3° sin(2πt/20 s), the lateral acceleration peaks at **0.04231 m/s²**. The bounded lateral track has **0.42856 m amplitude** about its cycle mean (0.857 m peak to peak). Relative to the previous 120 s steering cycle, acceleration is six times larger and the forcing has moved from 0.00833 Hz to 0.05 Hz, where the filters reject much less of it. An estimator can also amplify it through its coupled attitude/motion updates; its erroneous displacement is not bounded by the physical track amplitude.

`MovingVesselRao` in oceanography-waves-lib correctly differentiates the total track and wave response, including rotating-response terms. In this repository, `W3dSimulationRunner` subtracts the declared track from displacement and velocity truth, but retains total acceleration and body IMU measurements. Thus the scored displacement/velocity and acceleration columns intentionally belong to different decompositions. A desired wave-only error is meaningful, but it is not an inertial-integration error against one p/v/a history.

The filter output and scored displacement are world vectors. A boat-frame sway product would project **wave-relative** displacement onto the current lateral axis. That projection can legitimately change with heading; it does not remove physical turning acceleration.

Bias columns are sensor-body quantities in the harness basis `(x,y,z) -> (y,x,-z)`. In particular, plotted bias Y corresponds to the generator body X, whereas displacement Y is geographic west. Bias Y and displacement Y are not the same axis. The TFG gallery also had an independent presentation bug: its plot titled estimated bias read the true-bias columns. It now plots true and estimated body bias separately.

## Paired controls

All eight horizontal Kalman families were rebuilt from `5426da8a394b33f5f473f70faa05247f832f85fa`, using GCC 14.2 and the existing per-target host flags. Each control keeps the same 1200 s wave response, roll/pitch history, sensor-noise seed and initialization protocol. Scores use all 180,000 samples of the final 900 s. Both directions use the pinned, checksummed v1.4.0 inputs.

| Control | Translation | Yaw | Meaning |
| --- | --- | --- | --- |
| Released | Total cruise + waves delivered; wave-only p/v scored | Original 20 s cycle | Current production replay |
| Waves + yaw | Complete cruise p/v/a chain removed; body force updated | Original cycle and instantaneous gyro preserved | Separates physical turning force from rotation |
| Waves + fixed yaw | Same wave p/v/a as the preceding control | Mean heading; analytic yaw-rate contribution removed from gyro | Isolates rotation itself |

The latter two controls are deliberately different physical histories. They are not a deployable truth-based correction, and they do not regenerate the hydrodynamic response for a different course. Their p/v/a, specific force, attitude and instantaneous angular rate are checked independently in regression tests. No acceptance limits or filter defaults are changed.

### Y displacement RMS, 0.27 m incident JONSWAP [m]

| Filter | Following: released | Following: waves + yaw | Head: released | Head: waves + yaw |
| --- | ---: | ---: | ---: | ---: |
| OU-II | 0.12804 | 0.01230 | 0.07241 | 0.00617 |
| OU-III | 0.12162 | 0.01265 | 0.08928 | 0.00813 |
| TFG | 0.24759 | 0.01645 | 0.11033 | 0.00784 |
| CARMA | 0.04560 | 0.00447 | 0.02387 | 0.00231 |
| Shared QMEKF | 0.04532 | 0.00416 | 0.02387 | 0.00226 |
| TFG-CARMA | 0.04560 | 0.00447 | 0.02387 | 0.00231 |
| Shared TFG | 0.04532 | 0.00417 | 0.02387 | 0.00226 |
| Shared tangent | 0.04532 | 0.00417 | 0.02387 | 0.00226 |

### Larger-wave check, 1.5 m incident JONSWAP [m]

| Direction / filter | Released | Waves + yaw | Waves + fixed yaw |
| --- | ---: | ---: | ---: |
| following / OU-III | 0.63684 | 0.07730 | 0.07607 |
| following / Shared QMEKF | 0.43217 | 0.03741 | 0.03623 |
| head / OU-III | 0.68343 | 0.05297 | 0.05161 |
| head / Shared QMEKF | 0.41747 | 0.02799 | 0.02763 |

Keeping versus removing yaw changes Y RMS by at most **0.00136 m** in these cells once the cruise acceleration is absent. The large shared degradation is therefore attributable primarily to the cruise-turn forcing, not to the mere presence of yaw rotation. This is a deterministic attribution study, not Monte Carlo evidence for a new estimator.

## The bias ripple

The 10 s sinusoidal component of body-Y bias error falls by **94.5–99.6%** across the measured cells when cruise acceleration is removed but yaw is retained. It is coherent with the steering forcing and its second harmonic. Products of yaw angle and yaw rate contain a 2× frequency component; the experiment identifies the forcing, but does not uniquely identify which internal attitude, magnetic-reference or covariance pathway produces that harmonic.

Eliminating this ripple is not equivalent to solving bias estimation. For example, low-wave northbound OU-III body-Y bias RMS changes from **0.01116 to 0.01444 m/s²** in the waves-plus-yaw control even though its 10 s amplitude falls from **0.002110 to 0.000107 m/s²**. Fixed yaw can further worsen the slowly varying bias error by removing useful excitation. These effects are retained in the raw metrics, not hidden by detrending.

![OU-III southbound controls, 1.5 m incident waves](../reports/results/heading_kinematics/medium-head-wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.svg)

## Runtime separation

`MarineTurnSeparator` is shared by the OU-II, OU-III, independent CARMA,
shared-mode CARMA, TFG and tangent wrappers. It estimates the speed of steady
forward travel from gyro and accelerometer history; it does not read a speed,
track, steering period or simulated attitude.

For measured specific force f and angular rate omega, the transported derivative
`df/dt + omega cross f` cancels the gravity term. In the level bow frame, steady
forward speed U contributes jerk `U * (-r*r, dr/dt)`, where r is heading rate.
The observer regresses this vector against the transported sensor derivative.
Both sides pass through the same low-pass response; the wave measurement itself
is not notched. Centering rejects a constant nuisance offset. Excitation,
coherence, speed bounds and a settling interval prevent a constant heading or
noise-only regression from supplying a correction.

The learned turning acceleration is evaluated in the main estimator's gravity
plane and included in its physical accelerometer prediction, separately from
the wave acceleration. `MarineTurnModel` supplies the attitude and gyro-bias
Jacobians in the native QMEKF or TFG error coordinates. It uses the endpoint
gyro sample and the current gyro-bias estimate; attitude prediction retains
its existing interval integration. The current sensor sample can train only
the following sample's nuisance-speed estimate. The turning term has no
vertical component in that navigation frame. The raw measurement-only
startup, magnetic-reference and sea-state tuning paths retain their sensor
stream. Bias remains a sensor-body state; this correction is not a relabeling or
smoothing of the reported bias.

`get_world_accel()` is the wave-state acceleration. The new
`get_world_total_accel()` adds the estimated steering acceleration and is used
for the replay's existing total-acceleration columns and RMS. Position and
velocity remain wave-relative; every delivered IMU sample and every reference
column retain their original values.

`Config::turn_separator.forward_body` declares the bow in the calibrated body
frame (X by default). The replay explicitly uses Y because its body conversion
is `(x,y,z) -> (y,x,-z)`. This is a mounting convention, not navigation truth.
`Config::turn_separator.enabled = false` disables the model for applications
that do not satisfy its forward-travel assumption.

The IMU-only model assumes approximately steady forward travel with changing
heading and wave jerk that can be separated statistically from that turn.
Arbitrary sideslip, rapid speed changes, and wave motion coherent with steering
need independent navigation aiding; yaw alone does not make those motions
identifiable. Synthetic tests cover different turn periods and speeds, reverse
travel, mounting rotation, zero-speed yaw, constant heading and preservation of
the navigation vertical.

Finite differences check the shared force derivatives and the TFG group
retraction, including the gyro-bias columns. The TFG-CARMA sparse update is
also compared with a dense Joseph update using an independently differentiated
turn model. A zero speed supplies exactly zero force and derivatives.

This is a conditional motion model, not a joint navigation-speed covariance
state. Uncertainty and temporal correlation in the learned speed are not
claimed to be fully represented by the native Kalman covariance.

Full replay validation of the runtime change is in progress. Existing quality
limits are retained; the attribution tables above describe the original runtime
and are not acceptance results for this implementation.
