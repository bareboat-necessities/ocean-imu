# Rapid-motion failure: equations and causal controls

The shared symptom does not establish one shared implementation error. The
controls distinguish a qMEKF covariance implementation error, an invalid
gravity measurement assumption during acceleration, local antipodal correction
equilibria, and magnetic timing error. Hardware failure is not needed for any
of these reproductions. A reproduction is not identification of the cause of an
unrecorded hardware turn.

## qMEKF covariance used the wrong rotation

`QuaternionMEKF` stores a body-to-world quaternion and injects a right error:
`q_true = q_est Exp(delta_theta)`. For constant body rate `omega`, its error
dynamics are

\[
\dot{\delta\theta}=-[\omega]_\times\delta\theta-\delta b_g,
\qquad \dot{\delta b_g}=0.
\]

Writing `phi = omega dt`, the correct discrete transition is

\[
\Phi=\begin{bmatrix}
\exp(-[\phi]_\times)&-dt J_r(\phi)\\0&I
\end{bmatrix},\qquad
J_r(\phi)=I-\frac{1-\cos\theta}{\theta^2}[\phi]_\times
+\frac{\theta-\sin\theta}{\theta^3}[\phi]_\times^2.
\]

The old code instead used the upper-left 3x3 block of its **4x4 quaternion**
propagator. That block is `cos(theta/2) I - sin(theta/2) [u]x`. Its spin-axis
eigenvalue is `cos(theta/2)`, rather than one. In the bias-free layout, with
`P=I`, zero process noise, no observations, 720 deg/s and 5 ms samples,
five seconds of prediction gives
axial variance **0.3726473193**, rather than one. A deterministic rotation
cannot create that information. Transverse errors also rotate by the wrong
angle; the bias cross block must integrate that rotation rather than use
`-dt I`.

There was a second error: quaternion injection changed the error chart without
transporting covariance. After injecting `d`, the new error is
`Log(Exp(-d) Exp(d+epsilon))`; its derivative is `J_r(d)`. The implementation now
transports both attitude covariance and its bias cross terms with this reset.
The existing per-sample process-noise convention is unchanged.

`tests/ahrs/qmekf_geometry-test.cpp` differentiates the physical quaternion flow
and the finite reset independently of these matrix formulas. It covers rotated
initial attitudes, anisotropic covariance, both bias layouts, multiple spin
axes, small rates and 75 ms steps. Both prediction and reset tests reject the
old header and pass the corrected implementation.

These corrections **do not alone eliminate the offset-rotation flip**. With all
ordinary observations enabled, the 720 deg/s, 15 cm yaw case still reaches
179.86 degrees of total attitude error. The combined-axis case improves from
179.99 to 75.35 degrees, but that is not a satisfactory attitude estimate.

## The first false observation is specific force treated as gravity

For a sensor fixed at body offset `r` from an origin with world acceleration
`a_o`, the measurement equation is

\[
f_b=R^T(a_o-g)+\alpha\times r+
\omega\times(\omega\times r)+b_a+n_a.
\]

The gravity-only compass observation instead uses `g f_b / ||f_b||` as an
observation of `-R^T g`. Normalization does not remove either rotational term.
At 720 deg/s and 15 cm, the centripetal term is **23.687 m/s² (2.415 g)**.
For level yaw motion, treating that vector as gravity asks the observer for
**67.51 degrees of false tilt**. The gyro prediction can be exact while this
measurement correction moves the attitude away from truth.

False tilt can reverse north. For a world field `(H,0,Z)` and a pitch error
`-epsilon`, the computed north component is

\[
\hat m_N=H\cos\epsilon-Z\sin\epsilon.
\]

It changes sign at `epsilon = atan(H/Z)`. For the replay's `(20.3,0,46.8)` uT
field that is **23.45 degrees**. The unprotected qMEKF yaw replay exceeds
50 degrees of tilt error. Its heading loss is therefore not evidence that the
physical magnetic field reversed.

Matched qMEKF controls keep the same rate, specific-force history, noise
settings and correction schedule. Replacing the adaptive dip with the correct
fixed world field still reaches 179.68 degrees. Removing the gyro-bias state
still reaches 179.76 degrees, but ends below 0.001 degrees after the turn:
bias learning explains the persistent tail, not the initial flip. Removing the
known rotational force term from the measurement equation keeps the same turn
below 0.001 degrees with **every ordinary correction enabled**.

That last experiment is a known-truth causal control. It does not estimate an
unknown hand-motion pivot, and must not be presented as a deployed solution
that knows the user's acceleration. OU-II/III have an explicit configured
lever-arm model; their zero-offset configuration does not model an unknown
rotation pivot. Their acceleration states and OU/pseudo measurements, and the
TFG, NLO and PII acceleration assumptions, differ from the gravity-only compass.
They must not be described as sharing qMEKF's covariance bug.

The same noiseless 720 deg/s, 15 cm yaw experiment gives the following total
attitude errors in degrees. Every ordinary correction remains enabled in both
columns; neither protection nor reacquisition acts. Tail is seven seconds after
the turn. The force control uses the known physical term above.

| Backend | Measured force: peak / tail | Known rotational force removed: peak / tail |
| --- | ---: | ---: |
| qMEKF | 179.85678 / 13.77477 | 0.00005 / 0.00000 |
| Mahony | 179.98695 / 8.86270 | 2.92779 / 0.08002 |
| OU-II | 15.56995 / 4.54497 | 0.00040 / 0.00011 |
| OU-III | 30.74688 / 8.82418 | 0.00041 / 0.00012 |
| TFG | 54.62368 / 16.81752 | 0.00026 / 0.00016 |
| NLO | 49.31396 / 48.29376 | 3.60251 / 0.00012 |
| PII | 5.08808 / 2.92147 | 0.81299 / 0.42065 |

The remaining centered-control errors in the continuous observers show why
this experiment is not a claim of exact discrete propagation for every backend.
The complete rate, radius, combined-axis, noise, gap and latency controls are in
the [336-row result](../reports/results/rapid_motion/causal-controls.csv).

## Why a local magnetic correction can remain at 180 degrees

For predicted and measured unit horizontal north vectors separated by `e`,
their cross product is proportional to `sin(e)`. It is zero at both zero and
180 degrees. In a local Kalman vector update, at antiparallel vectors the
nonzero Cartesian residual is radial, while the attitude Jacobian spans the
tangent plane. With the symmetric quiet covariance in this test, its Kalman
attitude correction is also zero. Reducing noise or repeating the observation
does not supply a missing correction direction.

`antipodalInnovation()` uses an exactly representable wrong quaternion
`(w,x,y,z)=(0,0,0,1)`, quiet gravity and horizontal magnetic field. It contrasts
ordinary local updates with the shared finite-angle capture using identical
measurements. The test avoids `sin(pi)` roundoff, which would seed an eventual
escape and obscure the equilibrium. The existing tilted, inclined-field and
noisy near-180 tests separately exercise practical recovery.

All seven local controls retain exactly 180 degrees after three seconds; all
seven protected controls apply one capture and finish at zero degrees to the
printed precision. This establishes the equilibrium for the specified field and
covariance, not for every possible filter state or magnetic inclination.

The recovery operation solves attitude from two qualified vectors rather than
extrapolating the local cross-product innovation through its zero. This is a
separate correction to the capture behavior, not a claim that suspended
updates become informative.

## Timing is an independent forcing term

Using `m(t-delay)` as though it were `m(t)` violates the measurement equation
even at zero sensor offset. For constant rate the angular mismatch is
`omega * delay`: 720 deg/s and 40 ms gives 28.8 degrees before observer
feedback. The combined-axis, centered controls vary only magnetic delay and
also develop large heading errors. A host read timestamp does not identify the
BMM150 conversion epoch or each axis's measurement aperture.

## Inherited-state recovery is a separate, incomplete qualification

`--diagnose-recovery` runs ordinary aiding through the same noiseless
720 deg/s, 15 cm turn, then compares two continuations from identical histories:
ordinary aiding throughout, or protection enabled after the turn at 11 seconds.
Neither arm reconstructs the filter or injects an artificial heading error.
Bias, covariance, navigation, magnetic reference, delivery and guard histories
all continue. The CSV records the inherited attitude and bias errors, the tail
seven seconds later, actual captures and the final hold state.

The corrected qMEKF still fails the five-degree diagnostic recovery target:
its yaw-turn tail is 13.77477 degrees with ordinary aiding and 16.54773 degrees
with late protection. The latter performs one capture. Both inherit a
1.821243 rad/s gyro-bias estimate from the turn; the late-protection tail still
has 1.158784 rad/s of estimated bias although the physical gyro bias is zero.
The combined-axis late-protection tail is 17.24282 degrees. Preserving a
previously corrupted bias estimate is therefore not the same qualification as
preserving a healthy learned bias during attitude capture.

TFG's yaw-turn tail is also above target: 13.96484 degrees with late protection,
versus 16.81752 degrees with ordinary aiding. All fourteen late-protection cases
apply one capture; three remain above five degrees (both qMEKF cases and TFG
yaw). The other eleven tails meet this diagnostic target. The original 336
causal controls and all protected-regression metrics reproduce unchanged.

This deliberately disables prevention during the turn; it does **not** refute
the always-protected turn regressions or identify an actual device failure.
It does refute interpreting those regressions and artificial-angle capture
tests as sufficient evidence of recovery from every failed-turn state.
The diagnostic reports failures without converting them into passing acceptance
assertions. The source-bound results are in
[`recovery-controls.csv`](../reports/results/rapid_motion/recovery-controls.csv).
Recovery of corrupted carried state remains open; no runtime bias reset is
justified by these controls alone.

## Reproduction and scope

Run `python3 tools/rapid_motion_diagnosis.py --eigen <Eigen include directory>`.
The driver freshly builds each backend, runs 24 ordinary-aiding controls and
their 24 known-force-equation counterparts, four inherited-state recovery
controls, and the protection/recovery regressions. Source and executable hashes accompany
`reports/results/rapid_motion/causal-controls.csv`. Recovery counts record
operations actually applied, not decisions calculated by a disabled guard.

These finite synthetic controls preserve each production core's prediction,
measurement, covariance and bias operations through thin delivery adapters.
They do not exercise the entire device startup/calibration/tuner history, do
not qualify magnetic timestamps, and are not admitted OU-III theorem histories.
OU-III ordinary-aiding identities retain their original scope; no stability or
hardware qualification is inferred from these experiments.
