# OU-III proof: SLOW + FAST physical IMU qualification

Base main: `eef30627130f434eb14ea9f42f4df140a291d341` (2026-10-01).
Controlling derivation: [two-timescale IMU model](ou3-imu-two-timescale.md).
The previous ledger is preserved **verbatim** in
[the pre-migration history](ou3-proof-research-state-before-slow-fast.md).
Its suggested LF cutoff, 1-degree/60-second budget and stronger excitation
are research proposals, not established device qualifications or premises.

## Current hypothesis

Keep the single path construction -> capture -> magnetically informed H18 ->
refinement/release -> recurring A21 -> regional practical stability. Keep the
literal coupled (tau, sigma_aw, R_S, T_S), prediction/due-S/acc chronology,
full-covariance dissipativity and same-history LaSalle structure. No shipping
code, calibration behavior, MARINE MOTION, MAGNETIC SERVICE or gate is changed.

Physical calibrated errors are e_a=b_a,s+b_a,f and e_g=b_g,s+b_g,f.
Slow components have norm/rate bounds and linked increments min(2B_s,D_s*h).
Fast components have amplitude bounds plus signed accumulation on EVERY placed
window 0<T<=H: |integral b_f^hold|<=min(B_f*T,C), C<B_f*H. This is a simple
conditional qualification interface, not an empirical claim of cancellation.
Both H,C pairs are OPEN. The six inherited numerical amplitude/rate values are
explicit candidate deployment budgets; calibration-fit statistics, process
noise, typical sensor densities and estimator projections do not prove them.
One decomposition must satisfy the complete history, including all transitions.

## Evidence and derived inequalities

`imu_temporal.py` implements conditional same-history calculus. Distinct raw
sample epochs on complete hold cells have the global error-class supremum

    Delta_i(h)=min(2B_i,s,D_i,s*h)+c_i(dt_0)+c_i(dt_1),
    c_i(dt)=min(B_i,f,C_i/min(H_i,dt)).

The raw fast term can still be 2B_f. Averaging/cancellation never justifies a
smaller raw-point bound without the cell qualification. With K*(T) the tiled
window cap, shifted boxcar differences have fast bound
2 min(K*(L),K*(h))/L. Matrix Abel identities retain signed gain/frame variation,
the endpoint primitive and division by sample duration. They are support
bounds on one history, not a reason to independently maximize word loss and
forcing. Unknown temporal profiles fail closed in the physical validator and
in the sensor-supply helper. A finite prefix never certifies an all-time device.

The exact linked completed square remains valid:

    G=J0-M' JN M-gamma J0>0, z=M' JN b,
    chi=b' JN b+z' G^-1 z,
    V_N=(1-gamma)V_0+chi-(e_0-G^-1 z)'G(e_0-G^-1 z).

The required supremum is over the SAME reachable slow/fast histories that
produce M, b, covariance, tuner and scheduler. Physical BA/BG coordinates now
refer to the slow components; fast errors stay in the sensor forcing. The
H18 and active-A21 bias prediction identities keep their literal factors.
Every-prefix retention and chi/gamma<r_in^2 remain unproved.

## Counterexamples: classify by model, do not erase

The old V3 norm-only witness phi=.01 sin(t/2), with arbitrary n_a/n_g and a
constant BA offset, remains valid under its explicitly frozen OLD model.
Its all-slow realization requires .04903325 m/s^3 and .0025 rad/s^2: respectively
49.03325 and 250 times the inherited candidate slow-rate budgets. This excludes
only the all-slow realization. For ANY split, opposing half-windows of length
2*pi require continuous K_a>=.3725224 m/s and K_g>=.0198026 rad. Charging the
qualified 6-ms sampling defect gives conservative necessary requirements
K_a>=.3706550 m/s and K_g>=.0197081 rad. These are witness requirements, NOT
chosen sensor specifications. Missing H,C means current-model admission and
exclusion are both OPEN. The old sqrt(V)>=.4 bound cannot be transferred by
changing which part of the error is called physical slow BA.

The former phi=.001 sin(t/40)^3 construction has p=v=a=0 and therefore does NOT satisfy the amended MOVING displacement-span condition. It is retained only as a counterexample to the previous attitude-only MARINE contract. Charge the tiny gravity-representation
constant to slow accelerometer bias. Its slow amplitude/rate bounds obey the
inherited candidates, and the quiet execution retains MAGNETIC SERVICE.
It fits symbolic MARINE T_E=80*pi, theta_E=.002 on complete moving windows.
Thus that zero-translation construction no longer blocks the amended MARINE contract. Full moving compatibility must now include the nonzero displacement span together with attitude, SLOW+FAST IMU and magnetic service. This does
not refute every fixed deployment pair (T_E,theta_E), nor practical stability
relative to the compatible class. No EXCITED_MOVING strengthening is adopted.

## Retained facts and failed approaches

Retain the operation energy identities, linked completed square, signed
variation of constants, BA marginal elimination, source covariance comparisons
within their stated profile, and conditional radius-local field-axis theorem.
The previous .15 local result still requires its ACTUAL local-error premise.
Do not upgrade it to outer entry from physical tilt span alone.

**DEAD ENDS:** unrestricted residual boxes labelled fast; deriving deterministic
all-history caps from RMS/Allan/typical noise alone; fitting a cutoff to reject a
witness; deleting gyro or frame-weight variation; using only aligned-window
means; resetting the split at proof boundaries; inferring zero base innovations
from zero homogeneous action; S-only LIN contraction substituted for actual
interleaved H18; lower covariance bound substituted for full compactness; old
O1/O2 kernel-ceiling architecture used as the controlling route. Archived algebra
may remain correct as algebra or an explicitly amplitude-relaxed outer bound.

## Current limiter and alternatives

Numerical assembled-device slow/fast qualification is OPEN for both sensors.
The complete physical joint acc/gyro/magnetic reachable ambiguity must then be
compared with the EXISTING fixed MARINE constants, not a bias-only threshold.
H18 release compactness/BIBO for the literal interleaved word, outer retention,
physical-to-nominal transfer, linked finite supply, nonlinear/prefix retention,
regime composition and float32 totality also remain OPEN. All end-to-end theorem
flags remain false. No finite carried replay supplies the missing uniformity.

An observable-consistency-class practical target remains a legitimate option,
but is not substituted for the current theorem without stating the difference.
No behavior change or automatic STILL detector is justified here.

## Next falsifiable calculation

Independently qualify both delivered-stream fast horizon/cap pairs and the
slow budgets with residual/reference evidence. For a declared (T_E,theta_E),
intersect BOTH exact sensor compatibility equations, their transported gyro
integral, the physical kinematic continuation and actual MAGNETIC SERVICE.
Evaluate the linked chi and every-prefix bounds on that same reachable class.
Do not spend more enclosure effort on the refuted old norm-only .15 claim or
try to close a source theorem by multiplying independent norm suprema.

## Validation scope

See `imu-two-timescale-certificate.json` for exact rational scalar consequences
and `test_ou3_imu_temporal.py` for finite-window/cross-boundary algebra tests.
Finite-prefix and synthetic-test success is not device qualification, all-time
float32 validation, a uniform contraction certificate or theorem completion.
