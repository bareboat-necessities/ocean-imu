## Candidate FAST gyro removes the historical continuous field-axis obstruction

The old continuous obstruction used n_g=.02 sin(t)b and required about .04 rad signed accumulation over a half-cycle. The candidate C_g=.002 rad excludes it by a factor 20. In the exact 17-s field-axis Stieltjes telescope, the inherited slow charge is 2 P_max B_g,s/T + P_max D_g,s = .0129961764705882 m/s^2 and the candidate FAST primitive contributes at most P_max C_g/T = .000647058823529412 m/s^2, total .0136432352941176 m/s^2 before coefficient-rotation and exact SO(3) correction-jump/twist terms. Thus the historical .055 m/s^2 witness is no longer admissible. The controlling correlated remainder is now the source-uniform sum of coefficient rotation plus exact correction jump/twist action; do not fall back to independent gain or total-NIS norms.

## Signed AW reader: uniformity status

The carried six-profile stress family now passes the recovered AWB2 allowance: worst signed 16-s error .371143 m/s^2 versus .89553839501604595, leaving .52439539501604595 m/s^2. Compactness/continuity of the qualified regular 16-s same-history word class proves existence of a finite source-uniform B_AW=sup|Phi|. Strict numerical closure is reduced to one dependency-preserving finite-cover inequality L_Phi delta_F + E_nl + E_f32 < .52439539501604595. L_Phi and delta_F are not yet certified; independent gain/tuner/covariance boxes are invalid because they destroy the paired cancellation. This is now the controlling quantitative AW-transfer obligation.

## Recovered non-recurrent-acceleration bridge

The earlier useful mechanism is recovered explicitly. It does not assume recurring acceleration. Bounded physical velocity gives the signed physical mean charge 2 V_max/L, jerk/sample fidelity adds J_max h_acc/4, and Corollary A* needs only a nominal signed AW mean transverse to the magnetic field. At L=16 s the physical charge is 0.8375 m/s^2. Under the candidate H_a=60 s,C_a=.05 m/s profile, the FAST signed mean over 16 s is <=.05/16=.003125 m/s^2. Even conservatively charging the entire slow accel amplitude .22516660498395405 leaves 1.96133-.8375-.003125-.22516660498395405 = .89553839501604595 m/s^2 of field-separation margin at the sensor-history level. This is source-uniform and uses no recurrent-acceleration premise. The remaining bridge is specifically the literal time-varying acc/S gain loop: sensor signed mean does not automatically bound nominal AW signed mean because gain rectification is possible. The exact AW-loop/two-Abel identity must transfer this .895538395 margin without replacing the gain-weighted history by independent boxes.

## Candidate assembled-device FAST profile

Use H_a=H_g=60 s, C_a=.05 m/s, C_g=.002 rad as the current candidate qualification. On any 30-s MARINE window the signed FAST accumulation caps remain .05 m/s and .002 rad. The known phi=.01 sin(t/2) witness needs at least .370655073155 m/s and .0197081716368 rad after sampling charge, so this candidate excludes that witness by either sensor separately, with reserves .320655073155 m/s and .0177081716368 rad. This is not yet the general same-history gauge theorem; hardware validation of these candidate caps also remains required.

## User-qualified MARINE excitation constants (2026-10-01)\n\nFor the controlling theorem, use T_E=T_P=30 s, theta_E=2 deg, P_E=.03 m. These are theorem qualification premises supplied by the user; they are not inferred from RAO data. IMU fast H_a,C_a,H_g,C_g remain OPEN. Consequences: P_E/T_P=.001 m/s; theta_E/T_E=.00116355283466 rad/s. The generic signed acceleration chord lower max(0,P_E-T_P V_max) remains zero, so no pointwise acceleration floor follows.\n\n# OU-III proof: SLOW + FAST physical IMU qualification

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

## Exact joint gauge qualification

The controlling quantifiers are now written explicitly in [joint SLOW+FAST physical gauge separation](ou3-joint-gauge-separation.md). For fixed qualified MARINE tuple m=(T_E,theta_E,T_P,P_E), Q_IMU(m) is the set of fast horizon/cap tuples for which the infimum of the SAME-history stacked acc+gyro causal compatibility residual, restricted by actual MAGNETIC SERVICE, is strictly positive. Process/accelerometer columns are paired and physical p/v/S/a boundaries telescope before norms. Because the four MARINE excitation constants are also numerically OPEN on main, the current empirical qualification is Q_phys over both MARINE and IMU parameters; no four-dimensional numerical Q_IMU is asserted before m is fixed.

## Current limiter and alternatives

Numerical assembled-device slow/fast qualification is OPEN for both sensors.
The complete physical joint acc/gyro/magnetic reachable ambiguity must then be
compared with the EXISTING fixed MARINE constants, not a bias-only threshold.
Held-H18 LIN BIBO and captured-domain release compactness are CLOSED qualitatively/source-uniformly by the post-#640 certificate; outer retention,
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

## CI: information-ratio replay comparison

Failed quantity: `information_ratio_source_diagnostic --expect` on the
`theorem` job compared `diameter.relative_difference` (committed 1.03e-5 max,
fresh 2.3e-6 on the same case) with a 5e-3 relative tolerance. Classification:
implementation/CI. That leaf is the residual of an identity between two equal
quantities, so its value is float32 replay noise and is not reproducible
across toolchains. Invalidated hypothesis: every recorded leaf is well posed
for relative comparison. Retained facts: the identity contract
`relative_difference < 1e-3` is still enforced by `verify_diagnostic`; the
replay comparison now only requires both records below 1e-4; every other leaf
keeps the 5e-3 relative comparison. Current limiter and next falsifiable
experiment are unchanged from the sections above.

## CI: held-BA LIN word diagnostic did not execute

Failed quantity: `held_ba_lin_word_diagnostic` on the `theorem` job raised
`JSONDecodeError` (`"root_covariance":,`). Classification: implementation.
Its driver patch removed the root-covariance capture, compiled against the
untapped shipping header (zero recorded events), and recorded past first BA
activation, so the "last 3400 predictions" would not have been the held
window. Invalidated hypothesis: none; the diagnostic had not produced a
number. Fix: capture the root at the first Live sample, compile with the
`ag_readout_source_diagnostic.instrument` tapped header copy, and stop
recording permanently at the first `acc_bias_updates_enabled` transition.
Retained facts (finite carried evidence only): quiet live 18051 / active
24064, rho_LIN=0.03953, sigma_max=0.3258; wave live 6368 / active 24016,
rho_LIN=0.006467, sigma_max=0.03199. This does not establish source-uniform
held-BA LIN BIBO stability; the next falsifiable experiment above stands.
