# IMU BIAS: one carried SLOW + FAST model

Controlling proof-side correction to main `eef30627130f434eb14ea9f42f4df140a291d341`.
The estimator, calibration, frontend, tuner, scheduler, gates, MARINE MOTION and
MAGNETIC SERVICE are unchanged. This replaces the V3 monolithic bias plus
unrestricted residual interpretation. It does **not** assert new device evidence.

## SF1. Compact physical contract

At the actual calibrated delivered epochs, in the same body coordinates,

    y_a = f_true + b_a_s + b_a_f,
    y_g = omega_true + b_g_s + b_g_f.

For i=a,g, the slow component is locally absolutely continuous with

    |b_i_s| <= B_i_s,   |dot b_i_s| <= D_i_s,
    |b_i_s(t+h)-b_i_s(t)| <= min(2 B_i_s,D_i_s h).

The fast component has |b_i_f|<=B_i_f. Give it the following *conditional
qualification*, not a noise-color assumption. Let f_i be the timestamped
zero-order hold of its **delivered calibrated samples**, and U_i'=f_i. For
every placement s and every 0<T<=H_i,

    |U_i(s+T)-U_i(s)| <= K_i(T) := min(B_i_f T,C_i),
    0<H_i<infinity,   0<=C_i<B_i_f H_i.                 (SF1)

If B_i_f=0, use C_i=0. Units of C_a are m/s; units of C_g are radians.
No derivative is imposed on fast errors. This is a bounded signed accumulation
condition, **not** independence, white noise, zero mean, or a correlation time.
It permits persistent small leakage up to C_i/H_i; it does not permit arbitrary
persistent error at the whole instantaneous amplitude B_i_f. That distinction
is intentional. Neither a bounded-primitive assumption over infinite horizons
nor exact DC rejection is invented.

**H_a,C_a,H_g,C_g are all OPEN (null in constants.json).** The repository does
not yet support numerical values for this interface. Missing values mean that
physical qualification cannot be completed, NOT that the fast channel becomes
unrestricted. The definition is useful before qualification because it states
exactly which temporal facts each analytical estimate needs. No claim is made
that all assembled BMI270 devices satisfy a chosen, as-yet-unspecified profile.

One decomposition is selected for the *entire* history. Neither proof-word
boundaries nor hold/release/physical regime changes may reassign error between
slow and fast components, reset U, or choose a new bias origin. U itself is only
a proof bookkeeping primitive: its differences, not an absolute origin, matter.
Calibration scale/cross-axis/model residuals must fit these same two channels;
if their motion-correlated component cannot do so, qualification fails. There
is no third unqualified residual channel.

The hold is not a runtime change or an assumption of analog white noise.
Cancellation before hardware filtering, downsampling or calibration does not
qualify the delivered signal. Continuous physical motion remains continuous;
converting a continuous integral to delivered samples requires its existing
sampling/jerk/rotation error charge.

## SF2. What the existing constants actually establish

The six retained numbers are **declared candidate qualification budgets**.
They are not six newly measured device limits, and a bound on an old total
error does not prove separate bounds on an arbitrary slow/fast decomposition.

| New candidate | Value | Audited interpretation and limitation |
|---|---:|---|
| B_a_s | 0.22516660498395405 m/s² | Old sqrt(3)*0.13. Old prose says commissioned per-axis residual; no linked assembled capture establishing this slow envelope was found. Not a calibration covariance sigma. |
| D_a_s | 0.001 m/s³ | Old selected deployment limit. No assembled all-time rate evidence. |
| B_g_s | 0.02 rad/s | Inherited residual norm budget. No linked assembled slow-component certificate. |
| D_g_s | 0.00001 rad/s² | Inherited rate budget, not a process-noise setting or measured rate specification. |
| B_a_f | 0.3 m/s² | Old instantaneous proof residual envelope. Not an RMS/PSD-to-pathwise guarantee. |
| B_g_f | 0.02 rad/s | Old instantaneous proof residual envelope. Not the gyro noise density or sigma. |

Provenance audit:

* `src/imu_calibrate/CalibrateIMU.h`, `AccelCalFit.h`, `AccelCalCapture.h`,
  `GyroCalCapture.h`, and `src/AtomS3R/AtomS3R_ImuCalBlob.h` implement finite
  calibration/capture checks and the actually applied affine correction.
  Accelerometer norm residuals, held-out fit errors, cross-validation and
  parameter standard deviations are not all-time vector-error envelopes.
  For example, accel block/model floors .002/.003 m/s² and maximum per-axis
  bias sigma .015 m/s² are fit/information settings, not B_a_s,D_a_s or C_a.
  Gyro block scatter and the .0002 rad/s bin-level noise floor likewise do not
  establish a future residual-accumulation bound. The code specifically avoids
  shrinking correlated bin uncertainty merely by counting raw samples.
* `tests/imu_calibrate/{accel_cal-test,imu_calibrate-test,calibration_accuracy-test,
  calibration_safety-test,calibration_workflow-test}.cpp` and the simulation/
  replay helpers exercise known synthetic truth, errors, motion, fit rejection,
  persistence and finite captures. Synthetic white noise is a fixture input,
  not a theorem premise about a real assembled sensor. No all-time two-timescale
  qualification follows from passing these tests.
* The filter's sigma_a=.2 m/s², gyro noise densities, bias process covariance,
  bias OU time constant, estimate projection radii .4/.5 and all coupled tuner
  constants are **estimator/tuning** quantities. They are not physical bias
  amplitudes, physical drift rates or temporal fast-error bounds.
* Bosch's official BMI270 product specification lists component offsets,
  sensitivities and noise densities (including accel 160 micro-g/sqrt(Hz),
  gyro .007 degrees/s/sqrt(Hz)). These do not specify a deterministic
  accumulated residual bound after assembled-device calibration and sampling.
  Source: https://www.bosch-sensortec.com/en/products/motion-sensors/imus/bmi270
  (consulted 2026-10-01). No RMS/PSD conversion is used to fill the open cells.

The existing .3 m/s² **vibration detector-band RMS** envelope is retained as a
separate qualification of the actual frontend/Racc inflation. It is not implied
by the new fast amplitude: slow leakage and the actual filter response also
matter. It is not an additional physical bias category. Numerical timing,
Racc, noise and gate settings are unchanged.

## SF3. The actual two-epoch envelope, not independent endpoint noise

Define S_i as the complete slow-history set and F_i by SF1 on the same delivered
timestamp grid. The error-history set is E_i=S_i+F_i; its restrictions to epochs
are **projections of complete histories**, not products of per-epoch balls.
A fixed predecessor imposes additional intersections with already accrued
window constraints.

For two distinct complete hold cells of lengths d_0,d_1 at epochs separated by
h>=d_0>0, SF1 implies the exact single-cell cap

    c_i(d)=min(B_i_f, C_i/min(H_i,d)).

Hence the global error-class two-epoch supremum is

    Delta_i_s(h)=min(2B_i_s,D_i_s h),
    Delta_i_f(h;d_0,d_1)=c_i(d_0)+c_i(d_1),
    Delta_i(h;d_0,d_1)=Delta_i_s(h)+Delta_i_f(h;d_0,d_1). (SF2)

Proof: a constant cell contains windows of all lengths up to min(H_i,d).
Conversely, put opposite collinear fast values of magnitudes c_i(d_0),c_i(d_1)
on the two cells, zero elsewhere. Any admissible window includes at most the
positive or negative cell's capped integral in absolute value; opposite signs
cannot increase it. This realizes the fast supremum. A clipped slow ramp can
align its change with that difference and realize the slow supremum. Previously
fixed predecessors, physical measurements and other equations can only shrink
this unconstrained global class supremum. At h=0 the exact difference is zero.

In particular, when C_i>=B_i_f max(d_0,d_1), **Delta_i_f=2B_i_f**. A temporal
accumulation bound need not shrink raw sample differences: two brief opposite
peaks are possible. When H,C are unknown, SF2 is symbolic/OPEN; the older
min(2Bs,Dsh)+2Bf remains only an amplitude outer bound, never an admission test.
An endpoint whose next timestamp is absent has no certified cell length; do not
apply the complete-cell improvement to it.

## SF4. Averaging is where temporal qualification changes the comparison

For T=nH+r, n>=0, 0<=r<H, a valid long-interval envelope is

    K_i*(T)=n C_i+min(B_i_f r,C_i),   K_i*(0)=0.          (SF3)

This follows by tiling the *same* history, not resetting it. Actual reachable
integrals may be smaller; preserve the whole window set for sharper estimates.
For a boxcar of length L,

    |average_L f_i| <= K_i*(L)/L,
    |average_L f_i(t+h)-average_L f_i(t)|
        <= 2 min(K_i*(L),K_i*(h))/L.                   (SF4)

The second bound follows both by bounding the two L-windows and by cancelling
their overlap, leaving two h-windows. For continuous slow averages add
min(2Bs,Dsh). For averages of held slow samples on a mesh with maximum gap d_max,
a safe term is min(2Bs,Ds(h+d_max)) for h>0; alternatively charge the explicit
quadrature defect in each window. Do not differentiate held values as though
they were the continuously differentiable physical slow bias.

## SF5. Signed, rotating and actual-gain functionals

Every source functional uses its actual weights/frames on the same history.
For complete sample cells d_k, let A_k include all quadrature weights and
transported gains, q_k=A_k/d_k and U_{k+1}-U_k=d_k b_f,k. Then exactly

    sum_(k=m)^(n-1) A_k b_f,k
      = q_(n-1)(U_n-U_m)
        + sum_(k=m+1)^(n-1)(q_(k-1)-q_k)(U_k-U_m).    (SF5)

Its norm is bounded by the endpoint coefficient times K*(t_n-t_m), plus
sum |q_(k-1)-q_k| K*(t_k-t_m), or by the cell amplitude bound if smaller.
**Retain the endpoint term.** Dividing by the actual d_k before differencing,
and retaining rotating matrix coefficients, is essential. A signed unweighted
mean cannot bound arbitrary coning, switched-gain or adjoint-weighted effects.

For slow samples, w_j=b_s,j-b_s,j-1 gives exactly

    sum A_k b_s,k=(sum A_k)b_s,0
                   + sum_(j>=1)(sum_(k>=j)A_k)w_j,    (SF6)

with |w_j|<=min(2Bs,Ds d_j) and every partial sum inside the same Bs ball.
The exact functional support is the supremum over those linked constraints
and SF1. Bounding separate summands is permitted as a conservative inequality,
but choosing their extremes independently does not produce a physical witness.
`imu_temporal.py` implements SF2--SF6 and all-placed-window finite-prefix tests.
These float computations are diagnostics of analytical formulas, not validated
interval certificates or all-time device certificates.

## SF6. Joint accelerometer/gyro/magnetic ambiguity

For two physical histories producing the same delivered packets,

    R_1' (a_1-g_1) - R_0' (a_0-g_0) = e_a,0-e_a,1,
    omega_1-omega_0 = e_g,0-e_g,1.                     (SF7)

Thus the right sides belong to **E_a-E_a and E_g-E_g**, respectively, with both
histories' slow and fast dynamics retained. For an unspecified pair of histories
the two-epoch outer bound is twice SF2, not SF2. A fixed zero-error reference
allows a single-history charge. These sensor conditions must hold jointly with
one physical R,a,v,p,jerk,primitive continuation and the *unchanged actually
applied* MAGNETIC SERVICE rows. A comparison of gravity span alone that ignores
translation or gyro compatibility is not the desired theorem.

With Q=R_1 R_0', the exact relative-rotation derivative is

    dot Q = R_1 [omega_1-omega_0]_x R_0'.              (SF8)

So gyro accumulation enters a transported matrix functional, not an unjustified
bound by |integral b_g_f| alone. For a fixed-axis witness it reduces to a scalar
integral. A constant physical field leaves rotations about that field invisible
to the magnetometer; actual heading/axial-bias service alone does not exclude
this physical ambiguity. No zero *base* innovations follow from zero homogeneous
action. The compatibility line still needs its same-history physical-to-nominal
bridge in the LaSalle argument.

**Existing MARINE first:** its numerical T_E and theta_E are still OPEN. The
all-slow family phi=.001 sin^3(t/40), p=v=a=S=0, B=75 e_x, has gravity span .002
on every T_E=80 pi window, b_a_s=g(R'ez-ez) plus the constant literal-gravity
representation offset, and b_g_s=-dot(phi)e_x. Fast components are zero.
Its global slow bounds are

    B_a_s <= .00980665 + 1.61744e-7, D_a_s <= .00073549875,
    B_g_s <= .000075,               D_g_s <= .000005625.

They lie inside the inherited candidate budgets. Its identical quiet packets
retain the existing analytical service lower bound (>1). It therefore survives
**every** proposed nonnegative fast accumulation cap. This is a valid obstruction
to deducing a universal strict gauge-breaking margin from the presently symbolic
MARINE family, not a refutation for every fixed T_E/theta_E, nor of practical
stability of the measurement-compatible class. No stronger EXCITED_MOVING
assumption is added. First qualify the existing constants and evaluate SF7--SF8.

## SF7. Retest phi(t)=.01 sin(.5t), independently of model selection

For zero translation and fixed field along the rotation axis, the compensating
accelerometer change has maximum derivative g*.01*.5=.04903325 m/s³;
the gyro compensation has maximum derivative .01*.5²=.0025 rad/s².
These are respectively **49.03325 and 250 times** the candidate slow-rate
budgets. All-slow compensation is excluded under those candidate rates.
This does NOT exclude all mixed slow/fast decompositions.

A useful *necessary condition for any decomposition*, not a fitted device
constant, follows from opposite half-windows of length L=pi/.5=2 pi. The lateral
accel integral on a positive half-wave obeys

    I_a >= (g/.5)(2*.01 - 2*.01³/9),  I_g=2*.01.

Subtract two opposite windows. A continuous slow term contributes at most
L min(2Bs,Ds L) to their difference; the two fast windows contribute at most
2K*(L). Thus any continuous candidate split needs

    K_i*(L) >= I_i - (L/2) min(2Bs,Ds L).              (SF9)

The resulting lower requirements are about .37252243 m/s (accel) and
.019802608 rad (gyro). For the actual delivered hold, conservatively charge
L D_error d_max in each integral comparison and use Ds(L+d_max) for held slow
samples. The necessary lower requirements become

    K_a*(2 pi) >= .370655073155... m/s,
    K_g*(2 pi) >= .0197081716368... rad.               (SF10)

`imu-two-timescale-certificate.json` gives exact rational conservative bounds,
using a rational enclosure of pi and d_max=.006 s. If a **future independently
qualified** profile has either cap below its SF10 requirement, this witness is
excluded, including every possible slow/fast allocation. If not, these necessary
conditions do not establish admission; the complete reachable sets still matter.
At present both K profiles are unknown: **new-model admissibility is OPEN**.

The former norm-only counterexample is preserved in
`finite_residual_obstruction.py` and the archived ledger, explicitly labelled
V3. Its chosen constant slow bias cancelled by residual noise, and the resulting
sqrt(V)>=.4, cannot automatically be transferred to another decomposition.
The old hard-coded .25-Hz/1-degree/60-s LF proposal is historical exploration,
not an authoritative IMU or motion qualification. Its finite-capture tool no
longer emits an all-time `qualified=true` decision.

## SF8. Switch every stage of the existing proof, not the estimator

The proof bias error is now b_s-b_hat for both sensors. This is a proof-coordinate
interpretation; no estimator state, covariance, code path or tuning changes.
The same physical slow recurrence supplies one w_s to physical truth and error.
For accel prediction e_b^- = phi_e e_b^+ +(1-phi_e)b_a_s+w_a_s, with phi_e=1 in
H18 and the **literal** OU coefficient in A21; gyro uses the identity predictor.
Fast error enters the actual accel correction and gyro prediction operation;
it is neither another estimated bias state nor an OU law on physical truth.
Correction and projection defects stay separate in the full covariance metric.

For one realized word retain e_N=M e_0+b and the linked LS1--LS7 identity in
`ou3-linked-finite-supply.md`. Replace any independently selected physical input
box by the projection of E_a,E_g on that exact word, intersected with MARINE,
MAGNETIC SERVICE and the same generated frontend/tuner/covariance history.
In particular the controlling bound is

    chi_* = sup_{same reachable history}
      [b' J_N b + (M'J_N b)'(J_0-M'J_N M-gamma J_0)^(-1)(M'J_N b)],

not separate suprema for M, b, biases or fast residuals. SF5--SF6 bound its
sensor functionals only with the actual linked coefficients. Nonlinear,
reference, projection, OU mismatch, sampled physical motion and arithmetic
supplies are still present; they must not be relabelled arbitrary IMU fast error.

| Active use | New controlling treatment |
|---|---|
| Attitude/BA and gyro compatibility | SF7--SF8 over E_a,E_g, predecessor-linked slow increments and all fast windows; pair-history factor of two where required. |
| MARINE sampling and excitation | Existing jerk/velocity telescope plus actual weighted E_a support; SF4/SF5, not independent endpoint extrema. Existing Bs+Bf constants remain conservative one-step/mean outer bounds only. |
| Magnetic interaction | Unchanged applied-event information; retain frames and gyro kernels in SF5/SF8; do not confuse restricted service with full physical gauge exclusion. |
| Capture/Live/H18 | Carry both components, primitives and frontend state through the real chronology; held BA does not hold physical slow truth or clear fast accumulation. |
| H18 release | Same slow truth and fast history; estimator release changes permissions/covariance only. Held-H18 LIN BIBO and captured-domain release compactness are CLOSED qualitatively; general capture and sharp outer-entry supply remain OPEN. |
| A21 persistence | Homogeneous compatibility equations remain source-faithful; test attainable base histories using SF1--SF8, not zero-base-innovation substitution. |
| Outer-to-inner entry | Linked chi/gamma<r_in² plus every-prefix retention over the new reachable class; no model-only promotion. |
| Finite-error supply | Slow recurrence SF6 and fast functional SF5; all other existing defects remain linked. Unknown H,C blocks a numerical source-uniform temporal supply. |
| STILL/TRANSITION/MOVING | Physical decomposition and U differences survive switches; the amplitude-only quiet monitor is only a necessary screen, never a fast-history or STILL certificate. |

The coupled (tau,sigma_aw,R_S,T_S) word is unchanged. Existing conditional
radius-local field-axis, operation-energy and Schur identities remain usable
on a qualified subfamily. The .15-ball is not an entry certificate. No obsolete
O1/O2 scalar-kernel architecture controls this migration.

## SF9. Qualification still required

A future assembled-device qualification needs independently referenced true
acceleration/angular rate on the delivered timestamp grid, one explicitly
carried slow/fast decomposition, amplitude/rate checks, and signed window
accumulation at all placements covering the selected horizons. Finite data are
finite-prefix evidence only; reference uncertainty, calibration scope, unseen
operating conditions and the all-time continuation claim require separate
justification. Motion-correlated calibration errors cannot be discarded by
subtracting the observed mean independently from every window.

`BiasContinuationCertificate` requires both temporal profiles and evidence IDs;
missing profiles fail closed even when the older calibration/all-time booleans
are true. `audit_bias_trace` reports an OPEN prefix when temporal parameters are
missing. No numbers were chosen to reject either witness or to make the theorem
pass. End-to-end stability, release/outer entry and temporal device qualification
remain OPEN for the precise reasons above.
