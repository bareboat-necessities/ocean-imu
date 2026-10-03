# Quiet compatibility and the MOVING span calculation

## Status and controlling inequality

Regional practical stability remains **OPEN**. The quiet absolute-entry obstruction in
`ou3-quiet-inner-entry-obstruction.md` remains controlling. This note tests the
proposed implication from MOVING gravity span and persistent SLOW+FAST bias to
physical tilt/BA separation, before inserting a positive margin into

`V_N <= (1-gamma) V_0 + chi_gamma`.

The exact result below is **B for that two-epoch separation mechanism**. It is an
all-time physical sensor-compatibility identity, not an independently selected
nominal trajectory. All-time admission under *actually applied* MAGNETIC SERVICE
is still OPEN. Consequently this note does **not** claim a fully admitted shipping
MOVING counterexample, refute all MOVING entry theorems, or claim instability.

The analytical calculations are reproduced by `moving_quiet_compatibility.py` and
`moving-quiet-compatibility-certificate.json`. `moving_compatibility_diagnostic.py`
and `moving-compatibility-carried.json` are a separate finite float32 replay.

## Q1. An exact quiet fibre, not an arbitrary Euclidean ball

For a constant ideal magnetic record `y_m=R^T B`, constant delivered accelerometer
record `y_a`, and physical stillness, put `d_a(R)=y_a+g R^T e_z`. The existence of an
indefinite admitted SLOW+FAST decomposition of this constant total error is exactly

`||d_a(R)|| <= B_a,s + C_a/H_a`.

Necessity: average over any complete H_a window. The slow average has norm at most
B_a,s and the FAST average has norm at most C_a/H_a. Sufficiency: choose the
constant slow vector to be the radial projection of d_a onto the B_a,s ball and
put the remaining constant vector in FAST. Its norm is at most C_a/H_a<B_a,f and
it satisfies every placed short-window integral bound. This uses one persistent
split; nothing is reselected at a proof-word boundary. The analogous condition for
constant gyro record y_g is `||y_g||<=B_g,s+C_g/H_g`.

For `y_a=-g_m e_z`, `y_m=75e_x`, `B=75e_x`, the exact attitude fibre is

`R=R_x(theta)`,
`d_a=(0,g sin(theta),g cos(theta)-g_m)`,
`(g-g_m)^2+4g g_m sin(theta/2)^2 <= (B_a,s+C_a/H_a)^2`.

Here g_m is the exact binary32 shipping literal, not a replacement for physical g.
The zero-FAST subfibre uses B_a,s instead of B_a,s+C_a/H_a. The prior quiet witness
lies in that smaller subfibre. This is an exact statement for the constant ideal
magnetic subcase. The general quiet fibre retains the existing field/residual
histories and their magnetic-service constraints rather than setting their errors
to zero or independently varying the learned reference.

Physical p=v=a=0 on indefinite stillness about the equilibrium. The physical S
primitive retains its inherited value S_*; an S=0 estimator pseudo-measurement is
not a physical measurement of S_*. In the quiet construction S_*=0. On a later
quiet episode its allowed values come only from the preceding carried history.

To express a tube in shipping coordinates, let E(x,hat_x) be the actual 21-state
error map (including the local quaternion chart), let X_Q(y;S_*) be the constant
compatible representatives just defined, and use the **actual same-record** P_k:

`Q_still,k = {E(x_Q,hat_x_k): x_Q in X_Q(y;S_*)}`,
`dist_{P_k^-1}(e_k,Q_still,k) = inf_{q in Q_still,k} sqrt((e_k-q)^T P_k^-1(e_k-q))`.

No diagonal approximation to P^-1 is made. The local chart must be common to the
comparison states; the full nonlinear fibre is not silently replaced by its tangent.
Magnetic-reference, Mahony, guard, tuner, S and OU/AW processing are all generated
by the same delivered record. The definition does not authorize independent choices
of their states. Defining this fibre alone proves neither estimator boundedness nor
a time-uniform covariance-metric radius.

## Q2. Persistent quiet SLOW+FAST source charge

On a constant record, the FAST delivered hold equals `d-b_s,hold`. For a centered
window of length L contained in the quiet episode, L<=H, and maximum delivered gap h,

`||d-b_s(t_k)|| <= C/L + D L/4 + D h`.

Indeed, integrate `d-b_s,hold`, use the signed cap C, and compare each held slow
sample with b_s(t_k) by the same D-Lipschitz history. The centered distance integral
is L^2/4 and the hold mismatch contributes at most hL. There is no FAST derivative
assumption. Near a quiet boundary use an available one-sided window and replace
D L/4 by D L/2; a centered window cannot cross a moving episode.

Current constants give, with L_a=14 s, L_g=28 s and h=.006 s,

`eta_a = .00707742857142857 m/s^2`,
`eta_g = .000141488571428571 rad/s`.

The canonical constant accel split from Q1 differs from total d_a by at most
C_a/H_a=.000833333333333333. Thus its slow-BA comparison charge is at most
r_a=.007910761904761905. For zero gyro record the canonical slow gyro is zero and
r_g=eta_g. These are source-coordinate charges, not a fabricated estimator tube.
If J_bb is the full six-by-six bias principal block of **P_k^-1**, including its
accel/gyro cross terms, a valid same-prefix metric charge is

`r_Q(k)^2 = sup_{||u_a||<=r_a, ||u_g||<=r_g} [u_g;u_a]^T J_bb [u_g;u_a]`.

It bounds the difference to the constant compatible representative with the same
physical attitude and inherited LIN coordinates. A uniform r_Q still needs the
literal covariance bound at the relevant prefixes. This is a source projection
lemma, not proof that the estimate itself or the moving fibre is bounded.

## M1. The exact gauge mechanism

Let Q be a fixed rotation about the physical magnetic field B. Transform a physical
history by R_1=Q R_0, p_1=Q p_0, v_1=Q v_0, a_1=Q a_0 and

`b_a,1=b_a,0+R_0^T(Q^T-I)g e_z`.

Then delivered accelerometer and magnetometer records agree exactly, and the body
angular rates agree because Q is constant. The required bias shift has derivative

`d/dt [R_0^T d_Q] = -[omega_0]x R_0^T d_Q`, `d_Q=(Q^T-I)g e_z`.

It is constant whenever angular velocity is parallel to this body vector. Gravity
direction can still change: its derivative is `-[omega_0]x R_0^T e_z`, which need
not vanish. A scalar gravity-span premise does not forbid this geometry.

## M2. Explicit indefinite MOVING pair at the current constants

Set

`sin(theta)=400/40001`, `cos(theta)=39999/40001`,
`psi(t)=.02 sin(pi*t/10)`, `p(t)=.02 sin(pi*t/10)e_x`,
`R_+(t)=R_x(theta)R_y(psi(t))`, `R_-(t)=R_x(-theta)R_y(psi(t))`,
`b_a,+ = +g sin(theta)e_y`, `b_a,- = -g sin(theta)e_y`,
`b_g,+ = b_g,- = b_a,f = b_g,f = 0`, `B=75e_x`.

Take v=p', a=p'' and S=(.2/pi)(1-cos(pi*t/10))e_x; S'=p and S(0)=0. This is one persistent
continuous history in each case, not samples chosen independently. At every epoch,

`y_a = R_y(-psi)(a_x e_x - g cos(theta)e_z)`,
`y_g = psi'(t)e_y`, `y_m = 75 R_y(-psi)e_x`.

The complete delivered records are identical. In particular

`(R_+^T-R_-^T)g e_z = 2g sin(theta)e_y`

is constant, so its two-epoch difference is **exactly zero at every pair of times**.
Subtracting at more times or composing more moving windows cannot turn that zero
into a strict bias-variation charge. The slow biases do not track an oscillating
residual: both are exactly constant. Their norm is .09806404839879003, their rates
are zero, and both FAST primitives are identically zero in the real-arithmetic
source construction. No model mismatch has been dropped from the estimator error.

Every 30-second window contains a full 20-second period. Its displacement span is
.04 m>.03 m. Its gravity span is

`2 asin(cos(theta) sin(.02)) > .03999 rad > 2 deg`.

Rational checks use `sin(.02)>=.02-.02^3/6`; no floating angle threshold proves this
membership. The entire history also satisfies

`||p||<=.02 m`, `||v||<.006284 m/s`, `||a||<.001975 m/s^2`,
`||jerk||<.000621 m/s^3`, `||omega||<.006284 rad/s`, `||S||<.128 m*s`.

The field is fixed, horizontal and of admissible norm. Translation is the same in
both histories because rotations about e_x leave p unchanged. Isolated zero rates
do not split the indefinite moving episode. In native float32 replay, delivered
rounding is a separate sensor/arithmetic charge; the exact all-time identities are
statements in real arithmetic, not claims of an all-time native certificate.

## M3. Full covariance metric: a lower comparison without deleting cross terms

A deterministic shipping execution from the same initial state and identical
complete record has the same nominal state, P, K, Mahony, guard, tuner, committed
parameters, magnetic state and scheduler in both worlds. This conclusion follows
by induction through the entire algorithm, regardless of its internal coupling.

For the unmodified default BA prior, active prediction gives

`P_ba^-=phi_b^2 P_ba^+ + (1-phi_b^2)Sigma_ba`,

with stationary variance <=1/1600; held prediction leaves the marginal unchanged.
The initialization/release floor is smaller. A regular optimal Joseph correction
cannot increase a BA diagonal; held-BA zero gain rows leave it unchanged. Attitude
reset and AW covariance synchronization do not modify that diagonal. Consequently,
on this default frame/profile and regular real-arithmetic covariance branch,

`P_bay,bay <= 1/1600`.

For either physical error, full-metric Cauchy--Schwarz gives

`V_+/- >= ( +/- g sin(theta)-hat_bay )^2 / P_bay,bay`.

At every common finite SPD prefix, irrespective of any symmetry in hat_bay,

`max(V_+,V_-) >= 1600(g sin(theta))^2 = 15.386492141376374...`.

This already rules out simultaneous eventual *retained* .15 entry of the pair,
**provided both entire shipping continuations satisfy all theorem premises**.
It is not yet an unconditional refutation of that full-premise statement: actual
all-time magnetic service and execution totality have not been discharged.
The finite replay additionally has hat_bay=0 at every observed step, so both
physical errors satisfy the lower comparison there. That observed symmetry is
not used to promote the pair bound into an all-time float32 theorem.

## M4. Actually applied magnetic information: finite result and remaining gap

The native driver runs the literal full shipping wrapper, including construction,
H18, reference refinement, BA release, guard/Mahony, coupled tuner, OU/S chronology,
actual covariance/gains/Joseph updates, resets and continuous magnetic processing.
The fixed diagnostic profile is sigma_a=0.2, sigma_g=0.00135 and sigma_m=0.8.
The accelerometer value matches the current theorem and default wrapper; the
other two values match the world-frame fixture. That older fixture uses sigma_a=0.12,
so this replay is not its unchanged profile or every deployment's default.
No configured value was retuned after observing this result. Read-only temporary header taps
export F, H, K, the actual factorized innovation S and G. They do not replace the
estimator. Four homogeneous probe columns represent both true heading/axial-BG
pairs. Only these diagnostic columns restart at a service root; the estimator and
physical history never do.

In the 240-second carried run, Live begins at sample 17960, BA activation/reference
refinement at 24016, and 3756 magnetic corrections are accepted. On 102 disjoint
one-second windows after activation plus 17 seconds,

`lambda_min(sum Phi^T H_m^T (S_m^act)^-1 H_m Phi on E_hb) >= 7.0247646`.

Every checked window has 25 actually applied magnetic corrections. The minimum
observed BA-marginal storage lower comparison after t=200 seconds is 217.6418599.
This is **finite carried evidence**, not an interval theorem, not every placed
window, and not proof of an indefinitely admitted continuation. A finite successful
frontend or magnetic replay cannot discharge the missing all-time obligation.

The next falsifiable calculation is an all-prefix/all-time accepted-information
certificate for this *same* planar history, retaining the actual magnetic-reference,
frontend, tuner, covariance and scheduler recurrences. If that succeeds, M3 becomes
a lossless obstruction to universal retained MOVING absolute entry. If it fails,
identify whether a real premise is violated or only the enclosure is insufficient.
Do not classify this witness A until every admission item is closed.

## Preserved structures and scope

Structures preserved: the 21-state estimator, attitude/BG and BA forcing/projection,
LIN OU chain and S corrections, Mahony, guard, coupled adaptation, actual P/K/Joseph
and resets, startup/H18/refinement/release, continuous magnetic processing, all
physical constants and all SLOW+FAST/MARINE limits. No estimator or assumption is
changed. The analytical sensor identity needs no surrogate estimator. The finite
probe uses the complete literal same-history homogeneous factors; it is not a full
nonlinear extended-HistoryCell certificate.

Relaxations introduced: none in M1--M3's algebra or physical construction. Q2's
bias-charge product of norm balls is an explicit sufficient outer for a metric
upper bound, not an admission witness. Finite replay scope is not a relaxation
that can be promoted to an all-time theorem. No shipping instability counterexample
has been certified. The local qualified theorem, homogeneous loss identities and
forcing-aware completed-square algebra remain valid.
