# ALT same-history magnetic event composition

## Obligation and representation boundary

The controlling ALT obligation is the complete finite joint24/21-covariance
shipping word, followed by a coercive storage inequality. An asynchronous
magnetic call must not splice a new reference, gain, calibration state, physical
endpoint or private observer into that word. This note supplies the magnetic
product-state substitution and an all-history **finite-real calibration bound**.
Neither result certifies the complete runtime word, a useful storage margin,
startup capture, domain retention or deployment arithmetic. No storage search
is authorized by these results.

The supplying implementation is `src/tuner/ContinuousMagHardIronEstimator.h`
and `src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h`. The finite relations are
`finite_continuous_mag_runtime.py`, `finite_live_magnetic_word.py` and
`finite_live_interleave.py` under `tools/stability/ou3_alt_contraction`.
Shipping and the independent P2/P3/P4/P5 track are unchanged. The existing ALT
zero-wind-heel scope remains in force; continuous calibration is **not** disabled.

## One product state and one physical magnetic source

Carry the existing Live state (joint24, full 21-state covariance, physical
reference, private Mahony, WPE/band/statistics/tuner, guard, Racc, active OU/S
parameters, S scheduler and pending a_w floor), the tilt-watchdog state, and:

- one fixed physical field/hard-iron model and physical-history identity;
- continuous calibration sufficient statistics, solve clock and last estimate;
- the applied offset, latched calibration anchor and last application clock;
- the acquisition/refinement accumulator, clock and completion flags;
- the active reference generation and unchanged constructor Rmag;
- the inner update count, first-update time, lock/hold flags and call-schedule
  prefix clock.

At an actual physical endpoint the same qualified packet satisfies

`m_raw = R_true B_true + b_HI + n_m`.

`qualify_squared` checks exactly the existing MAG-BMM150-DET-v1 norm inequalities
by squaring their nonnegative thresholds. It introduces no new source bounds
and needs no rational square root of a rational vector's generally irrational
norm. The continuous estimator sees this **raw** packet. The acquisition and
MEKF consumers see precisely the corrected versions required by their stage.
The physical model is stored once, not supplied separately to each Live call.

Startup accumulation occurs ahead of the outer gravity/settle/admission gates.
Its memory is carried through the handoff, not initialized when Live begins.
`from_startup` requires the same carried private observer and checks that the
fresh nominal quaternion is the conjugate of the actual pending-gauge/proxy-tilt
seed. This is a conditional handoff composition; it does not prove the complete
startup history reaches that handoff.

## Literal Live call order

A represented Live call performs the following substitutions:

1. Apply the outer enabled/delay gate; qualify one physical raw packet.
2. Advance continuous statistics from the current private Mahony tilt. Read no
   MEKF state for this calibration. On repeated wrapper timestamps, use the
   configured fallback sample interval, as shipping does.
3. If refinement is due, reset its acquisition accumulator and clock once;
   accumulate the current corrected packet in the **private Mahony** tilt
   frame. This is what the executable source does, despite older explanatory
   comments referring to MEKF tilt. Rejections retain the updated sample clock.
4. On successful refinement, write the new reference, overwrite only nominal
   yaw using the same accepted mean, mark refinement done, and release the
   external BA hold. The internal count/time lock may remain closed. A yaw
   overwrite preserves every nonattitude nominal coordinate and the entire
   covariance; its finite attitude error is recomputed against the unchanged
   physical endpoint. It is not a Kalman/Joseph covariance reset.
5. Apply continuous hard iron, when eligible, and write its coupled reference.
   Both a refinement write and a continuous write may occur in one call; each
   increments the same generation and neither changes constructor Rmag.
6. Subtract the final applied offset from this same raw packet, perform the
   magnetic sanity/SafeLDLT/measurement/Joseph/reset/projection relation, then
   update the inner count/time/unlock control. Measurement rejection does not
   suppress that control update.

The resulting measurement discrepancy is derived, never supplied freely:

`nu_m = R_true (B_true - B_ref) + b_HI - b_applied + n_m`.

Here default startup hard-iron fitting is off, so `b_applied` is the continuous
applied offset. The code retains the total-offset expression for its explicit
state shape and rejects the unsupported nondefault startup fit.

An IMU event consumes the preceding product-state MEKF and covariance; a
magnetic event consumes that IMU endpoint and its exact private observer;
a hold event preserves the physical/frontend state and applies the literal
BA covariance/mode edge. Induction over any finite sequence proves equality
at every **represented** prefix, conditional on the supplying event relations.
Splitting the sequence at a prefix cannot restart P, a calibration accumulator,
the model, the first-mag clock, or the one-time Live S origin. The existing
watchdog composer is included, but its free transcendental/reset witnesses are
not silently qualified by this substitution theorem.

## Continuous sufficient-statistic identity

For accepted raw samples `m_i` and their proper private-tilt rotations `R_i`,
all sufficient statistics receive the same forgetting factor. With the
resulting positive weights `w_i`, define

`W = sum w_i`, `A = sum(w_i R_i)/W`,
`u = sum(w_i R_i m_i)/W`, `v = sum(w_i m_i)/W`,
`q = sum(w_i ||m_i||^2)/W`.

The exact recurrence implies these identities by induction from zero statistics;
rejected samples leave them unchanged. In particular, `||A||_2 <= 1` and

`M = I - A^T A = sum_i (w_i/W) (R_i-A)^T (R_i-A)`,

so `0 <= M <= I`. This is one weighted rotation variance, not an independent
entrywise matrix box. With shipping default ridge

`r = 0.0005 + 0.25 trace(M)/3`,

`M+rI >= 0.0005 I`; the ideal-real regularized solve is nonsingular regardless
of whether the separate **unridged** information gate succeeds. No excitation
is inferred from the ridge. On a successful eligible solve,

`(M+rI) b_fit = v-A^T u`, `B_fit = u-A b_fit`.

The residual sum is exactly

`SSE = Wq - 2 B_fit^T (Wu) - 2 b_fit^T (Wv)`
`      + 2 B_fit^T (WA) b_fit + W (||B_fit||^2 + ||b_fit||^2)`.

The finite evaluator uses the same moments for the equation, bias-fraction gate,
residual gate and later reference. An eigensolver failure is an explicit
conditional branch. A failed due solve **replaces** the prior estimate with
an invalid result; a not-due sample preserves it. The exact solver in this
Python layer is not an implementation of rounded Eigen LDLT.

## Uniform finite-real calibration bounds

Assume the named physical magnetic envelope at every source call, the actual
zero startup-offset initializer, fixed default calibration configuration, proper
normalized proxy rotations, and exact real arithmetic. These are premises of
this lemma, not conclusions from regression traces.

The source triangle inequality gives `||m_raw|| <= 75+5+2 = 82 uT`. Since `q`
is a weighted average of squared raw norms, `sqrt(q) <= 82 uT`. The shipping
acceptance gate implies

`||b_fit|| <= 0.35 sqrt(q) <= 28.7 uT`.

Invalid estimates cannot be applied. Between solve instants, a retained valid
estimate retains its bound. The default slew is a convex combination of the
old applied offset and this bounded target, because `alpha=1-exp(-dt/45)` lies
in `[0,1]`. Failed applications do not change the applied vector. Starting from
zero therefore gives `||b_applied|| <= 28.7 uT` indefinitely in the real graph.
This argument makes no claim that the fit is accurate or that it converges to
the true hard iron.

Let `C(x)=(sqrt(x_1^2+x_2^2),0,x_3)`. Reverse triangle inequality gives
`||C(x)-C(y)|| <= ||x-y||`. The same-statistics level reference is
`L(b)=u-A b`. Shipping latches the anchor before any continuous application;
therefore its anchor offset is zero. Its anchor reference `B_0` is the
provisional or completed-refinement reference, with `||B_0|| <= 82 uT`, since
continuous application is held off until refinement finishes. Successful writes
have the exact anchored form

`B_ref = B_0 + C(L(b_applied)) - C(L(0))`.

Consequently

`||B_ref-B_0|| <= ||A b_applied|| <= 28.7 uT`,
`||B_ref|| <= 110.7 uT`,
`||m_corrected|| <= 110.7 uT`,
`||nu_m|| <= 221.4 uT`.

If refinement is disabled, the same argument uses the provisional anchor.
The anchor and application clock latch before later reference-validity tests;
the finite graph retains those mutations even on failure. No repeated
accumulation of reference increments is hidden in the bound: each successful
write is rebuilt against the fixed anchor using the **current** common moments.

The last discrepancy bound is deliberately only a boundedness result. It is
not asserted small enough for a retained attitude chart or useful contraction.
The graph retains the sharper same-history vector expression for subsequent
storage work. No unknown main-filter state error is relabeled as sensor noise.
A finite-precision perturbation bound, uniform accumulator size, and applicable
clock/counter limits still have to be proved for deployment.

## Corrected count-and-time reachability

Let the first counted post-Live magnetic call be `t_1`, with
`t_1-t_L <= a`, and let every later gap be at most `g>0`. The admitted sequence
is locally finite and continues through unbounded physical time. It may include
multiple calls at an identical IMU timestamp. For count threshold `n`, the
n-th call occurs by `t_1+(n-1)g`, but it need **not** be more than one second
after `t_1`. An upper bound supplies no such lower bound.

If the n-th call is already strictly after `t_1+1`, it satisfies both guards.
Otherwise the first later call strictly after `t_1+1` exists and is no later
than `t_1+1+g`; its count is at least n. Thus the internal lock clears by

`t_L + a + max((n-1)g, 1+g)`.

For `a=g=0.04 s`, `n=250`, this remains `t_L+10 s`. A call exactly at `t_1+1`
does not qualify. The external hold can still leave H18 active indefinitely;
refinement or an explicit later release applies the separate H18->A21 edge.
`check_prefix` verifies elapsed deadlines and count continuity for finite
prefixes only; it cannot certify infinite call coverage.

This control theorem is over mathematical clocks and counters. Shipping uses
`int mag_updates_applied_` and increments it after every attempted update even
after unlock. Unbounded execution therefore also requires a target-qualified
integer-overflow argument; an unbounded Python integer does not provide one.
The floating-clock and signed-counter obligations remain open, not repaired by
restricting the physical source or by changing shipping in this proof PR.

## Validation and remaining obligations

Unit regressions check exact equations, source identity, rejection side effects,
nontrivial yaw overwrite, reference generation order, count/hold edges,
IMU->mag->IMU covariance substitution and anti-promotion guards. Their short
rational fixtures are component tests, not replay fitting, stability evidence,
or admission of all source histories. Continuous decay/eigenpair witnesses are
same-operand conditional relations; no rational test value is advertised as
an exact transcendental evaluation for deployment.

Still required: complete startup/ungauged paths, the literal tilt-angle and
preserve-yaw reset arithmetic, magnetic libm/casts/Eigen/nonfinite branches,
all remaining IMU/source/BIAS bindings, deployment lifetime arithmetic, and
source-uniform qualification of the entire 600-step master. Only then may the
existing finite-master guard authorize common joint24 storage feasibility,
every-prefix retention and the disturbance-dependent ultimate bound.

`ALT_STARTUP_PASS=false`, `ALT_LIVE_PASS=false`, `ALT_END_TO_END_PASS=false`.
