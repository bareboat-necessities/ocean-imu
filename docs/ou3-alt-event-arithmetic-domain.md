# ALT finite event arithmetic domain

The startup and all-event deployment-supply qualifications require a defined
machine successor before a bounded displacement from the exact word can exist.
Two different domain issues must be kept separate.

## Complete finite binary32 lattice

`finite_binary32_arithmetic.rn32` retains gradual underflow. For
`0 < |x| < 2^-126`, its value is

`sign(x) * 2^-149 * RNE_integer(|x| / 2^-149)`.

The ordinary binade calculation continues above that threshold. Both formulas
agree at the smallest normal. Midpoints use even integer lattice indices;
the midpoint between the largest subnormal and the smallest normal therefore
rounds to the smallest normal. Signed-zero encodings are identified only in
this rational value projection. Sign-sensitive expressions must use the
existing bit ledger, not infer an encoding from rational zero.

This yields the usual finite-result estimate
`|RN(x)-x| <= 2^-24 |x| + 2^-150`, including underflow. The absolute term cannot
be dropped from future event-supply inductions. Separate multiply/add and FMA
remain different operations: with `p=3*2^-149`, `target=0` and `alpha=1/2`, the
shipping EMA gives `2^-149` in separate mode and `2*2^-149` in fused mode.
The sigma/guard RNE-cell helper and scheduler predecessor use the same lattice.
An exact singleton interval at a midpoint belongs only to the even endpoint.

A positive lower bound on every nonzero sensor or frontend component is not
part of the declared source. Even a normal vertical input `2^-62` produces a
subnormal squared energy in the actual stillness relation. Excluding all
subnormal intermediates is therefore an invalid route to universal supply
qualification. This was a proof-model defect, repaired without changing
shipping or adding a source premise. Overflow, NaN handling, the target's
underflow mode, and target libm remain separate qualifications.

## Bounded raw residual is not an execution-domain theorem

The commissioned input contract in `ou3-alt-live-input-contract.md` now bounds
each raw gyro component by 35 rad/s and each raw acceleration component by
160 m/s² before execution. The following obstruction remains a regression
against omitting that input condition; it is outside the selected contract.

`finite_live_disturbance_overflow.py` constructs the same stationary physical
history with zero wave motion and zero BIAS0/1/2 bias. Startup uses zero gyro
residual and only float-conversion accelerometer error. The first post-Live
packet has gyro residual `(2^80,0,0)`; every later packet is quiet. Every raw
component is finite binary32, and the residual history has a finite amplitude
ceiling. The gyro pulse leaves accelerometer guard engagement at zero.

That input-side description does **not** imply a finite shipping execution.
For quiet aligned private Mahony, even the weak pre-pulse bound `q0 >= 1/2`
makes the literal Euler component

`RN((1/2) * RN(2^80 * RN(dt/2)))`

larger than `2^64`. Squaring it exceeds binary32's finite RNE overflow threshold
`(max_float + 2^128)/2`. Independently, the MEKF angular norm and rotation
expressions receive the same large rate. The exact report checks the threshold;
the native regression supplies the actual construction-to-Live execution.

The unchanged public wrapper reaches Live at sample 30,002 on the scalar host
profile. After the pulse, its MEKF attitude and private vertical acceleration
are nonfinite and stay nonfinite through the 600th post-Live call. The guard
stays dormant. The public private-observer quaternion getter returns identity
when its underlying quaternion norm is nonfinite; that fallback does not repair
the internal observer state or its vertical output.

This is a host execution counterexample to deriving finite execution from
arbitrary bounded raw residuals. It does not claim that the pulse is a possible
BMI270/MPU6886 hardware output, or falsify the commissioned startup profiles.
Their smaller residual caps apply during startup. The independent commissioned
post-Live input envelope now excludes this pulse. Hardware admission is not
inferred from the magnitude ceiling.

`BoundedHistory` bounds a finite forcing from an executed event.
`RestrictedForcing` for this nonfinite execution is not constructed or claimed.
Rejecting the raw history solely because no such finite result exists would
assume the totality that the prerequisite is supposed to prove. A theorem over
already-total executions is a different, conditional statement.

## Failure analysis and next choice

The input-side bound is now explicit. The corresponding scalar private-observer
induction closes on every finite initialized prefix using an integral rounding
barrier, and produces a WPE input bound below 322. The remaining all-event limiter is the full carried MEKF and
covariance/solver recurrence plus target correspondence. Larger symbolic ISS
bounds do not supply that missing recurrence.

The compared routes were a commissioned raw sensor envelope, a conditional
theorem over already-total executions, and a shipping recovery change. The
user-authorized MEMS envelope is the selected route. It is independently
checked before execution, so it avoids circular output-based admission without
changing shipping. The startup supply prerequisite still needs the remaining
seed/guard/frontend/tuner/Racc and target qualifications.

The next falsifiable experiment is to propagate the declared input contract
through the complete same-history machine covariance/solver and branch
relations. A bound on the actual input is now available, while a norm ceiling
on already-produced outputs remains insufficient. All ALT PASS and
storage-readiness gates remain false.
