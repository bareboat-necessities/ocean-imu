# ALT physical startup timeout witness

The current commissioned source does **not** ensure Live by the first
150-second clock crossing. The unchanged full wrapper has an admitted physical
input that is still in startup at samples 30,002 and 30,602. That same input
later reaches Live. This refutes the proposed finite horizon, not eventual
startup, and supplies no later universal deadline.

## One physical history

Let `P=8.09 m`, `omega=0.64 rad/s` and `A=P*omega^2=3.313664 m/s^2`:

```
p(t) = -P (cos(omega*t), sin(omega*t), 0)
v(t) = P*omega (sin(omega*t), -cos(omega*t), 0)
a(t) = A (cos(omega*t), sin(omega*t), 0).
```

True body attitude is `Rz(psi(t))`, where `psi(0)=0` and the yaw rate is:

| Physical interval | Yaw rate, rad/s |
| --- | ---: |
| 0–5 s | 0.20 |
| 5–125 s | 0.61 |
| 125–130 s | 0.41 |
| 130–135 s | -0.20 |
| After 135 s | 0.20 |

Yaw is continuous. Its derivative exists almost everywhere; the source has a
rate bound, without an angular-acceleration bound. Roll, pitch, lever arm and
true sensor biases are zero. Zero bias and driver solve BIAS0, BIAS1 and BIAS2.

The raw ideal specific force is `Rz(-psi)*(a-g*e3)`. The gyro adds the bounded
body residual `0.019*(cos(0.1*t),sin(0.1*t),0)`. The native regression emits its
actual binary32 API values. Exact rational Taylor/range-reduction enclosures
check **every** emitted accelerometer packet against this same physical history
with residual norm at most `0.0001 m/s^2`; an exact squared-norm check bounds
every gyro residual below the smaller `0.02 rad/s` profile. These residuals
therefore fit both commissioned sensor profiles, including the fast-error
allowance. This packet audit verifies restrictions of an already specified
analytic source; sampled extrema do not define source membership.

## Continuous temporal and physical bounds

Write `D=sqrt(g^2+A^2)`. Before the tiny audited accelerometer conversion error,
the world measured-direction error has the exact decomposition

```
m = (0,0,g/D-1)
xi(t) = (-A*sin(omega*t), A*cos(omega*t), 0)/(D*omega)
r = m + d(xi)/dt.
```

A continuous bounded extension of the audited conversion residual changes the
normalized direction by at most `2*0.0001/(D-0.0001)`. Put that difference in
`m(t)` while retaining the same `xi`.

| Quantity | Certified bound |
| --- | ---: |
| Position norm | 8.09 m < 8.1 m |
| Velocity norm | 5.1776 m/s < 5.5 m/s |
| Acceleration norm | 3.313664 m/s² < 8.8 m/s² |
| Centered position primitive | 25.28125 m s < 1100 m s |
| Body rate | 0.61 rad/s < 35 degrees/s |
| World direction mean, including conversion allowance | < 0.053 < 0.10 |
| World direction primitive | < 0.501 s < 1.5 s |

There is no displacement DC or source restart. The translation frequency lies
inside the retained physical envelope. No pre-Live magnetic call is required
by the current service schedule, so the regression makes none. Its quality
handoff cannot bypass the failed timeout predicate through a north lock.

## Literal wrapper result and scope

The host test uses default construction, public update calls and no installed
state. The actual guard remains transparent: maximum detector RMS is below
`0.000131 m/s^2`, and engagement remains exactly zero. The proxy is initialized.
A parallel evaluation of the literal 12-second gravity LPF from its public
proxy quaternion gives:

| Update | Wrapper Live | Gravity LPF z |
| --- | --- | ---: |
| 30,002 | false | +2.84663 m/s² |
| 30,602 | false | +4.35638 m/s² |
| 33,447 | first true | approximately -0.00313 m/s² |

The last row corresponds to 167.235 seconds of update duration. The test checks
both the wrapper status and the entire input ancestry; it does not replace the
wrapper decision with the parallel LPF. Compiler-specific numerical output is
native evidence, not qualification of every ESP32/Eigen/libm target.

The failure invalidates treating the 30,602-update budget as a universal
startup-plus-600-step-word horizon. Increasing that budget to this particular
recovery time would be fitting a witness, not proving a universal bound.

## Eventual-capture architecture

An exact conditional bad observer invariant also exists in the named scalar
binary32 graph. With `c=16748919/16777216`, take
`q=(0,c,0,0)`, integral `(0,0,-RN32(0.2))`, quiet translation, true yaw rate
`RN32(0.2)`, raw gyro equal to that rate and ordinary gravity acceleration.
All feedback cross products vanish, gyro plus integral cancels exactly, and
the inverse-square-root normalization reproduces the same quaternion. The
normalized projected acceleration is `+g*e3`; a positive gravity LPF remains
on the unaligned branch. `conditional_inverted_axis_fixed_point` verifies the
exact operation ledger and invariant modulo elapsed time.

This invariant is **not** the first-sample seed and has not been reached from
an admitted startup history. The finite timeout witness still has horizontal
integral norm about 0.123 rad/s near 150 seconds, exceeding the gyro-residual
budget needed to cancel it directly. A valid indefinite obstruction needs
actual steering of that integral and attitude while preserving physical
position, velocity, centered primitive and guard conditions. Conversely, a
capture proof must show the true startup histories avoid such invariant sets
and produce alignment by a source-uniform later time. Neither conclusion is
established here; all universal startup and rho gates remain false.

A larger conditional target exists in the continuous observer. Set
`h=(20/101,0,-99/101)`, true yaw rate `-0.6 rad/s`, raw specific force
`(20g/99,0,-g)`, and integral `I=0.73*h+0.6*e3`. Measured down is `-h`,
feedback is zero, and corrected gyro is `0.73*h`, so the observer down remains
fixed while projected force stays positive. The true acceleration rotates at
0.6 rad/s, with position norm below 5.504 m, velocity below 3.303 m/s and
centered primitive below 18.345 m s. Its direction mean is `2/101` and
primitive `100/303 s`. `conditional_circular_antialignment` checks this exact
rational algebra. Thus horizontal integral above the gyro-residual cap alone
does not exclude a bad invariant: physical yaw can support a tilted circular
one. A source-compatible transition into this different circle and binary32
preservation are both unproved; the construction cannot replace startup
reachability.

## Constructive docking feasibility

`startup_antialigned_docking.cpp` fixes a second physical circle with
acceleration amplitude 2.9 m/s² and angular frequency 0.6 rad/s for the entire
history. Its position, velocity and centered-primitive norms are respectively
`145/18 m`, `29/6 m/s` and `725/27 m s`. Thirty explicitly stored rational yaw
rates steer the initialized observer; every rate is below 0.592 rad/s. At
150.01 s, physical yaw changes to the circle's 0.6 rad/s, making ideal body
specific force constant. This does not splice physical position or velocity.

The unchanged full wrapper reaches sample 30,002 with LPF z=+2.534785 m/s².
The observer is nearly antiparallel to measured down, and physical yaw cancels
the transverse part of its accumulated integral. An explicit bounded residual
controller then sets acceleration along the observer's current down and holds
that down near its docked value. Its computed Mahony feedback is below
`3.73e-8` in the 40,000-update native regression; every integral increment
rounds to zero. The shadow observer used to construct these packets is
initialized from the same first packet and matches the wrapper's public proxy
quaternion exactly at every update. No filter state is installed.

The finite tail remains non-Live, with LPF z above 2.53 m/s², acceleration
residual below 0.00705 m/s², gyro residual below `4.2e-6 rad/s` and guard RMS
below 0.000135 m/s². The test audits every actual binary32 packet against the
same analytic physical circle and yaw history. `docking_source_bounds` checks
the physical bounds and the continuous direction budget conditional on the
stated residual bound. These are finite constructive feasibility results.
An all-time invariant for quaternion direction, literal integral deadzone,
residual bounds and the dormant guard is still required; the native prefix
does not prove indefinite failure to capture.
