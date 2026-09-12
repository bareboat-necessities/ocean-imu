# ALT finite physical-source continuation

## Scope and controlling proof obligation

`finite_source_bound_live_word` carries the physical continuation alongside the
full finite filter/covariance, frontend, calibration and control state. Before
an IMU event, `SOURCE.append` checks the next physical segment; that exact
segment and its raw packet then reach the existing finite predictor. Magnetic
and hold events preserve this physical endpoint and do not restart its moments.

The source ports entering the future inequality

`V_next <= rho V + w'Gamma w + c'Beta c`

include the physical `J0/J1/J2` moments in
`finite_physical_prediction.linear_prediction_polynomials`. They must not be
three independent disturbance vectors. The result below constructs their
same-prefix source graph; it is not a storage or contraction certificate.

The legacy `Qualified*` class names denote **checked necessary finite outer
constraints**, not membership of a complete BRMM or BIAS generating history.
`certified_root` validates a declaration; it cannot certify a newly named history
by checking its strings. This distinction remains explicit in readiness reports.

## Exact moment projection and concatenation theorem

On a physical interval of duration h>0, define the three vector moments

`J0 = integral a(s) ds`,
`J1 = integral (h-s) a(s) ds`,
`J2 = integral (h-s)^2 a(s)/2 ds`.

For one spatial axis, the Gramian of these three kernels is

```
G(h) = [ h       h^2/2   h^3/6  ]
       [ h^2/2   h^3/3   h^4/8  ]
       [ h^3/6   h^4/8   h^5/20 ].
```

It is positive definite for every h>0: its quadratic form is the integral of
the square of a nonzero quadratic polynomial. Put

`E(J,h) = sum_axis J_axis' G(h)^(-1) J_axis`.

Orthogonal projection gives `E(J,h) <= integral ||a(s)||^2 ds`. Equivalently,
with `x=(J0/h,J1/h^2,J2/h^3)` and

```
G(1)^(-1) = [  9   -36    60 ]
           [ -36   192  -360 ]
           [  60  -360   720 ],
```

`E(J,h) = h * sum_axis x_axis' G(1)^(-1) x_axis`.
The sum is across all three axes, not three independent acceleration budgets.
The existing `ou3_brmm_acceleration_moment_iqc` supplies this exact inverse.

For adjacent durations a,b>0, let

```
T(b) = [ 1       0   0 ]
       [ b       1   0 ]
       [ b^2/2   b   1 ].
```

Translation of the integration kernels gives, componentwise across the three
spatial axes,

`J_ab = T(b) J_a + J_b`,
`G(a+b) = T(b) G(a) T(b)' + G(b)`.

Both identities are polynomial; coefficient tests cover arbitrary durations,
and concatenation associativity is checked with all 27 moment coordinates and
both shift durations as indeterminates, not sampled trajectories.

For each spatial axis set `u=G(a+b)^(-1) J_ab` and

`r_a=J_a-G(a)T(b)'u`, `r_b=J_b-G(b)u`.

Expansion, using the two identities, proves the exact nonnegative loss formula

`E(J_a,a)+E(J_b,b)-E(J_ab,a+b)`
`= sum_axis [r_a'G(a)^(-1)r_a+r_b'G(b)^(-1)r_b] >= 0`.

Consequently, by induction on any finite number of adjacent segments,

`E(J_prefix,H) <= sum_i E(J_i,h_i) <= A_max^2 H`.

The prefix budget is computed from those segment moments. It is not a free
input supplied to the next event. The physical endpoint identities remain

`v_H=v_0+J0_prefix`,
`p_H=p_0+H*v_0+J1_prefix`,
`S_H=S_0+H*p_0+H^2*v_0/2+J2_prefix`.

The same one-time Live origin is retained in every `PhysicalSegment`. Checking
an endpoint alone does not reset S or prove a new wave potential. The universal
argument is over real moment vectors; exact-rational executions are algebraic
checks, not a restriction of the physical source to rational-valued histories.

At the fresh Live origin there is intentionally no predecessor segment.
`QualifiedPhysicalOrigin` checks that actual endpoint directly: its time and
Live origin coincide, centered S is zero, the carried history/BIAS labels agree,
and the same physical vector and true-bias caps hold. This lets an asynchronous
magnetic event at sample zero use the real endpoint while leaving transition 1
unconsumed. It is a necessary outer-source check, not a proof that the endpoint
extends to an admitted generator history.

## Executable necessary physical constraints

The finite source constructor checks the existing declared physical vector caps:
`||a||<=8.8 m/s^2`, `||v||<=5.5 m/s`, `||p||<=8.1 m`, and
`||S_L||<=1100 m*s`, at both endpoints. These constrain physical source values,
not filter errors. It also checks the coupled moment IQC on every segment.

For projective quaternions q0,q1, the necessary rotation constraint is

`4*(1-<q0,q1>^2/(||q0||^2 ||q1||^2)) <= (Omega_upper*h)^2`.

The left side is `4 sin^2(theta/2) <= theta^2`. The declared 35 deg/s rate
is conservatively enclosed by `Omega_upper=35*(22/7)/180=11/18 rad/s`.
The bound `pi<22/7` follows from the positive integral
`integral_0^1 x^4(1-x)^4/(1+x^2) dx=22/7-pi`.
This predicate is invariant under quaternion sign and projective scale. It is
necessary, not proof of one continuous angular-rate realization. The raw packet
check bounds physical `omega_sample`, not the bias-corrected estimator rate or
the gyro residual; neither unknown estimation error nor sensor noise is clipped.

BIAS0, BIAS1 and BIAS2 are invoked separately. The continuation retains one actual
selected phi at fixed 5 ms duration, not merely a common parameter-token string.
Every segment checks driver component/norm caps and true-beta component/norm
caps. Existing outward floating endpoints are converted with `Fraction.from_float`
exactly, avoiding an inward decimal conversion. BIAS2 retains phi=1 exactly.
A recurrence and these necessary envelopes do not establish the full generating
function, thermal/strain history, or BIAS1 sinusoidal-parameter membership.

## Unclosed boundary

No complete physical-source membership is inferred. The remaining obligations
include generator/potential realization and Q/O continuation, the full BIAS
parameter/driver graph, source-to-runtime exp/trig/solver/roundoff correspondence,
the sample-zero magnetic endpoint, startup/ungauged capture, and continuation
beyond the currently represented 600-transition container. Native clock and
counter lifetime, every-prefix filter-chart retention, common joint24 storage,
and the ultimate bound remain separate obligations.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
No storage search is authorized by these source checks.
