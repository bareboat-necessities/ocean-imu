# ALT bounded-input WPE arithmetic supply theorem

`finite_wpe_uniform_bounds.py` proves an inductive bound from literal reset
for every finite WPE input sequence with `|x|<=32 m/s²`, default construction
constants, `dt=RN32(0.005)` and the named scalar RNE transcendental relations.
It covers each local permitted FMA/separate moment outcome, early returns,
independent log production and the one-way usable latches. It is not a trace
bound or a contraction result.

## Connection to the startup source

For the declared sensor profiles, `finite_startup_sensor_contract.py` proves
that a **defined** scalar Mahony update produces this bounded vertical input.
The inverse-square-root proof covers every finite nonnegative norm word using
16 complete mantissa cells for each of the 255 finite exponent fields. The
partition includes zero and subnormals. It encloses the actual integer seed
and Newton operations, not an ideal reciprocal square root.

The scalar normalization enclosure is at most 1.1101404428482056 for
`RN32(RN32(n*y)*y)`. Charging the four-term norm accumulation and four
quaternion multiplications gives `||q||² < 1.112`. The exact quadratic
Mahony down row satisfies `||d(q)||=||q||²`; its source-order rounding and
the same packet's norm bound imply `|vertical_accel|<32`. This does not need
a tilt-accuracy hypothesis. It does need the preceding seed/feedback/Euler
operations and norm sum to be defined, which remains a separate obligation.

The profile is carried from startup construction through every executed
guard/observer/WPE event and pending boundary. goLive retains that history.
The qualified path checks the source packet, dormant guard, same machine
vertical output, WPE bounds and log/exp argument ancestry. These checks do not
prove that all physical histories execute that path. An arbitrary post-Live W
cannot be substituted for the bounded startup profile.

## Uniform induction

Let `u=2^-24`, `eta=2^-150` and `R(z)=(1+u)z+eta` for nonnegative z.
This is an absolute binary32 rounding envelope, including gradual underflow.
The fixed leak coefficient is derived from rigorous exp enclosures at the
actual rounded argument; both enclosure endpoints select the same RNE value.
The default lambda retains the shipping float multiplication of 2π and 1/50.

| Persistent state | Absolute upper bound |
| --- | ---: |
| Previous input | 32 |
| First high-pass state and its previous value | 2^17 |
| Second high-pass state | 2^29 |
| Velocity V | 2^33 |
| Elevation E | 2^37 |
| Velocity/elevation mean | 2V / 2E |
| Velocity/elevation second moment | 2V² / 2E² |
| Weight | 2 |
| Elapsed time | 2^17 |
| Initialized canonical log-period | 48 |

Substitute these bounds in each literal recurrence using R at each rounded
node. All margins are strictly positive. Examples are >18 for the first
high-pass state, >75,000 for the second, >2.7 million for velocity and
>43 million for elevation. Weight has margin >2.7·10^-5. The moment envelope
is affine in alpha, so checking its two endpoints covers every rounded horizon
between the default 20 and 180 seconds. In particular, second moments retain
the literal `(alpha*v)*v` multiplication order. Elapsed time is monotone and
`RN32(2^17+RN32(0.005))=2^17`, proving its bound from reset.

The literal positive weight/variance/omega gates bound every division that is
actually executed. Clamping a second moment minus a nonnegative square cannot
increase that second moment. These facts give

`2^-57 < raw_period < 2^16`.

No positive variance or eventual period production is assumed. Every
nonproducing branch retains its existing log and usability state. Range
reduction and the positive-tail atanh series enclose log on the complete
interval, giving `-40 < log(raw_period) < 12`. The actual smoothing recurrence
preserves `|log_period|<=48` with margin >0.0011. All exp arguments lie within
the proved range [-60,60], and all arithmetic intermediates have magnitude
below 2^120. The conservative numbers prove finite supply, not useful gain.

Every consumed period/frequency exp is checked against its own stored-log
argument, including the eager frequency getter before the usable latch.
Frequency and period are independent rounded exp calls, not exact reciprocals.
Raw log is checked against the period derived from the same machine moments.
Neither a detached finite witness nor a nonfinite return qualifies this model.

## Parameterized Live supplies

The same induction covers input envelopes `32·2^k`, for integers `k=0,...,5`.
Linear states scale by `2^k`, squared moments by `4^k`, and the raw-period lower
bound by `2^-k`. The raw-log lower bound becomes `-40-k`; the canonical stored
log still satisfies `|log_period|<=48`. Every induction, logarithm and overflow
margin remains strictly positive. At `k=5` the derived largest-intermediate
enclosure is `3.022315e38`, below binary32's largest finite value by `3.805086e37`.
For `k>=6` this enclosure cannot prove totality. That is a limitation of this
bound, not physical source rejection or evidence about rho.

The declared MEMS input contract admits raw body acceleration up to
`160 m/s²` per axis and gyro up to `35 rad/s` per axis, before executing the
filter. Under the retained transparent guard scope, the source-locked scalar
observer proof gives `|vertical_accel| < 321.169 m/s²`. Thus the fixed `k=4`
WPE envelope (`512 m/s²`) covers every admitted packet. This is the raw Mahony
vertical value that directly drives WPE; the tracker LPF is a separate
statistical descendant. The ISS disturbance parameter W is not capped or
replaced by a startup sensor residual limit.

The observer proof uses the exact all-time integral barrier `|I_i|<=4096`:
the largest admissible increment is below half an ULP at that boundary, so
RNE cannot cross it from construction zero. Its finite Euler norm and the
normalization invariant consequently hold on every initialized finite prefix.
No 30,602-sample startup deadline or reset of the Live integral is used.

Every bounded startup WPE history automatically receives this physical Live
envelope through `begin` and the canonical admitted startup constructor. The
configuration, moment states and logs retain object identity; latches and
counters retain their values. Only the monotone proof-envelope metadata widens.
Each ensuing IMU edge admits the same raw packet before lower execution and
checks transparent guard output, zero weight and the dormant removed-RMS
predicate. The default enabled guard is permitted. MAG and HOLD retain the
same envelope and machine history; complete words require all 600 WPE updates.

The auxiliary `live_iss_certificate(W)` also proves the conditional inequality
`||raw_residual|| <= sqrt(2) W <= 1.5 W`, using the separately retained
`nu_acc` and thermal blocks. It is not used to impose a numeric W threshold on
the physical source word.

## Library-error envelope

The persistent `libm_profile` can select the strict RNE enclosure or the named
`exp-sqrt-relative-2^-20-log-absolute-2^-14` relation. The latter allows exp and
sqrt relative error `2^-20`, log absolute error `2^-14`, and the final binary32
rounding of enclosing endpoints. It does not assume correct rounding of libm.
Every actual exp, log and sqrt operand remains bound to the same machine state
and raw period. The profile is preserved across startup, goLive and every
subsequent WPE transition.
The canonical startup constructor selects the error-inclusive profile;
standalone component constructors can still select the strict RNE relation.

The exact induction recomputes all leak, gain, moment-alpha, log-alpha and root
bounds for this error relation. At physical `k=4`, its high-pass margin exceeds
`291`, log margin exceeds `0.0005648`, and all intermediate values are at least
`2.647e38` below binary32 overflow. Both separate and permitted FMA outcomes are
retained, including the target's fused `ratio - lambda*lambda` operation. The
`target_error_profile_certificate` connects the pinned newlib approximation
proofs directly: exp relative error is below `3.339e-7` on `[-60,60]`, and log
absolute error is below `5.249e-6` on `[2^-62,2^16]`. Both fit the WPE allowance
with a strict margin. Scalar/division/sqrt instruction semantics and the final
firmware's selection of those objects remain separate deployment premises.

## Remaining composition premises

The WPE supply result consumes the declared input domain, the source-locked
initialized observer induction, a valid startup seed, the transparent guard
scope, and the scalar arithmetic profile. It establishes neither startup
capture nor the complete MEKF/covariance/reset word. Those premises and actual
target correspondence remain visible in their own gates. Bounded WPE supplies
do not authorize storage search, rho estimation or ALT promotion by themselves.
