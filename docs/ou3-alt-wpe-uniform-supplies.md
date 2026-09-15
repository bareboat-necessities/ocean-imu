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

## Open qualifications

The named RNE relation is a theorem arithmetic hypothesis, not a certificate
for a device libm or compiler. Full source-to-observer totality, actual target
correspondence, Live source continuation and the complete machine word are
still required. `WPE_source_uniform_machine_supply_bounds` therefore remains
open in the top-level master. The bounded-input theorem is a checked component
available for that composition; storage search and rho estimation remain false.
