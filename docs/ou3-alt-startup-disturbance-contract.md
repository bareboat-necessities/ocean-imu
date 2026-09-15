# ALT startup disturbance contract

The post-Live ISS quantifier does not certify startup. In
`finite_admitted_imu_disturbance.py`, W bounds the forcing produced by an
already-executed Live event; its ordinals are 1..600 after a supplied goLive
result. It neither supplies a raw startup sensor bound nor proves that such a
goLive result exists. Racc and configured sensor covariances are not hard
pathwise bounds.

## Exact obstruction to an arbitrary-bounded startup extension

Take a level boat at rest: p=v=a=0, constant identity attitude, and zero physical
gyro and accelerometer biases. Zero translation satisfies corrected COMPLETE-BRMM
including every centered primitive bound. Zero bias with zero driver satisfies
BIAS0, BIAS1 and BIAS2. Choose one constant accelerometer residual

`epsilon=2^-11; n_a=(0,0,g-epsilon); g=9.80665 m/s^2`.

The same physical sensor identity gives:

| Quantity | Exact value |
| --- | --- |
| Physical specific force | `(0,0,-g)` |
| Raw measured acceleration | `(0,0,-epsilon)` |
| Residual norm | `12551887/1280000 = 9.80616171875 m/s^2` |
| Measured direction | Exactly the correct up direction |
| Binary32 norm | `1/2048` |
| Shipping binary32 seed threshold | `8589935/8589934592` |
| Norm minus threshold | `-4395631/8589934592` |

This is a finite residual, not a claim about BMI270/MPU6886 hardware noise.
It belongs to an arbitrary-bounded-residual extension and must be excluded by
a qualified startup contract if that extension is not intended. The original
small-disturbance theorem and its source assumptions are not changed here.

## Constant-input induction

The first guard sample seats its four low-pass stages at the raw vector,
the first detector stage at that vector and the second at zero. At the default
5 ms step, exact rational exp enclosures identify one RNE value for each guard
coefficient. All named separate/FMA outcomes of the constant stage recurrence
equal the same input. Hence the detector output, mean square, RMS and engagement
are exactly zero and the guard returns the raw vector unchanged. The proof
explicitly includes `sqrt(0)=0`.

`VerticalAccelComplementary::update` then returns before seeding because the
norm is below `1e-3f`. No private observer state, integral, elapsed time or
vertical output changes. Its exact successor equals reset, so induction gives
an uninitialized proxy on every such sample. The guard's mathematical sample
counter is bookkeeping, not a shipping state needed by the induction.

Both `ready_by_quality` and `ready_by_timeout` require `proxy_ready`. Granting
every other favorable predicate cannot make either true. Magnetic acquisition
does not initialize this private observer. The argument therefore excludes a
Live handoff on any defined finite prefix of this constant input, independently
of increasing the timeout. It makes no all-target libm or indefinite arithmetic
totality claim.

`finite_startup_disturbance_obstruction.py` checks the exact source identity,
guard fixed point, observer fixed point and handoff predicate, and binds the
argument to the reviewed shipping source hashes. The native regression starts
the unchanged public wrapper from construction, with and without 25 Hz valid
magnetic service. No internal state is installed. One subsequent ordinary
gravity sample checks that the missing seed magnitude caused the obstruction.
Native finite-prefix checks support correspondence; the induction supplies
the no-handoff implication for arbitrary finite prefixes in the named graph.

## Required qualification

A sufficient first-sample magnitude condition is

`g - A_max - B_max - N_max - E_max > RN32(0.001)`.

Here A bounds physical acceleration, B true accelerometer bias, N the raw
residual and E API conversion, guard displacement and norm-evaluation error.
Reverse triangle inequality proves this implication. Each term must be
source-qualified on the same startup history. This condition alone supplies
neither the seed-angle/Mahony invariant nor the later gravity-alignment,
magnetic-frame or target-arithmetic qualifications.

No numerical residual cap, informative-sample premise or hardware admission
has been selected. The 11 master prerequisites remain open. The master retains
this obstruction so that a post-Live W bound cannot be promoted into startup
observability. A startup sensor/observability contract must be specified and
proved before universal reachability and rho estimation can proceed.
