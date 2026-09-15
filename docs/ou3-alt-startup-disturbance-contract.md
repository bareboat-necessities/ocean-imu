# ALT startup disturbance contract

The theorem uses deterministic errors throughout startup, with the two
commissioned profiles below. The numerical source is
`tools/stability/ou3_alt_startup_sensor_domain.json`. These are conditional
installation requirements selected using device evidence, not hard bounds
guaranteed by a datasheet or inferred from covariance. No calibration that is
not actually applied is credited. The independent P2/P3/P4/P5 route is unchanged.

## Declared input profiles

All sensor residual limits are Euclidean **vector** norms, not per-axis limits.
They include fast error, residual scale/cross-axis/nonlinearity and acquisition,
conversion, delay and modeling errors. True offset and slow drift remain in the
same BIAS0/BIAS1/BIAS2 history; they are not removed or counted as zero.

| Startup requirement | BMI270 / AtomS3R | MPU6886 / AtomS3 |
| --- | ---: | ---: |
| Accelerometer residual | 0.30 m/s² | 0.50 m/s² |
| Gyroscope residual | 0.02 rad/s | 0.03 rad/s |
| Commissioned residual scale/cross-axis operator norm | 0.005 | 0.01 |
| Fast accelerometer error allowance | 0.12 m/s² | 0.15 m/s² |
| Acquisition/conversion/delay/model allowance | 0.08 m/s² | 0.16 m/s² |
| Sum of accelerometer allowances at specific force 18.60665 m/s² | 0.29303325 m/s² | 0.4960665 m/s² |

The physical acceleration cap remains 8.8 m/s², gravity 9.80665 m/s²,
body rate 35 degrees/s and sample period 5 ms. Retain the 0.13 m/s²
component envelope for true accelerometer bias, so its norm is at most
sqrt(3)·0.13 < 0.226 m/s². True gyro bias has initial norm at most
0.01 rad/s, norm at most 0.02 rad/s throughout startup and rate at most
10^-5 rad/s². All three BIAS families retain their existing dynamics and
shared history. Their admission must cover the actual uncompensated sensor.

Finite, fresh, correctly scaled, unsaturated samples are required from the
first filter update, after sensor startup transients have settled. There are
no unreported missing samples. The existing zero-heel, zero-lever-arm and
dormant, transparent guard scope remains. A device qualification record must
identify actual temperature, supply, mounting and applied calibration. A
board that fails these requirements is outside this conditional theorem.

## Datasheet basis and limits of the evidence

The [Bosch BMI270 datasheet](https://www.bosch-sensortec.com/media/boschsensortec/downloads/datasheets/bst-bmi270-ds000.pdf),
revision 1.6, March 2026, pp. 11–15, reports 1.51 mg RMS accelerometer noise
at 200 Hz in the stated mode and 0.07/0.09 degrees/s RMS gyro noise.
It also lists offset, temperature, cross-axis and mounting effects. Its
min/max convention is statistical (±3σ), not an absolute pathwise guarantee.
The chosen peak limits provide headroom above those RMS figures and reserve
explicit error allowances; they are commissioning requirements.

The [InvenSense MPU6886 datasheet hosted by M5Stack](https://m5stack.oss-cn-shenzhen.aliyuncs.com/resource/docs/datasheet/core/MPU-6886-000193%2Bv1.1_GHIC_en.pdf),
DS-000193 revision 1.1, pp. 9–10, lists 1.0/1.7 mg RMS acceleration noise
and 0.04/0.1 degrees/s RMS gyro noise at 100 Hz bandwidth. Several limits are
characterization values. Board-level accelerometer offset can reach ±100 mg
per axis, exceeding the retained bias envelope. Therefore the larger residual
profile does not automatically admit an uncommissioned board. Neither sensor
profile converts a sigma multiple or finite recording into a universal bound.

## Temporal gravity information is a separate premise

Let R be the true body-to-world rotation, e3 world down and a_raw the same
measured acceleration used above. Define the total direction **error**

`r(t) = R(t)(-a_raw(t))/||a_raw(t)|| - e3`.

Require one same-history decomposition `r=m+d(xi)/dt` almost everywhere,
with `||m||<=0.10` and an absolutely continuous `xi` satisfying
`||xi||<=1.5 s`. These budgets include true bias and residual noise, not just
ideal physical acceleration. They bound persistent false gravity and its
oscillatory integral without requiring a derivative bound on the noise.
They are additional source hypotheses; a datasheet, the peak caps and BRMM
translation alone do not establish them. Packet checks verify necessary peak
and ancestry restrictions, not complete temporal membership.

The requirements hold throughout actual startup. They can be reused on a
Live word only if they continue to hold there. An arbitrary post-Live ISS
disturbance bound is not silently restricted by these startup constants.

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

## Proved seed magnitude and remaining capture obligation

A sufficient first-sample magnitude condition is

`g - A_max - B_max - N_max - E_max > RN32(0.001)`.

Here A bounds physical acceleration, B true accelerometer bias, N the raw
residual and E API conversion, guard displacement and norm-evaluation error.
Reverse triangle inequality proves this implication. Each term must be
source-qualified on the same startup history. This condition alone supplies
neither the seed-angle/Mahony invariant nor the later gravity-alignment,
magnetic-frame or target-arithmetic qualifications.

`finite_startup_sensor_contract.py` proves the reverse-triangle bound and the
scalar RNE conversion, nonnegative sum of squares and square-root margins.
On the retained transparent-guard branch, **every admitted sample passes the
literal seed magnitude test**, without selecting a favorable startup sample.

| Consequence | BMI270 profile | MPU6886 profile |
| --- | ---: | ---: |
| Exact raw norm floor | 0.48065 m/s² | 0.28065 m/s² |
| Computed norm floor, rounded down for display | 0.48064878 m/s² | 0.28064880 m/s² |
| Margin above compiled 0.001 threshold, rounded down | 0.47964878 m/s² | 0.27964880 m/s² |
| Physical seed tilt bound, before seed arithmetic | approximately 71.99° | approximately 76.26° |

The seed-angle expression is `asin((8.8+0.226+N)/9.80665)`. It is a
geometric bound, not qualification of Eigen's near-antiparallel SVD branch.
Even with its original, smaller temporal budgets, the retained Mahony metric
would require seed levels about 2.05463/2.27142, above both its outer level
1.7689 and chart ceiling 1.87267. Its original conditional certificate remains
valid under its own premises; it does not cover these new sensor profiles.
Raising that level alone cannot close capture. A new capture argument is needed
for actual gravity alignment and timeout production.

The same module covers all nonnegative finite norm words of the literal Mahony
inverse-square-root algorithm, including zero and subnormal words. After a
defined scalar Mahony update, its quaternion norm and the actual down-row
projection give `|vertical_accel|<32 m/s²` for both sensor profiles, without
assuming accurate tilt. This supplies the conditional WPE induction in
`ou3-alt-wpe-uniform-supplies.md`. Totality of the preceding observer operations
and target/compiler correspondence remain open.

The master retains conditional components and 11 open top-level prerequisites.
The seed and WPE sublemmas do not close universal startup, complete
source-uniform arithmetic or the complete 600-step word. Rho estimation and
all ALT PASS gates remain blocked.
