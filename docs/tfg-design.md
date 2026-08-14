# The two-frame Lie group: conventions and structure

This is the reference for `src/lie/` and the filter built on it
(`src/kalman_tfg/`). It exists because the single most common way a
Lie-group estimator gets implemented backwards is by stating one convention
choice and leaving two implicit.

Everything here is pinned by `tests/kalman_tfg/lie_group-test.cpp` and
`tests/kalman_tfg/convention-test.cpp`. If the prose and the tests ever
disagree, the tests are right.

## 1. The convention contract

Three independent choices. Never state one without the other two.

| # | Choice | This filter |
|---|---|---|
| 1 | rotation direction | `R = R_bw`, body to world |
| 2 | error definition | right-invariant, `eta = Xhat o X^-1` |
| 3 | correction side | left, `Xhat+ = Exp_G(dxi) o Xhat-` |

**"Right-invariant" describes the error, not the side the correction lands
on.** A right-invariant error is invariant when *both* states are multiplied
on the right by the same group element, and such an error is applied by **left**
multiplication. That is not a contradiction. `convention-test.cpp` asserts both
halves, including the negative cases — the right-invariant error must *not*
survive left multiplication, and left and right retraction must differ for a
generic element. Without the negatives, a degenerate error definition would
satisfy the positives trivially.

The binding identity between (2) and (3): if `Xhat = Exp_G(xi) o X`, then
`eta = Xhat o X^-1 = Exp_G(xi)` exactly, so `Log_G(Xhat o X^-1) = xi` with no
approximation. Apply the same `xi` on the right and this fails.

### Why right-invariant here

The states that distinguish this filter — `v`, `p`, `S`, `a_w`, gravity, the
magnetic reference, and the `S = 0` regularizer — are all world-frame
quantities. A right-invariant error puts all the translational errors in a
common world frame, which is why the accelerometer and magnetometer Jacobians
come out built from *fixed* references rather than from the current estimate.

A left-invariant error (`eta = X^-1 Xhat`, corrected on the right) would put
translational errors in a body frame. Gravity would appear as `Rhat^T g`, north
as `Rhat^T B_w`, the `S = 0` constraint would need transforming, and the
NED-anisotropic `R_S` would rotate with the estimate. It has real advantages —
body-frame sensor noise needs no rotation, and body biases are natural — but
they do not pay for the world-frame structure this filter is built around. That
trade would only flip if the state were redesigned around body-frame quantities.

### The OU-III bridge

OU-III makes a different choice on (1). It stores a **world-to-body** quaternion
`C = R_wb` and injects on the left:

```
C+ = Exp(dtheta) C-
```

This filter stores `R = C^T = R_bw`. Transposing:

```
R+ = (C+)^T = (C-)^T Exp(-dtheta) = R- Exp(-dtheta)
```

The same physical rotation is a **right** multiplication on `R`. So "OU-III
corrects on the left and TFG corrects on the left" is not a statement of
agreement — the rotation directions differ.

To express it as this filter's left-multiplicative correction, conjugate through
`R Exp(u) R^T = Exp(R u)`:

```
R+ = Exp(-R- dtheta) R-        i.e.   phi = -R- dtheta
```

**The correspondence is not `phi = +-dtheta`.** It is a rotated, negated
`dtheta`. This matters because both forms agree to first order at `R- = I`, so a
test that only exercised small rotations near identity would never catch the
difference. `convention-test.cpp` therefore asserts the correct bridge *and*
asserts that both naive correspondences fail, gated on `|dtheta| > 0.2` and
`|Log(R-)| > 0.2` so the check is only made where the forms genuinely diverge.

### One more trap, in the sign of `Exp`

`ocean_imu::lie::Exp` is the standard right-handed exponential: `Exp(phi)`
rotates by `+phi` about `phi/|phi|`. `ou_detail::rot_and_B_from_wt` returns its
**transpose**, because OU-III propagates a world-to-body quaternion. The two are
not interchangeable; swapping them silently inverts every rotation in the
filter. `convention-test.cpp` asserts the transpose relation directly, so the
day someone "fixes" one of them the build goes red.

`ou_detail::quat_from_delta_theta(dtheta)`, on the other hand, *does* equal
`Exp(dtheta)` — also asserted.

## 2. The group

```
X = (R, X_w, B_b)     R   = R_bw in SO(3)
                      X_w = [v p S a_w] in R^{3x4}, world frame
                      B_b = [b_g b_a]   in R^{3x2}, body frame

(R1,X1,B1) o (R2,X2,B2) = ( R1 R2, X1 + R1 X2, B2 + R2^T B1 )
identity                = (I, 0, 0)
(R,X,B)^-1              = (R^T, -R^T X, -R B)
```

The asymmetry in the last slot is the whole point. World vectors are carried by
`R`; body vectors by `R^T`, from the other side. A direct product
`SE_4(3) x R^6` cannot express it, and that is exactly the difference between
Level 1 and Level 2 in `Kalman3D_Wave_TFG`'s `two_frame_bias` flag.

Right-invariant error in local coordinates:

```
dR     = Rhat R^T = Exp(phi)
rho_x  = xhat - dR x       for each world column x
beta_b = R (bhat - b)      for each body column b
```

Note that `rho_v` is **not** `vhat - v`, and `rho_p` is not `phat - p`. The
attitude error and the translational errors share one frame, which is what makes
a correction rotate `v, p, S, a_w` coherently instead of nudging the quaternion
and leaving the rest behind.

### Tangent layout

Deliberately matches `Kalman3D_Wave_OU_III.h:40-54`, so the 12x12 world block
stays contiguous and the scalar-loop propagation and Joseph-update block
structure port over at the same offsets:

```
xi = [ phi | beta_g | rho_v rho_p rho_S rho_a | beta_a ]
       0     3        6     9     12    15      18
```

`TwoFrameGroup<T, NW, NB>` is templated on the column counts; `NB = 1` drops the
accelerometer bias and gives 18 states, matching OU-III's
`with_accel_bias = false`. The offsets are `static_assert`ed in the test.

### Group exponential

```
E_R = Exp(phi)
E_X = J_l(phi)  [rho_v rho_p rho_S rho_a]
E_B = J_l(-phi) [beta_g beta_a]
```

**The `J_l(-phi)` on the bias block is derived, not chosen.** The group law
induces the one-parameter ODE

```
b'(s) = beta - [phi]x b(s)
```

whose solution at `s = 1` is `Exp(-phi) J_l(phi) beta`, and
`Exp(-phi) J_l(phi) = J_r(phi) = J_l(-phi)`.

Because Lie++ has no oracle for four world vectors plus six body-frame biases,
this block is pinned by the **one-parameter subgroup property**,
`Exp_G(s1 xi) o Exp_G(s2 xi) = Exp_G((s1+s2) xi)`, which only holds when the
exponential genuinely integrates that ODE. Changing `J_l(-phi)` to `J_l(+phi)`
fails 182 checks.

### Retraction

```
Xhat+ = Exp_G(xi) o Xhat-

Rhat+ = E_R Rhat-
Xhat+ = E_R Xhat- + E_X
Bhat+ = Bhat- + (Rhat-)^T E_B
```

### Adjoint

`X Exp(xi) X^-1 = Exp(Ad_X xi)`, with

```
phi    -> R phi
rho_i  -> R rho_i + [x_i]x R phi
beta_j -> R beta_j + R [b_j]x phi
```

Checked at finite `xi` (so a first-order-only adjoint fails) and against central
finite differences.

## 3. How the group math is validated

Two independent references, because finite differences alone share their
assumptions with the code under test:

1. **Lie++** (`third_party/Lie-plusplus`, `group::SEn3<double,4>`) — a separately
   derived implementation by the authors of the equivariant filtering literature
   this filter follows. Covers the world block to `1e-12` on compose, inverse,
   `exp`, `log`, and `Adjoint`. See `third_party/Lie-plusplus/VENDORED.md` for
   the convention correspondence and why it is test-only.
2. **Structural identities** for what Lie++ cannot reach: the one-parameter
   subgroup property, exp/log round trips in both directions, the exact adjoint
   relation, `J_l = int_0^1 Exp(u phi) du` by Simpson quadrature, and
   `retract`/`local` as an inverse pair.

The whole battery runs at `double` (tight tolerances, the real check) and at
`float` (loose, proving the firmware instantiation compiles and keeps its
small-angle branches).

### On testing the small-angle branches

Every branch in `SO3Jacobians.h` switches at `|phi| = 1e-2`. It is tempting to
test continuity by comparing `f(seam - eps)` against `f(seam + eps)` — **this
does not work.** Those are different arguments, so the difference is dominated
by the genuine slope of `f`, and the check fails spuriously at any useful
tolerance.

What has to hold is that *both* branches are accurate at the seam; if each is
within `tol` there, the function is continuous across it to within `2 tol`. So
the test compares each branch against a closed form carried in `long double`,
which has the headroom to absorb the cancellation in `(1-cos t)/t^2` and
`1 - (t/2)cot(t/2)` at `t = 1e-2` and still leave about 14 digits.

### Mutation results

Every check in these suites has been confirmed to fail when the code is
deliberately broken. A test that cannot fail is not evidence.

| Mutation | `lie_group-test` | `convention-test` |
|---|---|---|
| bias block uses `J_l(+phi)` | 182 fail | 50 fail |
| adjoint bias coupling sign flipped | 50 fail | passes (does not exercise `Adjoint`) |
| `retract` applied on the right | 200 fail | 350 fail |
| `compose` body block not transposed | 529 fail | 150 fail |
| `Exp` handedness flipped | 253 fail | 516 fail |
| `J_l` 0.1% coefficient error | 228 fail | 50 fail |

## 4. Measurement updates

### Sign convention

One convention, stated once, used everywhere:

```
xi     is the error OF THE ESTIMATE:  eta = Xhat X^-1 = Exp_G(xi),  P = E[xi xi^T]
r      is the innovation, measured minus predicted
r      ~ H xi + noise                      (plus, not minus)
K      = P H^T (H P H^T + R)^-1
xihat  = K r                               the estimated error
Xhat+  = Exp_G(-xihat) o Xhat              remove it
```

The correction is injected **negated**, because `xi` is the error rather than
the correction. Half the sign bugs in filters of this shape come from mixing
the two conventions -- deriving `H` against `r ~ -H xi` and then injecting
`+K r`, which cancels out only if both mistakes are made together, and leaves
a filter that works until someone fixes one of them. Here `H` is literally
`d r / d xi`, which is what `tfg_jacobians-test` measures by finite
differences, so there is nothing left to get backwards.

### The three residuals

| Update | Innovation | `H` (phi, db_g, rho_v, rho_p, rho_S, rho_a, db_a) | Noise |
|---|---|---|---|
| Accel | `r_a = Rhat(f_m - bhat_a) - (ahat_w - g)` | `[ [g]x, 0, 0, 0, 0, -I, -Rhat ]` | `Rhat R_a^body Rhat^T` |
| Mag | `r_m = Rhat m_m - B_w` | `[ -[B_w]x, 0 ... ]` | `Rhat R_m^body Rhat^T` |
| Integral | `r_S = -Shat` | `[ [Shat]x, 0, 0, 0, -I, 0, 0 ]` | `R_S`, world frame |

Neither `H_a` nor `H_m` contains the estimated attitude or the estimated wave
acceleration. They are built from **fixed** references -- gravity and the world
magnetic vector. That is the structural gain over OU-III, whose
`J_att = -[f_cog_b]x` is rebuilt from the current specific-force estimate on
every step. `tfg_jacobians-test` asserts it directly: rotate the estimate by a
large angle, move `a_w`, and the Jacobians must not change at all.

### The accelerometer-bias column is `-Rhat` at Level 1

Not `-I`. At Level 1 `delta_ba` is an additive **body-frame** error, so it has
to be rotated into the world frame the residual lives in. It collapses to `-I`
only at Level 2, where `beta_a = R(bhat_a - b_a)` is already world-referred.

This is a concrete statement of what the two-frame bias geometry buys, and it
is a trap: writing `-I` at Level 1 is wrong by a full rotation and still looks
right at zero attitude. The test pins the value *and* asserts the test
attitude is far enough from identity for the two to be distinguishable.

### What is approximated, exactly

Differentiating the accelerometer residual exactly gives

```
d r_a / d phi = -[Rhat (f_m - bhat_a)]x + [ahat_w]x = [g]x - [r_a]x
```

because `Rhat(f_m - bhat_a) = ahat_w - g + r_a` by definition. So `[g]x` is the
exact derivative **only where the residual vanishes**, and the neglected term
is exactly `-[r_a]x` -- first order in the innovation, vanishing as the filter
converges. The magnetometer has the same structure: `-[B_w]x - [r_m]x`.

The term is dropped deliberately. Keeping it would put the measurement and the
estimated attitude back into the Jacobian and forfeit the one property this
formulation exists for. It is the standard invariant-EKF choice.

The tests treat it accordingly, and this is worth copying elsewhere: finite
differences are checked at a **consistent state** where the identity is exact,
and a **separate** test asserts that the gap at an inconsistent state equals
`-[r]x` to `2e-8`. Asserting only that the gap is "small" would pass equally
well if the gap were some other small thing, and would stop discriminating the
moment a genuine error of similar size appeared. Pinning its exact structure
means any additional error shows up immediately, however small.

### Covariance

Rank-3 Joseph form, `P+ = (I-KH)P(I-KH)^T + K R K^T`, in the invariant tangent
frame. There is deliberately **no** MEKF-style covariance reset afterwards:
`ou_detail::apply_left_error_reset`'s `G = I + 0.5[dtheta]x` exists because
OU-III zeroes its error state after injection and has to transport `P` into the
new linearization point. Here the update is already expressed in the
post-update tangent frame, so applying that transport on top would
double-correct.

### Two testing notes

**Finite-difference step size.** Central differences have error
`~ h^2 f''' + eps_machine |f| / h`. Reaching for a smaller `h` past the optimum
makes things worse, not better: at `h = 1e-9` the rounding term is
`2.2e-16/1e-9 = 2.2e-7`, which swamps the quantity being measured. `1e-6` is
the right order here.

**Asserting that something did not happen.** Use raw state comparison, not
`error_from`. `error_from` computes `Log(Rhat R_ref^T)`, and the floating-point
product of a rotation with its own transpose is `I +- 1e-17` rather than
exactly `I`, so it reports about `1e-16` of "error" between two bitwise
identical states. For inertness checks that is noise; compare `R`, `X` and `B`
directly.

## 5. The two bias geometries

`two_frame_bias` selects how the bias error is parameterized. Level 2 is the
default; Level 1 is kept as a build flag because it isolates group-retraction
defects from bias-geometry defects, and because it is the ablation that
answers whether the bias geometry earns its keep.

| | Level 1 | Level 2 |
|---|---|---|
| error | `delta_b = bhat - b`, body frame | `beta = R(bhat - b)`, world frame |
| injection | `Bhat + delta_b` | `Bhat + R_pre^T J_l(-phi) beta` |
| `Phi_bias,bias` | `I` | `Exp(w_world h)` |
| `Phi_phi,bias` | `-Rbar h` | `-J_l(w_world h) h` |
| `H_a` bias column | `-Rhat` | `-I` |

### They are the same filter in different coordinates

Related by the exact linear map `beta = Rhat delta_b`:

```
T     = blkdiag(I, Rhat, I, I, I, I, Rhat)
xi_L2 = T xi_L1,   P_L2 = T P_L1 T^T,   H_L1 = H_L2 T
```

Under **pure propagation this is exact**, and that is the strongest test in
the suite. Level 1 was already verified against finite differences, so
asserting the relation holds validates the entire Level-2 derivation -- `Phi`,
`Q`, and the bias self-dynamics -- against an independently checked reference.
Finite-differencing Level 2 against itself could never do that; it would only
confirm the code computes its own formula.

Under measurement updates the two agree only to first order, because the
Level-2 injection carries a `J_l(-phi)` the Level-1 one does not.

### The conjugation is time-varying

`Phi_L2 = T(h) Phi_L1 T(0)^-1`, with `T` evaluated at *different* times on
each side. That is precisely why `Phi_beta,beta = Exp(w_world h)` where
Level 1 has `I`: the bias error is carried along by the vehicle's rotation
instead of sitting still. Derived directly,

```
betadot = Rdot delta_b + R delta_bdot = [R w]x beta - R w_b
```

This term has no Level-1 analogue and is the substantive new content of the
two-frame geometry.

### Where the noise lands

The bias random walk is body-frame, and `beta` carries it to the end of the
step:

```
Q_beta = int_0^h Exp(w_world (h-s)) Rhat(s) Q_b Rhat(s)^T Exp(w_world (h-s))^T ds
```

Because `Exp(w_world (h-s)) Rhat(s) = Rhat(h)` identically, the integrand is
the constant `Rhat(h) Q_b Rhat(h)^T`. So the Level-2 bias blocks of `Q` use
the rotation at the **end** of the step, not its start. Using `Rhat(0)` is
wrong at order `|w| h` -- the same order, and for the same reason, as using
`Rhat` rather than `Rbar` in `Phi`. This was caught by the equivalence test,
not by inspection.

### A note on making a test discriminate

The equivalence test originally used a physical gyro-bias random walk of
`1e-10`, three orders below the accelerometer's. An error in the *gyro* bias
block's frame convention then landed just under tolerance, and mutation
testing showed the test passing on a deliberately broken build. The
equivalence under test is algebraic, not physical, so the fix is to inflate
the input until both blocks are exercised equally. Worth remembering
generally: a tolerance calibrated against the largest term in a matrix says
nothing about the smallest.

## 6. r_S is a standard deviation

`Kalman3D_Wave_TFG::set_RS_noise` takes a **standard deviation** and squares it
internally, matching `Kalman3D_Wave_OU_III::set_RS_noise` exactly. The tuning
law produces `r_S` in those units — the article calls it "the standard-deviation
scale whose square forms the integral pseudo-measurement covariance".

This is recorded because getting it wrong was expensive to diagnose. An earlier
version of the orchestrator passed `r_S` into a setter that took a variance, so
the applied covariance was `r_S` rather than `r_S^2` — roughly a factor of ten
of over-tight regularization at a typical operating point.

**It did not look like a units error.** The `S = 0` constraint is a high-pass,
and over-tightening a high-pass *overshoots* rather than attenuates, so the
symptom was a heave amplitude 37% too large and a 46% RMS error that read like a
modelling problem. Two things made it findable:

1. **Decomposing the error instead of staring at the total.** Fitting `p_z` to a
   sinusoid at the known frequency split the 46% into gain, phase, DC and
   residual. The residual was 0.002 m — so it was not noise, not drift, and not
   the missing detrender. Almost all of it was gain, which points at the
   regularizer rather than at the integration chain.
2. **Comparing against OU-III on identical input.** That is also where the first
   comparison misled: the harness passed `r_S^2` to OU-III, which squared it
   again, making OU-III's regularizer far too weak and its gain look better
   (1.05) than it is. Once both filters were given a standard deviation, OU-III
   measured *worse* than this filter (1.23 vs 1.19).

`tfg_orchestrator-test::test_rs_units_are_a_standard_deviation` guards it now.

### What remains after the fix

About 21% RMS on a 0.15 Hz monochromatic sea, and it is the regularizer's
high-pass overshoot, shared with OU-III:

| | 0.08 Hz | 0.10 Hz | 0.15 Hz | 0.20 Hz | 0.30 Hz |
|---|---|---|---|---|---|
| TFG gain | 1.54 | 1.39 | 1.19 | 1.11 | 1.04 |
| OU-III gain | 1.62 | 1.46 | 1.23 | 1.13 | 1.05 |

Monotonic in frequency, and equally monotonic in `r_S` (1.39 at `r_S = 3` down
to 1.05 at `r_S = 100`). TFG is better than OU-III at every point measured.

None of that is a performance claim. The `R_S_coeff = 0.35` law was calibrated
against broadband JONSWAP and PM spectra, and a pure sinusoid concentrates all
its energy at one frequency, so this synthetic case is harsher than what the
tuning was designed for. The real measurement belongs in the paired study
against recorded waves, not here.

## 7. Yaw must be anchored at magnetic lock, not merely referenced

`initialize_from_acc` levels the filter but leaves yaw at zero, because gravity
says nothing about heading. The magnetic reference is captured later, once the
attitude has settled.

Capturing it as `B_w = Rhat m` at that moment is wrong, and wrong in a way that
is easy to miss. It bakes the arbitrary startup yaw into the reference, and the
filter then holds yaw consistent with *its own starting frame* rather than with
the field — a constant gauge offset. Every other channel looks healthy while
yaw carries a fixed error.

Measured on a JONSWAP H1.5 record: **12.2 degrees RMS yaw**, against roll 0.58
and pitch 0.17. Anchoring first — rotate about world z until the measured
field's horizontal component points north, *then* capture — gives **2.37
degrees**, in line with OU-III's yaw gate as it then stood, with every other
channel unchanged. OU-III's gate came down to 1.068 degrees on the strength of its continuous
hard-iron correction; this filter now carries that correction too (section 9),
and its own gate is 1.536.

Capturing the canonical `(B_north, 0, B_down)` is also what makes `H_m` a
genuinely fixed reference rather than a frame-dependent one, which is the
property section 4 relies on.

This is worth stating because the symptom points away from the cause. A large
yaw error with clean roll and pitch reads as a magnetometer or observability
problem; it was neither. The diagnostic that separates them is whether the
error is a constant offset or a wander — a gauge offset is constant, and a real
heading problem is not.

## 8. Adaptation policy, and why each part is not a free parameter

The tuner does not know which filter it is driving, so everything OU-III
measured about it carries over unchanged. This filter had drifted from that
policy in six ways; all six are now aligned, and each is guarded by a test in
`tfg_orchestrator-test`.

### Exogeneity is a timing property, not only a signal choice

Feeding the tuner from the complementary observer keeps its *inputs*
independent of the filter. Necessary, not sufficient. If the schedule smoothed
during step `k` were also applied during step `k`, the covariance the MEKF uses
at `k` would depend on `y_k` — the filter would be weighting a measurement
against a covariance that measurement had just moved.

So the smoother runs every sample, the cadence tick marks a candidate, and the
commit happens at the top of the *next* `update()`, before `y_{k+1}` reaches the
MEKF. `test_schedule_is_exogenous_to_the_current_sample` delivers one violent
sample and asserts the applied `r_S` does not move within that same update.

### Cadence

Committing every sample is not "more adaptive" — it couples the schedule to the
measurement stream 200 times a second. The commit runs on a 0.1 s tick
(OU-III's `ADAPT_EVERY_SECS`). The EMA still sees every sample, so the
trajectory is unchanged apart from being held piecewise constant between ticks.
`test_adaptation_is_cadenced` counts distinct committed values over 10 s and
requires ≤ 120 against 2000 samples.

### The `r_S` floor is 0.15, not 0.4

On OU-III the 0.4 floor was never a safety limit — it was the binding
constraint on every low-motion sea. The schedule asks for 0.24 m·s at the
calibrated H_s = 0.27 m point, so the floor clipped it, and a full sweep of
`c_tau`, `c_sigma` and `c_R` left those scenarios constant to three decimals.
Dropping it recovered −8.3% on the worst sea, −9.5% on PM-Stokes 0.27 m, −2.5%
on the eight-sea mean, and cut the near-still H_s = 0.05 m case from 27.0% to
17.6% of H_s. This filter had inherited the old 0.4 and the same clipping.

### The S=0 cadence is self-similar in tau, and `r_S` must follow it

One pseudo update has covariance `r_S^2`; at one update per `T_S` seconds the
continuous-equivalent information rate goes as `1/(r_S^2 T_S)`. Scaling `T_S`
with tau while holding `r_S` fixed therefore changes the regularization
strength with sea state as a side effect. The cadence is
`T_S = (0.015/1.1)·tau`, clamped to [5 ms, 250 ms], and the filter input is
renormalized by `sqrt(T_0/T_S)` to hold the information rate. That turns the
base `sigma_aw·tau^3` schedule into an effective `sigma_aw·tau^(5/2)` one.

The renormalized value is deliberately **not** re-clamped. The smallest-sea
point sits on the base floor and must be allowed below it once `T_S > T_0`, or
the clipping the previous item removes comes straight back.

### The sigma channel is a wave-band quantity

`sigma_aw` is meant to be the acceleration the OU process has to carry. A
broadband variance is not that: it also holds the sensor floor, the engine band,
and whatever drifts below the sea. This filter measured it on the raw
complementary-observer output and subtracted a constant floor.

It now goes through `AdaptiveWaveBandPass`, whose two corners are fixed
multiples of the tuning frequency (0.5 and 4.0), so away from the absolute
clamps the transfer shape is fixed in `f/f_tune` — the condition the JONSWAP
similarity argument for `sigma_aw` needs. The pre-band noise floor is referred
through that same band's white-noise variance gain rather than subtracted raw,
because the band rejects part of it and the band moves.

Measured over the eight records, at the deployed coefficients: vertical RMS
−1.1% mean and −1.9% worst, 3D −0.5% mean and −1.2% worst, yaw −4.0% mean and
−14.3% worst, accelerometer bias −5.7% mean. `TFG_WAVE_BAND=0` restores the
broadband channel for ablation.

### The a_w marginal is re-aligned on the adaptation cadence

`commitTune_()` writes a new stationary scale into the OU process covariance,
which is how the schedule reaches the *increments*. It does not touch the
marginal the filter is already holding, so after a sea-state change the
posterior `a_w` variance can sit well under the scale the schedule now claims,
and the accelerometer update is weighted against a confidence nothing supports.

`Kalman3D_Wave_TFG::synchronize_aw_covariance_to_stationary()` stages a
correction that is applied inside the next propagation, and adds only the
positive-semidefinite part of the shortfall — so it can lift a marginal that has
fallen behind and can never remove what the measurements established. It runs on
its own clock rather than inside `adaptMekf_`, so the fixed and frozen ablation
arms get the identical policy and stay matched controls.

This is the smallest of the six by some distance: −0.1% on 3D mean, +0.8% on
accelerometer bias mean, nothing measurable on vertical. It is here for the
consistency argument rather than for the number. `TFG_AW_COV_SYNC=0` ablates it.

### The coefficients were re-fitted afterwards

Every item above changes what the tuner coefficients multiply, so leaving them
at their previous values would have measured the features against a schedule
fitted for a filter that no longer exists. The sweeps are one-at-a-time over the
eight records:

| coefficient | was | now | shape of the sweep |
|---|---|---|---|
| `tau_coeff` | 1.0 | 1.0 | clean minimum; 0.92 costs +2.4% vertical, 1.10 costs +2.5% |
| `sigma_coeff` | 1.0 | 1.0 | 0.9 buys 0.4% vertical at 4% on 3D and 4% on bias |
| `R_S_coeff` | 0.35 | 0.28 | band-passed sigma reads lower, so `r_S ~ sigma tau^3` over-regularized; flat over 0.26–0.30 |
| `S_factor` | 1.87 | 1.20 | −5% bias, −0.05° pitch, no vertical cost |
| `R_S_xy_factor` | 1.0 | 1.15 | better on vertical, 3D, yaw and bias at once; stops being free above ~1.2 |

All five were later re-swept paired over three IMU seed triplets per record
rather than on one realization, which confirmed four of them and took
`S_factor` the rest of the way down, from 1.20 to **1.00** — the horizontal
acceleration magnitude the eight records actually contain is 0.997 of the
vertical one. The table above is what this change set;
[`docs/tfg-adaptation-refit.md`](tfg-adaptation-refit.md) carries the current
values, the measurement and the yaw sentinel it moved.

Four more were swept and left at OU-III's values because the eight-record result
is flat across the range: `adapt_tau_sec` (1.0–3.0 moves vertical RMS by 0.01%),
`adapt_RS_mult` (2–4, likewise), the pre-band noise floor (0.06–0.20 — the band
refers it, so the schedule barely notices), and the sigma band ratios themselves
(0.5/4.0 is at or within noise of the best of five shapes tried). Fitting those
to eight records would have been fitting noise.

## 9. Startup policy: MahonyProxy, and the magnetic acquisition that makes it work

`StartupInitPolicy` selects who solves the attitude the filter starts from.
**`MahonyProxy` is now the default**, as it is in both OU families. `StagedMekf`
remains as the ablation against the previous behaviour.

Under `StagedMekf` the MEKF runs from the first sample and levels itself while
degraded (linear block off, bias frozen, `Racc` inflated), and the magnetic
reference is captured in that attitude's frame. Under `MahonyProxy` the private
Mahony observer owns levelling and magnetic acquisition, and the MEKF is seeded
once and goes straight to Live.

### Why the proxy policy used to be worse here

The seeded tilt error was absorbed by the horizontal accelerometer bias, and did
not come back.

Measured on jonswap H1.5 at the time, the accelerometer-bias error reached 744%
of the true bias on the Y axis — 0.34 m/s², which is `atan(0.34/9.81) = 2.0°` of
tilt. The scored pitch error on the same record was 2.08°. They are the same
quantity. Tilt and accelerometer bias are only weakly separable, so a seed that
is two degrees off is indistinguishable from a bias and gets parked there.

OU-III's own measurements say where a 2° seed comes from: the observer's mean
tilt error is 2.76° on the worst record at 7–20 s, 0.85° at 40 s and 0.74° at
120 s. Handing off before it converges seeds exactly that error.

Three defects in the handoff itself were fixed earlier — each the same mistake,
an over-confident prior on something the seed does not know:

| defect | effect |
|---|---|
| `v, p, S, a_w` seeded to zero with the constructor's `1e-4` prior | asserts position known to 1 cm at a random wave phase |
| `b_a` seeded to zero with σ = 0.01 m/s² against a true bias up to 0.059 | filter refuses the correction it needs |
| accelerometer-bias learning unlocked immediately at handoff | the seed transient is absorbed into bias before tilt can settle |

What remained missing was the acquisition itself.

### The magnetic acquisition is windowed and two-stage

The previous version learned the world reference from **one** magnetometer
sample, in the proxy's tilt frame, the first time the gravity gate was happy —
one reading of a noisy vector at one phase of one wave, frozen for the run.

`MagAutoTuner` now averages the field in the proxy's yaw-stripped tilt frame
until a window closes. The window is held in **seconds**, not samples, because
in waves the tilt error is periodic and what the window has to buy is whole wave
periods: 128 samples at a 25 Hz mag ODR is 5.1 s, short enough to lock in the
phase it started on rather than cancel it. Default `mag_min_window_sec = 15`.

That reference is still provisional. It is re-learned at 90 s over a 30 s window
(`mag_refine_*`), once the observer has converged — a device that reports no
heading for 105 s is not a device, so the first stage buys availability and the
second buys accuracy. Both stages average in the **proxy's** tilt frame, never
the MEKF's: the MEKF has been steering to the reference the refinement exists to
replace, so its tilt carries that reference's error and the refinement would be
self-confirming. The refined stage replaces both the reference vector and the
heading gauge; the yaw write is a one-time step through
`Kalman3D_Wave_TFG::set_attitude_yaw_absolute`, which rewrites heading only and
leaves the world columns alone — this is a state correction, not a gauge change,
so `v`, `p` and `S` must not rotate with it.

Accelerometer-bias learning is held shut until the refinement lands, for the
same reason it is held through the seed transient.

The refinement's own start time and window length were swept (45/60/90/150 s and
15/30/60 s) and the eight-record metrics move **non-monotonically** across both
— `accb` worst goes 134 / 167 / 137 across the three window lengths. That is the
phase sensitivity of a single one-shot re-acquisition, not a trend, so the
values stay at OU-III's rather than being fitted to eight records.

### Hard iron is tracked continuously

The MEKF has no magnetometer-bias state. A body-fixed offset is therefore
heading error one-for-one against the horizontal field: an offset component
`b_h` perpendicular to north on a field with horizontal component `B_h` is
`atan(b_h/B_h)` of yaw. **This was most of this filter's standing yaw error.**

`ContinuousMagHardIronEstimator` never closes its accumulation. It is driven by
the proxy's tilt quaternion and the **raw** magnetometer — never the corrected
stream, or it would be shown data with its own answer already subtracted — so no
loop closes through the MEKF and the exogeneity argument of section 8 is
untouched. The applied offset and the magnetic reference move together, out of
the same statistics, so the filter is never subtracting one offset while
steering to a reference that belongs to another. The correction is a change of
measurement-model parameters only: no attitude state is written, and the
magnetometer update walks yaw to the corrected heading over its own time
constant.

Measured over the eight records, ablating it alone: **yaw 2.35° → 1.08° mean and
3.30° → 1.53° worst**. Nothing else moves by more than 3%, and accelerometer
bias moves the other way by 2.8% — the offset was buying part of its heading
error back through the bias state. It is the single
largest effect in this whole change, and it is what took TFG's yaw gate from
2.938 to 1.536.

### Evidence

Both simulators, byte-identical records, last 900 s, eight JONSWAP and PM-Stokes
records. Against the previous TFG default (`StagedMekf`, one-sample reference,
no hard iron, broadband sigma, old coefficients):

| channel | before | after | change |
|---|---|---|---|
| vertical RMS %H_s, mean | 4.638 | 4.369 | −5.8% |
| vertical RMS %H_s, worst | 5.209 | 4.778 | −8.3% |
| 3D RMS %, mean | 21.04 | 19.74 | −6.2% |
| 3D RMS %, worst | 25.78 | 21.03 | −18.4% |
| yaw RMS deg, mean | 2.149 | 1.077 | −49.9% |
| yaw RMS deg, worst | 2.923 | 1.528 | −47.7% |
| roll RMS deg, mean | 0.698 | 0.308 | −55.9% |
| pitch RMS deg, mean | 0.306 | 0.340 | **+11.3%** |
| acc-bias RMS % of true, mean | 144.1 | 92.7 | −35.7% |
| acc-bias RMS % of true, worst | 398.2 | 166.7 | −58.1% |

Pitch is the one channel that moves the wrong way, and it should be stated
rather than averaged away: the loss is concentrated on the two large PM-Stokes
records (0.55° → 0.73° at H4.0, 0.34° → 0.60° at H8.5) and comes with the
refinement, which halves roll on the same set. Combined attitude RMS
`sqrt(roll² + pitch²)` goes 0.762° → 0.459°, so the trade is net favourable, but
it is a trade.

The result holds off the tuning set. Repeating the whole comparison under two
further sensor-noise and bias realizations (`W3D_SEED=20260813` and
`W3D_SEED=777`) improves **every** channel in the table on **both** — worst-case
vertical −30.7% and −9.9%, worst-case accelerometer bias −68.7% and −71.2% —
pitch included, at 0.558° → 0.369° and 0.674° → 0.345°. So the pitch regression
above is a property of the default seed's realization rather than of the change,
which is also why it is not worth tuning against.

Against OU-III on the same records, TFG's yaw is no longer attributable to a
missing feature (1.53 worst against OU-III's 1.068 gate), 3D on JONSWAP sits at
OU-III's own bar and PM-Stokes below it, and accelerometer bias remains the
outlier at 167% against 98.4 — still above 100% of the true bias on two of the
eight records, and still unattributed.

## 10. What this does not claim

The bias-free, frozen-parameter skeleton is group-affine. The complete
estimator is **not** shown to be, and must not be described as an InEKF or as
having log-linear error propagation. It carries estimated body biases, OU decay,
an online `tau`, a tuner driven by estimated signals, staged startup, yaw
locking, resets, and a synthetic `S = 0` observation. Verifying the group-affine
property of the fixed-parameter continuous dynamics is a prerequisite for any
such claim, and has not been done.

The accurate name is **right-invariant two-frame Lie-group error-state EKF**.

A Lie-group formulation also does not cure the physical observability problem
between tilt, horizontal acceleration, and accelerometer bias. It can represent
that uncertainty more consistently; it cannot create information the sensors do
not provide.
