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

## 4. What this does not claim

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
