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

## 6. What this does not claim

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
