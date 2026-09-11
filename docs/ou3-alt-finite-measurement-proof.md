# ALT finite physical measurement identity

## Scope and role in the master inequality

The state is `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, with
`c=2 tan(theta/2) u` and `E=R_true R_hat^T`. These are **finite** errors,
not perturbations. This note proves the accepted, finite, real-arithmetic,
zero-lever measurement map used in a fixed-trajectory dissipation argument.
The actual estimator supplies covariance, geometry, tuning, gain masks and
measurement conditioning at every event. Their dependence on the trajectory
is retained; there is no comparison with a second, frozen-gain estimator.

The controlling expression is `V(z_N)-rho V(z_0)-physical/bias supply` over the
complete literal word. A product of pointwise Jacobians is not an identity
for `z_N`. In particular the physical S event does not fix zero error when
`S_phys != 0`. Before evaluating that expression, the word must either compose
finite maps with all reference offsets, or have a proved anchored mean-value
representation on the entire required domain. The construction below gives
finite measurement descriptors directly, without a Taylor remainder or an
unproved radial-segment bridge. It is a necessary map-construction step, not a
local contraction certificate or permission to refine an infeasible metric.

The code is `tools/stability/ou3_alt_contraction/finite_measurement_graph.py`.
Its rational interface evaluates the identities exactly. Universal coefficient
checks expand all monomials, rather than treating a finite collection of
sampled states as proof. It intentionally returns no source-coverage or
stability PASS.

## 1. Exact residual secant

Write `C=[c]x`, `D=4+c^T c`, and

```
U = 4 I + 2 C + c c^T,
T = U/D,
E = I + (4 C + 2 C^2)/D.
```

Using `C^2=c c^T-(c^T c)I` and `c^T C=0` gives

```
(I-C/2) U = D I,       E-I = T C.
```

Consequently, for every vector f, **exactly**

```
(E-I)f = -T [f]x c.
```

For an accelerometer event, let `f_hat=R_hat(a_hat-g)` and let `nu_acc` be the
actual conditioned sensor/model residual relative to `R_true(a_true-g)+beta`.
Temperature-centering differences belong in that actual model residual unless
the physical beta definition already includes them. Then

```
r_acc = -T [f_hat]x c + E R_hat e_aw + e_ba + nu_acc.
```

The coefficient of e_ba is I in BOTH H18 and A21: holding the estimate does not
remove a physical bias error from the measurement. For magnetometer data,

```
r_mag = -T [m_hat]x c + nu_mag,
nu_mag = measured_conditioned_m - E m_hat.
```

Thus a learned-reference error or hard-iron/de-heel discrepancy is not silently
set to zero. Its bound and source ancestry remain required in the full theorem.
For every due S event,

```
r_S = e_S - S_phys.
```

All three cases have the exact form `r=Hbar(z,xi)z+nu`. Hbar is a **secant
factor**, not the derivative used to compute the gain. The physical S value
retains its single Live-origin primitive; nu is not permission to create an
independent bounded S input. All source relations stay in the master.

## 2. The actual masked gain and held covariance

Use the full physical 21-state covariance even in held mode. In symmetric real
arithmetic and the accepted unrepaired positive-innovation branch, let H be the
shipping linearization, H_g its column-masked version and D_h the row mask:

```
Sigma = H P H^T + R,
N = D_h P H_g^T,
Sigma q = r,             delta = N q.
```

In A21 both masks are identity. In H18 the BA columns of H_g and BA rows of N
are zero, but H used for Sigma still contains BA. This is exactly the distinction
in `measurement_update_acc_only`; it is not an unmasked information-form
replacement. The raw Joseph covariance polynomial remains

```
P_after = P - K N^T - N K^T + K Sigma K^T.
```

**Held covariance reduction.** Conditional on the existing held invariant
`P=diag_block(P18,B0)` with zero BA cross blocks, the accelerometer operands are

```
Sigma_H = H18 P18 H18^T + R_acc + B0,
N_H = [P18 H18^T; 0].
```

For arbitrary P18, this follows by block multiplication, not by fitting a
trajectory. Substitution into the actual Joseph polynomial leaves BA cross
blocks zero and B0 unchanged; the motion covariance is exactly the reduced
18-state covariance update with **effective noise R_acc+B0**. Magnetometer and
S residual Jacobians have no BA column and need no such addition. Prediction
and reset preserve the block invariant when their held BA blocks are identity,
BA process covariance is zero, and their other blocks do not mix BA. The
queued a_w covariance addition has no BA rows or columns.

This is a covariance computation reduction, NOT an 18-state error/storage
reduction. The finite residual still contains e_ba and the storage still has
24 coordinates and all cross terms. A kernel using bare R_acc in the reduced
H18 gain needs this bridge and corrected covariance ancestry before its word
can stand for shipping. Merely fixing the final gain does not fix earlier P.

The existing release-floor lemma uses the same B0 at enable, so
`max(B0_ii,B0_ii)=B0_ii` exactly. Its interval implementation must enclose the
real seed product and apply max by endpoint selection; an additional numerical
widening of max is not a mathematical change of the covariance fixed point.

## 3. Exact finite quaternion injection

At a normal measurement boundary the nominal local attitude-error slots are
zero (initialization and every successful injection clear them; prediction
propagates the reference quaternion rather than assigning those slots).
Let `d=delta[0:3]`. Shipping constructs the normalization of `(w,k d)`, with
`u=d^T d` and the following real-arithmetic branches:

```
||d|| < 1/100:
  w = 1-u/8+u^2/384,
  k = 1/2-u/48+u^2/3840;
otherwise:
  w = cos(||d||/2),
  k = sin(||d||/2)/||d||.
```

These are the two branches of `quat_from_delta_theta`. They remain hard graph
relations. Neither the small polynomial nor its derivative replaces the other
branch. The point `||d||=1/100` belongs to the axis-angle branch.

The estimator injects on the left, so the physical relative attitude satisfies
`E_after=E Q(d)^T`. Multiplying homogeneous quaternions `(2,c)` and `(w,-k d)`
gives

```
W = 2w + k c^T d,
V = w c - 2k d - k c cross d,
c_after = 2V/W.
```

For every c,d,w,k with nonzero W, normalization cancels and subtraction yields

```
W(c_after-c) + k U d = 0,
c_after = c - L(c,d) d,
L(c,d) = k U/W.
```

This is a polynomial identity after clearing the denominator, valid for all
values of its eight scalar indeterminates, and therefore for both deployed
quaternion branches. No small-angle replacement `c_after=c-d` is used. W=0 is
not accepted by the finite chart evaluator. A uniform separation from that
pole is still an every-prefix/chart-retention obligation, not a supplied fact.

## 4. Projection and the full finite descriptor

Lift N into a thin 24-by-3 factor B: its attitude rows are `L N_theta`, its
remaining estimated-error rows are the corresponding rows of N, and its beta
rows are zero. The pre-projection error is exactly `z_pre=z-Bq`.

Set `b_pre=beta-e_ba_pre`. The radial projection has the exact graph

```
alpha=1                    when ||b_pre|| <= R_b,
0<alpha<1,
alpha^2 ||b_pre||^2=R_b^2   when ||b_pre|| > R_b.
```

Then `e_ba_after=alpha e_ba_pre+(1-alpha)beta`. Let P_alpha be identity except
for precisely those BA/error-truth blocks. The full mean map is

```
z_after = P_alpha (z-Bq),
Sigma q = Hbar z + nu,
d = N_theta q.
```

The same projection runs after injection in held mode. Taking alpha=1 there
requires the actual held estimate already be inside the projection ball; it
must not be assumed from a covariance-consistency claim.

The implementation uses `chi=(z24,q3,nu3,h)` with h=1. Its descriptor equalities
include both `N_theta q=d h` and `z[0:3]=c h`, so the finite-correction parameters
cannot be detached from the state and innovation variables. Quaternion and
projection branch graphs, P/H/R ancestry, conditioned-sensor relations and the
physical source constraints remain additional hard constraints. Matrices
computed at one point alone do not certify that family.

## 5. Exact storage contribution using thin factors

For any symmetric common 24-by-24 M, put `M_alpha=P_alpha^T M P_alpha`.
Expanding the finite map, with no omitted cross terms, gives

```
V_after-V_before =
  z^T(M_alpha-M)z
  -2 z^T M_alpha B q
  +q^T (B^T M_alpha B) q.
```

M_alpha is computed by the three affected bias columns and rows, not a dense
24-state congruence. The measurement work uses `M_alpha B` (24-by-3) and
`B^T M_alpha B` (3-by-3). This is exact thin-factor arithmetic on full storage;
it does not freeze any gain, delete true bias, or assert low rank of the whole
word. A later derivative of this graph must still retain derivatives of B and
all its endogenous coefficient factors; rank three of its value correction does
not bound the rank of that derivative. The contribution can be summed inside
the one physical word with its hard constraints before deciding the sign. There is no claim that individual
measurements must be contractive.

## Remaining bridge and non-promotion

These identities close the finite accepted-measurement algebra. They do not
close continuous physical attitude-versus-sampled-gyro prediction forcing,
all covariance/frontend successors and repairs, the asynchronous branch family,
or a complete finite endpoint/prefix certificate. The current
`physical_word.compose_endpoint_lineage` produces Jacobian products; those
must not be submitted as the finite map in a storage search. The executable
finite-master guard rejects that representation before the coarse endpoint
or metric calculation starts.

Binary32 row-wise LDLT solves, raw innovation asymmetry, scalar accumulation,
normalization, projection and covariance updates need their actual numerical
residual relations. The real formulas above do not bound them. BIAS family
admission, a useful ultimate bound, startup capture and all final ALT gates
remain separate. Original P2/P3/P4/P5 machinery is unchanged.
