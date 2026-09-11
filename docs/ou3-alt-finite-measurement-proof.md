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

## 6. Finite physical prediction, not a sampled shadow truth

`finite_physical_prediction.py` supplies the physical predictor. Let q_t and
q_h be the physical and estimated world-to-body quaternions at the same time,
and c their Cayley error. Let

```
Q_phys = q_t(t+h) conjugate(q_t(t)),
Q_nom  = actual shipping quaternion increment from the conditioned gyro,
(W,V)  = Q_phys (2,c) conjugate(Q_nom).
```

Nonzero projective scaling is harmless. Exactly, `W c_next = 2 V` when W is
nonzero. To prove this, substitute `q_t ~ (2,c) q_h` in the two attitude updates.
Their relative quaternion is

```
q_t_next conjugate(q_h_next)
  ~ Q_phys (2,c) q_h conjugate(q_h) conjugate(Q_nom)
  = ||q_h||^2 Q_phys (2,c) conjugate(Q_nom).
```

The common scalar cancels in Cayley coordinates. The test expands this equality
for all 15 independent polynomial variables (c, q_h, Q_phys, Q_nom), rather than
sampling attitudes. Thus the relation is finite, not its derivative, and its
value at zero error is not silently assumed to vanish.

The physical increment is NOT a second call to the estimator's quaternion
propagator. For a source-defined sampled physical increment Q_sample, retain

```
D_omega = Q_phys conjugate(Q_sample),
Q_phys ~ D_omega Q_sample.
```

The quaternion product ordering is checked by a second all-coefficient identity.
D_omega retains continuous angular integration, conditioning/frame changes and
sampling/model defects. It is not set to identity. If a proposed decomposition
of D_omega contains an unknown estimation error, that term must remain in the
state graph; it is not an independently bounded disturbance. The source/gyro
relation must also retain `omega_hat=omega_sample+e_bg+n_g`. The earlier
`finite_prediction_deployed_step` evaluates a sampled shadow and, by itself,
does not provide this physical increment.

For one translation axis, define the SAME physical integrals

```
J0 = integral_0^h a(s) ds,
J1 = integral_0^h (h-s) a(s) ds,
J2 = integral_0^h (h-s)^2 a(s)/2 ds.
```

Writing va, pa, Sa and alpha for the actual shipping mean coefficients gives

```
e_v_next  = e_v + va e_aw + J0 - va a0,
e_p_next  = e_p + h e_v + pa e_aw + J1 - pa a0,
e_S_next  = e_S + h e_p + h^2 e_v/2 + Sa e_aw + J2 - Sa a0,
e_aw_next = alpha e_aw + a1 - alpha a0.
```

These follow by subtracting the actual mean prediction from the physical
kinematic equations. All 17 polynomial indeterminates are checked coefficient
by coefficient. In particular, pa and Sa are the **shipping** coefficients,
including `safe_phi_A_coeffs`' small-step polynomial branch, not silently the
ideal exponential integrals. The physical discrepancy remains in the formula.
The moments are not independent inputs: the global corrected COMPLETE-BRMM
relation, including `D_S<=1100 m*s`, must still constrain them jointly.

`PhysicalSegment` enforces positive segment duration and one persistent
Live origin, consecutive physical v/p/S moment equations, and the shared true
bias recurrence. The joint24 predictor additionally checks that beta in z is
that predecessor's beta. It retains

```
e_bg_next = e_bg + bg_true_next - bg_true,
beta_next = phi_true beta + u_b,
e_ba_next = phi_hat e_ba + (phi_true-phi_hat) beta + u_b.
```

H18 sets phi_hat=1; A21 uses its actual estimator factor. BIAS2's phi_true=1
endpoint is retained. Neither this algebra nor a consistent finite segment
qualifies BIAS0, BIAS1, BIAS2 or an all-time physical history. In particular a
finite moment-consistency check alone cannot enforce the all-time primitive
bound. The same beta goes on to the measurement and radial projection graph.

## 7. Branch-correct rank-three covariance arithmetic

On the symmetric real accepted branch, the actual inverse-free gain relation
`K Sigma=N` implies, for arbitrary N,

```
P-K N'-N K'+K Sigma K' = P-K N'.
```

This is a direct polynomial substitution. It does not require `N=P H'`, an
optimal gain, or an information-form replacement, so the actual masked H18
numerator is allowed. `solved_joseph_covariance` checks the solve relation
exactly before applying the thin expression. An all-coefficient entrywise
proof and exact full-21 covariance regressions cover held and active masks.
With a numerical solve defect E=K Sigma-N, the omitted term would be E K'; the
routine rejects nonzero E instead of pretending to enclose deployment arithmetic.

`benchmark_finite_rank3.py` compares time and peak Python allocation, and checks
zero **exact rational** entrywise difference using the same P/N/Sigma/K operands.
It does not rectangularize any coefficient. Source-uniform enclosure width,
subdivision count and certificate margin have not been benchmarked; no claim
about their improvement follows from a faster rational evaluation.

## 8. What the composed implementation checks establish

`shipping_finite_identity.py` builds both the unchanged runtime and a passive
observation overlay. Erasing only the added observation statements recovers
the shipping source modulo whitespace. The two builds must produce bit-identical
**recorded** sample state, including the frontend audit coordinates. Neither
build overwrites a tracked shipping header or sets a synthetic gain, covariance,
mode, scheduler, or startup state.

One analytical oscillatory sensor history drives the real startup and runtime.
The checker obtains a 600-step H18 window, a 600-step A21 window and a 600-step
window across the actual release guard. It verifies physical prediction, finite
sensor/S residuals, masked N and full Sigma, inverse-free solves, Joseph,
quaternion injection, covariance reset, same-beta projection, queued a_w floor,
S due/not-due service credit, and the observed H18->A21 setter. Consecutive
observed core states/covariances are checked against each other; the Live S
origin is never restarted. The first H18 window checks fresh centered e_S=0.

These are implementation correspondence regressions, NOT universal word or
stability evidence. Binary64 evaluations of formulas are compared to host
binary32 values using explicitly labelled regression tolerances; those
numbers are not rigorous roundoff bounds. The finite polynomial identities
above are the algebraic theorem contribution. The finite traces neither
qualify the analytic source family nor fit a metric.

## Remaining bridge and non-promotion

The complete source-uniform finite word remains OPEN. The source-independent
algebra can be composed by substitution once **every** event's runtime operands
and guard have been bound to the same physical predecessor and successor. A
collection of locally correct events or a recorded execution is not that
binding. In particular the finite relation must still materialize all
frontend/tuner/covariance successors, angular-defect source constraints,
innovation-repair/rejection and watchdog branches, asynchronous magnetic
history, and every legal H18/A21 edge. Finite core observations are not a proof
of complete literal runtime-prefix coverage.

No high-precision feasibility diagnostic, common metric search, worst rho,
retention radius or startup capture bound is authorized by these regressions.
`physical_word.compose_endpoint_lineage` remains a rejected pointwise-Jacobian
representation. Do not submit either that representation or these captured
operands to a source-uniform storage search.

Binary32 row-wise LDLT solves, innovation asymmetry, accumulation, normalization,
projection, floors, scheduler and frontend arithmetic still need their numerical
residual relations. BIAS admission, a useful ultimate bound and startup capture
remain separate. All final ALT gates remain false. Original P2/P3/P4/P5
machinery is unchanged.
