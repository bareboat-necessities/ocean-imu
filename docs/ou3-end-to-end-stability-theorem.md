# OU-III end-to-end shipping stability theorem

## The target

The goal is not "P4 passes". It is the end-to-end statement

> for every motion admitted by the declared BRMM source, every admitted
> BIAS0/BIAS1/BIAS2 history and every admitted sensor disturbance, the shipping
> implementation started from the declared startup uncertainty reaches a
> certified Live basin in finite time through its actual Mahony/proxy startup,
> and then stays there with practical ISS.

Symbolically,

```
declared startup uncertainty
  --[P5 capture, finite time]-->  R_P4  --[P4 invariance]-->  practical ISS
```

P3 supplies the covariance/observability properties along the admitted
execution; it is not the stability statement. P4 proves "inside this basin the
shipping filter is stable". P5 proves "the actual Mahony/proxy startup gets you
into that basin in finite time". Only together do they answer the question.

The executable holder of the composition is
`tools/stability/ou3_end_to_end_stability_gate.py`. It is fail-closed:
`END_TO_END_STABILITY_PASS` is the exact conjunction of `P5_CAPTURE_PASS` and
`P4_INVARIANCE_PASS`, every open obligation carries a declared failure class,
and no stage may be skipped.

## The basin is maximised, never shrunk

Every bit of width certified in P4 is width the capture stage does not have to
achieve. A small basin does not make the theorem easier; it moves the work into
P5 and can demand unrealistic settling times. So the design rule is

> maximise the certified P4 basin subject to the actual nonlinear filter
> remaining provably stable inside it,

and the basin is a **correlated polytope**, not a product box. Coordinate-wise
worst cases that the implementation cannot jointly generate are not part of the
theorem's obligation, but removing them requires a reachability argument, never
an assertion.

`tools/stability/ou3_p4_basin_frontier.py` materialises the trade-off. Scaling
entry ball `j` by `s_j` scales exactly its term of the retained per-prefix
subadditive sum, so the chart hypothesis is the half-space

```
sum_j s_j a_j <= B,     B = chart_cap / R_attitude - t
```

with `a_j` the retained single-ball attitude reach and `t` the template reach.
For H18 the declared product box overshoots the chart by 5.0388, and the
frontier gives

| Frontier point | attitude | gyro bias | velocity | position | integral | latent | accel bias |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| declared box | .5359 | .01 | 5 | 20 | 300 | 2.942 | .4 |
| largest uniform inflation (.19793) | .1061 | .00198 | .9896 | 3.959 | 59.38 | .5823 | .07917 |
| maximum-volume point | .1382 | .0550 | 5.625 | 7.565 | 11.03 | 8.969 | 1.600 |
| others at declared radius | .5359 | .01 | 5 | 20 | **4.566** | 2.942 | .4 |

The frontier points are not comparable coordinate-wise, which is the point: the
uniform inflation buys a 59.38 m*s integral radius by paying for it everywhere
else, while the maximum-volume point keeps velocity above its declared radius,
multiplies gyro bias by 5.5 and latent acceleration by 3.0, and still admits
11.03 m*s of integral. The volume optimum is the larger basin by measure, and
being free to choose among such points is exactly what "correlated, not
Cartesian" buys.

Two caveats. These are frozen-map binary64 diagnostics on one captured word, so
they define the CANDIDATE maximal basin the outward certificate has to close.
And the half-space models the attitude chart constraint ONLY: the accelerometer
bias radius of 1.600 at the volume optimum is chart-admissible but not
physically attainable, because the deployed radial projection keeps
`|e_b| <= R + B_true = .6252`. Every coordinate still owes its own deployed cap.

## What the deployed startup actually hands to Live

`tools/stability/ou3_p4_live_entry_reachability.py` settles the entry-set
question from the shipping sources rather than by choosing radii. Three parity
facts, all checked against `SeaStateFusionFilter_OU_III.h` and
`Kalman3D_Wave_OU_III.h`:

1. `updateFrontEnd` runs the whole front end with `drive_mekf=false`, and every
   MEKF drive call inside `updateCore_` sits inside an `if (drive_mekf)` block.
   The MEKF translational block, the latent-acceleration state and both bias
   states are therefore **never propagated before `goLive`** and stay at the
   constructor value, which `xext.setZero()` sets to zero.
2. `goLive` calls `initialize_from_attitude`, which seats attitude on the Mahony
   proxy quaternion and calls `zero_AL_cross_cov_once_()`, zeroing the entire
   attitude-block to linear-block cross covariance.
3. `enterLive_` reseats the latent OU covariance on the committed operating
   point.

The audit is over the whole fusion source, not one function. Every MEKF entry
point that writes the error state is enumerated, and only three occur outside
the guard, each with its own argument: `initialize_from_attitude` IS the entry
event; `initialize_from_acc`, reachable pre-Live through the public
passthrough, writes `xext.head<3>()` only and itself zeroes the cross; and
`measurement_update_mag_only`, reachable pre-Live through `updateMag`, has zero
linear gain rows because the constructor writes only block-diagonal covariance
seeds and the only propagator that could create `P_lin,att` is the guarded
`time_update`. Any fourth occurrence fails validation instead of being absorbed.

So at the Live entrance instant `T`,

```
e_v(T)  = -v_true(T)      e_p(T)  = -p_true(T)      e_S(T)  = -S_true(T)
e_aw(T) = -a_w_true(T)    e_bg(T) = -b_g_true(T)    e_b(T)  = -b_true(T)
P_theta,v = P_theta,p = P_theta,S = P_theta,aw = 0,   K_theta,S = 0
```

with attitude the only coordinate the startup stage has to earn. Every other
entry coordinate is the negated physical truth of ONE admitted BRMM history at
ONE instant. They are not independent, and `e_S` is the exact running integral
of `e_p` along that same trajectory.

That is why the 300 m*s integral ball was never the right object. The reachable
integral entry radius is not an estimator accumulation at all; it is the BRMM
bounded-integral-displacement primitive `S_m`. Likewise the velocity and
position entry radii are `V_m` and `P_m`, and the accelerometer-bias entry
radius is the BIAS0/1/2 true-bias envelope .2252, not the .4 projection radius.

`V_m`, `P_m` and `S_m` are declared in the physical BRMM contract's
`unfrozen_physical_constants` and are all `None`. Instantiating them from the
declared source class is the class-E qualification the theorem needs; inventing
values here would be the fitting the contract forbids.

## What the certificate admits for S_m

Two independent constraints bound `S_m`, and they agree in order of magnitude.

`tools/stability/ou3_p4_coefficient_dependency_cover.py` materialises the
adaptive dependency the theorem forbids discarding. The deployed pinnings are
exact: `T_S = clamp(pseudo_ratio*tau)`, `R_S_target = SpectralMSE(tau,sigma)`
and the applied anisotropic `R_S` is the ray `(.72r, .72r, r)`. So the reachable
target quadruple is a two-parameter surface inside a nominally four-dimensional
rectangle, and on the non-saturated branch a `T_S` strictly inside its rails
fixes `tau`, hence the tuning frequency. The clamped target box is also forward
invariant for the candidate and active schedules, because every EMA step is a
convex combination and the staged commit copies the candidate; no reachable-set
argument is needed for that. `candidate_ema` does not clamp, so a convex
combination stays inside only if it starts inside: the invariance is conditional
on initial membership, and the producer checks that premise against the deployed
initial schedule rather than leaving it implicit.

Evaluating the joint `S=0` residual scale `(T_S+dt)/(.72 R_S)` on cells where
cadence and applied covariance come from the SAME schedule gives .9709, against
1.4352 at the independent `(T_S_max, R_S_min)` corner an independent rectangle
would license: a 1.478 over-approximation in scale and 2.185 in energy.

Turning that into an admissible `S_m` requires dividing by the transverse
attitude variance, and this is where the retraction in the next section bites.
Under the retracted prior-independent cap the same-cell correction ceiling
admitted `S_m <= 7.3537 m*s`. Under the conditional cap that actually holds, the
allowed residual scale falls from 87.10 to 1.639 while the dwell term from the
declared pre-entry position envelope is 19.43, so the headroom is NEGATIVE and
the magnitude-only route admits no `S_m` whatever. The producer reports the
budget as absent rather than as a number.

The chart frontier still admits 4.566 m*s with all other coordinates at their
declared radii, and 11.03 m*s at the maximum-volume basin point, so the geometry
is not what fails. What fails is the magnitude-only Cauchy-Schwarz accounting,
and the recorded fallback is the one that survives: retain the same-cell
attitude/`S` cross-covariance direction instead of the `lambda_max` product.

## The attitude covariance cap, and a retraction

An earlier revision of this document, and of
`tools/stability/ou3_p4_blocker_falsification_classification.py`, asserted a
prior-INDEPENDENT cap on the attitude variance transverse to the specific force:
for a scalar angle measurement with gain `|f|` and noise variance `r`,

```
P^+ = P^- r / (P^- |f|^2 + r) <= r / |f|^2    for every prior P^-,
```

giving 1.18634e-3 rad^2 per axis at `R_acc = .04`, `|f| >= 5.80665`, executed at
every valid IMU sample by the declared Normal-Live invariant.

**That claim is false for the deployed filter and is retracted.** The scalar
identity holds only when attitude is the sole state in the residual.
`measurement_update_acc_only` builds `J_att = -skew(f_cog_b)` alongside a
latent-acceleration block and an accelerometer-bias block, so a transverse
attitude error and an `a_w` error produce the SAME residual; one update cannot
separate them, and the attitude MARGINAL of the full-state posterior is not
capped by anything. `tools/stability/ou3_p4_attitude_measurement_cap.py`
exhibits the refutation directly: on the deployed residual structure at a
latent-acceleration prior of 1e8 the transverse attitude marginal is 1.0291e6
rad^2, exceeding the claimed cap by 2.4743e9. The verdict is taken outside the
covariance-cancellation noise floor, and the degenerate attitude-only case is
checked to still satisfy the scalar identity, so the refutation is the residual
structure and not an arithmetic artifact.

### What holds instead

Through the variational form of the posterior,

```
c^T P^+ c = min_k [ (c - H^T k)^T P (c - H^T k) + k^T R k ],
```

take `c = (e,0,0)` with `e` a unit vector orthogonal to `f`, and choose
`k = -(f x e)/|f|^2`. Then `H_theta^T k = e` exactly and `|k| = 1/|f|`, leaving
`c - H^T k = (0, -R^T k, -k)`, so for every prior and every such `e`

```
e^T P^+_theta,theta e  <=  (sigma_a^2 + 2 lambda_max(P_(a_w,b_a))) / |f|^2 .
```

The bound is CONDITIONAL on the joint latent/bias covariance block, taken from
the certified BRMM covariance ceiling rather than assumed:
`lambda_max(P_(a_w,b_a)) <= 56.4621`, using
`lambda_max([[A,B],[B^T,D]]) <= lambda_max(A) + lambda_max(D)` for PSD blocks.
That gives 3.35035 rad^2 per axis. It is tight rather than merely true, and the
evidence is reproducible from the repo: `sweep()` in the cap module searches the
full nine-state residual with random specific force, random carried rotation and
prior condition numbers over many decades, runs on every build, finds no
violation, and attains about 91% of the bound. Cases where the variational
identity fails to reproduce `c^T P^+ c` to a relative 1e-8 are rejected rather
than counted either way, since at those condition numbers the posterior
subtraction loses the answer outright. A larger out-of-repo sweep of 188609
validated cases attains 99.75%.

### Consequence for the correction/reset blocker

The retracted cap made the same-cell `S^{-1} = (H P H^T + R)^{-1}` accounting
appear to close the accelerometer correction inside the reset utility domain.
Under the conditional cap it does not: the ceiling is 3.38698 against the reset
utility limit 3.0, and the `R^{-1}` relaxation gives 180.024.

That near miss is not luck, and it is sharp. Retaining the same-cell `S^{-1}`,

```
ceiling(P)^2 = 2 P E_acc r / (f^2 P + r),      P = transverse variance per axis,
```

is strictly increasing in `P` with supremum `sqrt(2 E_acc r / f^2) = 3.38758`.
The supremum EXCEEDS 3.0, so the route does not close for free by saturation;
but because it is finite and increasing there is an exact threshold,

```
P* = C^2 r / (2 E_acc r - C^2 f^2) = 4.31271e-3 rad^2 per axis   (C = 3.0),
```

i.e. 6.5671e-2 rad, 3.7627 deg one-sigma. So this obligation closes through
this route if and only if the transverse attitude variance is bounded by `P*`,
and the conditional cap is 776.85 times too loose. Propagating `P*` back through the
cap turns the blocker into one scalar target on the certified latent/bias
ceiling:

```
lambda_max(P_(a_w,b_a)) <= (P* f^2 - sigma_a^2)/2 = 5.27063e-2,
```

against the 56.4621 now certified: a shortfall of 1071.3. All of these are
rounded outward by the producer, so the stated `P*` and `lambda*` are safe to
aim at rather than optimistic.

That target is attributable to one block. On the certified ceiling the
accelerometer-bias half contributes 6.25e-4 and the latent-acceleration half
56.4615, so the bias block already meets the target on its own by a factor of 84
and the latent acceleration owes the entire shortfall. The open obligation is
therefore not "tighten the covariance ceiling" in general -- it is a bound on the
LATENT ACCELERATION posterior variance, which is precisely the state the
accelerometer residual cannot separate from attitude in one event. The circularity
is the finding: the attitude bound needs the `a_w` bound, and the `a_w` ceiling is
loose for the same confounding reason. Breaking it needs the window, where
attitude, `a_w` and `b_a` separate because their dynamics differ -- attitude drifts
with the gyro bias, `a_w` mean-reverts at rate `1/tau`, and `b_a` is constant. That target is the
concrete statement of what remains, and it is achievable in principle rather
than excluded.

The separate observation that the retained endpoint-referenced envelope
(3.99983e8 rad^2 per axis) exceeds the trace the reset domain admits by 6.4485e11
is unaffected by the retraction -- it compares the envelope against
`C^2/E_acc`, not against the retracted cap -- and it is why the blocker stays
class C. What the retraction removes is the claim that the one-shot measurement
event by itself repairs the envelope. It does not, and it cannot: the bound has
to come from the uniform observability/detectability machinery, where the
latent-acceleration and bias states are separated over a WINDOW instead of at a
single event. Rotation about the specific force is not accelerometer-observed at
all and still needs the asynchronous magnetometer plus gyro-bias transport.

## Obligation status

### P5 capture, from the deployed state machine

The shipping stages are `Cold -> TunerWarm -> TunerReady -> (external bootstrap
calls goLive) -> Live`.

| Obligation | Stage | Status | Class |
| --- | --- | --- | :---: |
| elapsed-time warmup | Cold -> TunerWarm | **established** (pure elapsed-time comparison, 10 s) | -- |
| tuner + usable period attained | TunerWarm -> TunerReady | open | E |
| bootstrap tilt and north acquired | TunerReady -> goLive | open | E |
| certified proxy attitude radius | goLive | open | E |
| entry state inside the certified basin | goLive -> Live | open (needs V_m, P_m, S_m) | E |
| held-bias to A21 release covered | H18 -> A21 | open | E |

### P4 invariance

| Obligation | Status |
| --- | --- |
| all three bias families closed | **established** |
| conditional binary32 additive ISS | **established** |
| motion practical ISS | open: six gate blockers, classes C and E |

No class A counterexample and no rigorous class B infeasibility exists, so
nothing here supports saying the theorem is unprovable on the declared domain.

## Honest qualifications

"For all admissible BRMM and BIAS scenarios" is only as strong as *admissible*.
If the declared source class permits a history under which the implementation
can never acquire sufficient vector information, or the startup admits unbounded
initial attitude or bias error, finite-time capture is simply impossible. The
theorem therefore needs explicit, physically defensible admission assumptions --
sensor availability, magnetic and gravity geometry, finite initial uncertainty,
scheduler operation, BRMM and BIAS bounds. Those assumptions must stay broad and
representative of the intended deployment, and the necessary excitation and
observability assumptions must be stated rather than hidden.

**END_TO_END_STABILITY_PASS = false. P4_PASS = false. P5_MAY_START = false.**
