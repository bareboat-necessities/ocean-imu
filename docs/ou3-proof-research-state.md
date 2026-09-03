# OU-III proof research state

This file is the current research ledger required by the root `AGENTS.md`. Keep it short and replace stale state rather than accumulating PR history.

## State

**REPLAN COMPLETE — NEXT EXECUTION AUTHORIZED**

Current PR research question: make canonical P3 reach one rigorous numerical PASS/FAIL verdict at the unchanged `1e-18` gate; only if P3 passes, run one non-promoting structure-exact P4 complete-word feasibility diagnostic. P5 is out of scope.

## Failure classification

The current canonical run is **not a numerical P3 theorem FAIL**. It is blocked earlier by a translation-covariance enclosure/conditioning failure:

`Loewner prediction lower lost strict SPD; split x cell`

The previously identified `x=0.01` branch-partition defect is already fixed on the live PR branch: branch clamping is now applied only to the actual boundary subcell. It was an implementation defect and is not the remaining mathematical limiter.

## Evidence

For source node 137 / gap 13, the first uncertain prediction on the widest initial small-`x` cell has diagnostic center `lambda_min ~= 3.405e-10`, absolute row-radius `eps ~= 1.142e-5`, and `eps/lambda_min ~= 3.35e4`. Existing binary subdivision improves that ratio by about 2x per level, but at depth 12 it is still about `8.26`; roughly four additional binary levels would be needed merely to make this first 5-ms prediction plausible. That is not a viable source-complete architecture.

The theorem-relevant 13-sample point covariance is much healthier. On source 137 its diagnostic `lambda_min` is about `2.02e-5 ... 2.83e-5` across the tau endpoints. Across the ten physical tau cells, using the actual low-sigma/strongest-`R_S` source in each cell, the 13-sample point floor remains strict; the longest-tau cell still has `lambda_min ~= 1.20e-6`.

Point diagnostics also show the 13-sample covariance increasing in Loewner order as `x=h/tau` increases on all ten tau cells: a 101-point grid produced positive minimum eigenvalue increments, with the smallest observed increment about `4.3e-9`. This is feasibility evidence only, not a proof of monotonicity.

Two segment-level interval experiments have now failed:

1. **Post-hoc congruence normalization of the old natural interval Riccati recursion.** Dependency has already exploded before normalization. Even with 256 tau subcells, the normalized common lower remains strongly non-SPD (worst diagnostic normalized `lambda_min - eps` about `-80`).
2. **Natural-interval derivative/monotonicity propagation.** Propagating `P(x)` and `dP/dx` as ordinary interval boxes inherits the same recursive dependency loss; subdivision improves widths but does not produce a usable derivative-SPD certificate.

A third, independent experiment confirms the same conclusion from the opposite
direction. `one_step` imposed three strict-SPD obligations per sample (after the
prediction and after each scalar measurement update), 39 per 13-sample subcell,
where the theorem makes one: P3 consumes the segment endpoint floor, certified
downstream by `common_boundary_floor` via `certified_rho(posterior) > 0`.
Removing the 38 surplus obligations collapses the adaptive split tree and drops
the segment module from 1191 s to 434 s, and the failures stop being
`split x cell` exceptions. But the propagated endpoint lower is then **not SPD**
on the same cells, while the point diagnostics above put the true 13-sample
`lambda_min` at `2.02e-5 ... 2.83e-5`. So the per-step SPD demand was never the
limiter, and the collapse architecture loses the entire floor between the true
value and the propagated lower.

Therefore the rejected mechanism is now broader than the original `C-eps I` step: **recursive natural interval covariance boxes that forget the repeated scalar source parameter at every Riccati operation are frozen as a dead end.**

Canonical P3 still has not reached translation/H/A margin calculation, so none of these diagnostic values are P3 theorem margins.

## Endpoint-only certification succeeds; the rejection was against the wrong quantity

The dead-end call on recursive interval covariance boxes was measured against
the **intermediate per-step** SPD demand, where the first uncertain prediction
has `lambda_min ~= 3.4e-10` and about four more binary levels would be needed.
The theorem does not impose that demand; `common_boundary_floor` certifies the
**segment endpoint**, where `lambda_min ~= 2e-5`, about 5.9e4 times larger.
Against the endpoint quantity the existing machinery succeeds.

Measured at 64 uniform x subcells with the intermediate gates removed, 18 of 18
tested (node, gap) pairs certify a strictly positive endpoint floor:

| node (tau,sigma,R_S) | gap 13 | gap 20 | gap 26 |
| --- | --- | --- | --- |
| 0 (0,0,0) | 64 / 9.50e-7 | 64 / 1.18e-5 | 64 / 2.00e-5 |
| 137 (1,5,7) | 64 / 8.74e-7 | 64 / 9.04e-6 | 64 / 1.58e-5 |
| 399 (4,7,9) | 64 / 7.75e-5 | 128 / 4.72e-4 | 128 / 4.93e-4 |
| 555 (6,7,5) | 64 / 8.43e-5 | 128 / 3.70e-4 | 128 / 4.05e-4 |
| 729 (9,0,9) | 64 / 8.69e-8 | 64 / 5.47e-7 | 128 / - |
| 799 (9,7,9) | 64 / 6.08e-5 | 64 / 7.14e-5 | 128 / 2.68e-4 |

Subdivision now has the scaling argument it previously lacked: halving the cell
width flipped 33/33 failing to 65/65 passing, and the endpoint margin improves
roughly linearly in width. This is a different quantity from the rejected one,
not a third refinement of it.

**The blocker is now cost, not feasibility.** `common_boundary_floor` is 800
nodes x 14 gaps = 11200 segment scans; at the observed 2 s to 33 s each the full
run extrapolates to about 103 h against a 120 minute CI budget.

The floor is provably monotone in the sigma and R_S cell index, so node
`(tau,0,0)` dominates its whole tau cell and the scan collapses to 10 nodes:

- `Q` scales as `sigma_lo^2` and the Riccati map is Loewner monotone in `Q`;
- `d/dR [P - Pe(e'Pe+R)^-1 e'P] = Pe(e'Pe+R)^-2 e'P >= 0`.

This is domination among actual reachable nodes within one tau cell, not the
Cartesian tau/sigma/R_S extrema rejected above. Measured: for tau=1 the sigma
sweep gives 8.74e-7 for indices 0..5 (identical, because the committed filter
sigma is clamped at the 0.05 floor and the code uses `sigma.lo`), 2.47e-5 at
index 6 and 6.75e-5 at index 7. Reducing 800 nodes to 10 brings the estimate to
roughly 25-77 min, which fits the budget only just, and the dominating nodes are
the slow ones.

## The S=0 pseudo-measurement is inert over a segment

Across a 23-decade override sweep at node 80 the S update is live and strong
when `R_S_z` is small -- at `1e-4` it collapses `P[2][2]` from 4.277 to 8.95e-5
and moves the floor from 2.83e-5 to 4.62e-6 -- but it saturates by `R_S_z ~ 1e6`
and the deployed value is `7.46e11`, six orders past saturation. The control
sweep confirms the accelerometer update is in its active region and is doing the
work. So this is a real regime fact, not a defect.

The reason is the zero start: after 13 samples `P[2][2] = 4.277` in `z_S = S/h^3`
units, i.e. an `S` standard deviation of `2.6e-7 m.s`, against a pseudo
measurement standard deviation of `0.108 m.s`, a factor of 4.2e5. A triple
integral accumulates little uncertainty in 65 ms.

Consequence to test, not yet a claim: P3 translation observability is the
four-`S` spread argument over a 0.765 s window inside the 3.17 s word, but
`common_boundary_floor` builds its floor from a zero start over a single 13-26
sample segment and therefore never sees that mechanism. The bound stays valid --
`rs.lo` is the strongest measurement in the cell, which is the correct choice for
a lower bound -- but the segment floor may be a lossy proxy for the true relative
injection margin, and that is a candidate explanation for `delta` landing near
1e-18.

## Measured: the margin deficit is about 9 orders, and neither enclosure nor scalarization explains it

`delta` was computed directly against the real `Sigma_upper` from
`HIST.endpoint_phase_upper`, at phase 0, for three source nodes. The target is
identical across them:

`Sigma_upper (v,p,S,a_w) = [1.119e9, 3.970e9, 7.179e9, 1.109e3]`

| node | delta with `rho*I` (global rho 8.69e-8) | delta with the true anisotropic segment floor | gain | binding group |
| --- | --- | --- | --- | --- |
| 0 | 1.890285e-31 | 1.610690e-27 | 8.5e3 | S |
| 137 | 1.890285e-31 | 1.518656e-27 | 8.0e3 | S |
| 729 | 1.890285e-31 | 1.607829e-28 | 8.5e2 | S |

The canonical gate is `1e-18`. The current isotropic floor is **13 orders**
short; repairing it to the full anisotropic floor still leaves **9 orders**.

Two consequences.

First, the isotropic `rho*I` collapse is a real defect worth about 8.5e3, but it
is not the limiter. An earlier estimate in this ledger put the gain at ~5e6 from
the anisotropy ratio `lambda_max/lambda_min ~ 4e6`; that was wrong. `delta`
comes from a matrix PSD binary search, not the diagonal bracket: for node 0 the
S per-group ratio `d_i^2 Pz_ii/upper_i` is 9.2e-24 while the certified `delta`
is 1.6e-27, about 5700 times lower.

Second, and more important: **no enclosure or scalarization work can close a 9
order gap.** Taylor models, tighter collapses, finer subdivision and the
sigma/R_S cost reduction are all worth factors of 10^3 or less against a deficit
of 10^9. They should not be pursued as a route to the gate.

### Where the deficit plausibly comes from

The S group binds everywhere, and its ratio is a tiny floor over a very large
ceiling. `Sigma_upper` for S is `7.18e9 (m.s)^2`, a standard deviation near
85000 m.s, against a declared entrance box of 300 m.s -- the ceiling is the
same-history covariance propagated over the word, some 8e4 times the entrance.
The floor is `h^6 * Pz[2][2] = 1.5625e-14 * 4.25 = 6.6e-14`.

`common_boundary_floor` discards all older process covariance and restarts every
segment from `P = 0`, which is what makes the segment bound easy to state. For a
triple-integrated OU state the accumulated floor grows steeply in horizon; over
the 635 sample word versus a 13 sample segment that is roughly
`(635/13)^7 ~ 1.4e12`, the right order to cover a 10^9 deficit. This is the same
conservatism that makes the S=0 pseudo-measurement inert over a segment.

Hypothesis to test before any further enclosure work: the zero-start
per-segment floor, not its enclosure, is what costs the margin. The falsifiable
version is to propagate a floor across a horizon long enough to contain the
four-`S` spread and recompute `delta` against the same `Sigma_upper`. If
`delta` moves by orders, the segment formulation is the limiter and the fix is
structural. If it does not, the theorem formulation itself is marginal and
should be reconsidered rather than re-enclosed.

## Physical-reality check on the delta comparison target

`Sigma_upper` is the ceiling `delta` is measured against. In physical units it
is far outside anything a marine IMU in ocean waves can reach:

| channel | std(Sigma_upper) | declared entrance box | ratio |
| --- | --- | --- | --- |
| v | 33459 m/s | 5 m/s | 6691x |
| p | 63000 m | 20 m | 3150x |
| S | 84730 m.s | 300 m.s | 282x |
| a_w | 33.3 m/s^2 (3.4 g) | 2.94 m/s^2 | 11x |

A covariance ceiling asserting 33 km/s of velocity uncertainty and 63 km of
position uncertainty describes a diverging filter, not this one. It is a
bounding looseness, not a physical claim.

Substituting the declared physical box as the ceiling moves the binding S ratio
from 9.25e-24 to 7.38e-19 and the matrix `delta` from 1.61e-27 to about
1.28e-22. So the 9 order deficit decomposes as roughly **5 orders from the loose
ceiling and 4 orders structural**, the latter being the segment/word horizon
mismatch: a 65 ms zero-start segment accumulates 2.6e-7 m.s of `S` uncertainty,
which against any realistic `S` ceiling is intrinsically near 1e-18 at best.

`Sigma_upper` may not simply be replaced by the entrance box -- `delta` is a
relative Riccati injection margin and the ceiling must stay a valid upper bound
on the filter covariance, not on the error. But tightening that bound so it
reflects covariance the filter can actually reach is worth about 5 orders,
which exceeds every enclosure avenue combined and is not currently being
pursued by any recorded plan.

Conditions that are already correctly excluded, so they are not the problem:
impacts and slams out of normal-Live scope, lever arm disabled, vibration guard
dormant and transparent, non-gravitational CoG acceleration 4.0 m/s^2 = 0.41 g.

Entrance realism: the declared attitude entrance is 45 deg and the canonical P4
gate hard-requires 0.8 rad = 45.8 deg, while the operating domain's own search
list is `[30, 25, 20, 15]` deg. A vessel starting on flat water has about 1 deg
of tilt error; the 45 deg figure is really a heading assumption, and the domain
separately declares a 10 deg internal heading gauge error. The frozen gate
therefore forces the widest and least realistic sector against the domain's own
stated intent.

## Whole-chain physical-realism audit

Checked every declared bound against a marine IMU in ocean waves, and traced
what each costs. `Sigma_upper` is reproduced exactly by
`BASE.translation_upper` at the full declared tuner ranges, so that function is
the source.

**Correctly excluded already, not the problem:** impacts and slams out of
normal-Live scope, lever arm disabled, vibration guard dormant and transparent,
non-gravitational CoG acceleration 4.0 m/s^2 = 0.41 g.

**Rotation rate is realistic.** `body_rate_norm_upper_deg_s = 30`. A 30 degree
roll at a 6 s period gives a 31 deg/s peak, so the bound sits at genuine
heavy-sea motion rather than being inflated.

**The tuner cells are unphysical but nearly free.** `sigma_aw` reaches
6.0 m/s^2 = 0.61 g and `tau` falls to 0.333 s, i.e. 0.61 g of process noise
decorrelating in a third of a second. That is vibration, which the domain
excludes elsewhere, and it exceeds the domain's own 0.41 g acceleration
envelope. Cost of restricting to a physics-consistent `tau >= 3 s`,
`sigma <= 4 m/s^2`: `a_w` improves 4.04x, and v, p, S improve by 1.002x. So the
unphysical corner is not what inflates the translation ceiling.

**`R_S` upper drives the ceiling, but tightening it is not enough.**
`Sigma_upper` scales as `R_S,max^2` through the four-`S` inversion:

| rs.hi | Sigma_upper v | std |
| --- | --- | --- |
| 0.15 | 2.086e6 | 1444 m/s |
| 100 | 7.190e7 | 8480 m/s |
| 400 (declared) | 1.119e9 | 33459 m/s |

The declared 400 costs 538x on v, p and S. But even the tightest possible S
measurement leaves a 1444 m/s velocity ceiling, still three orders beyond any
boat. **The ill-conditioning is intrinsic to inverting the triple integrator
v -> p -> S over a roughly 0.8 s observation span**, not a consequence of any
tuner or domain choice.

### Deficit accounting on the binding S channel

| lever | worth | delta after |
| --- | --- | --- |
| current | - | 1.89e-31 |
| repair the isotropic `rho*I` collapse | 8.5e3 | 1.6e-27 |
| tightest `R_S` | 538 | ~8.7e-25 |
| both | 4.6e6 | ~4e-25 |
| **canonical gate** | | **1e-18** |

Every realism lever pulled together still leaves about **7 orders**. Physical
realism accounts for roughly 3 of the 9 order deficit; the remainder is
formulation. Do not spend further effort on domain realism as a route to the
gate.

The two structural items that remain are the intrinsically ill-conditioned
four-`S` triple-integrator observability route, and the zero-start 65 ms segment
floor compared against a whole-word ceiling.

## Master P3 quantity

The new backend is relevant only if it produces the complete-segment matrix lower used by canonical P3:

`P_segment(x) >= L_segment > 0`.

That lower enters the existing translation gate directly through

`D_h L_segment D_h^T - delta * Sigma_upper > 0`,

which is what `_certified_delta` tests before the H/A precision join. Improving unrelated one-step quantities is not useful unless it sharpens this complete-segment comparison.

## What the failures invalidate

- A useful absolute point Loewner lower after every uncertain 5-ms prediction.
- Blindly increasing `MAX_ADAPTIVE_X_DEPTH`.
- Post-hoc normalization of a covariance interval after natural Riccati dependency has already exploded.
- Natural-interval `P,dP/dx` recursion as a monotonicity proof engine.

## What the failures do not invalidate

- The `1e-18` canonical P3 usefulness threshold.
- The same-history P2-V1 source language.
- Existence of a rigorous complete 13--26-sample translation floor.
- H=18/A=21 full-state joining, which has not yet received a translation margin from this backend.
- P4 feasibility or infeasibility, because P4 remains blocked by canonical P3.

## Critic pass and alternatives

Assume the interval covariance recursion architecture is wrong. Its strongest defect is loss of the fact that the entire segment depends on **one repeated scalar** `x=h/tau`. Ordinary interval arithmetic replaces that one-dimensional curve by a new Cartesian matrix box after every operation.

Qualitatively different alternatives:

1. **Univariate centered Taylor-model / polynomial enclosure of the complete 13-sample Riccati map.** Preserve the same scalar symbol through all prediction and measurement operations and bound only a final remainder.
2. **Verified polynomial collocation/Bernstein or Chebyshev enclosure of each complete-segment matrix entry.** Again preserve one-dimensional dependence rather than recursively hulling matrices.
3. **Analytic complete-segment Gramian/Riccati lower in a different information/covariance representation.** Avoid intervalizing the covariance recursion itself. Note before attempting this: the natural candidate `P_k >= (G_c^-1 + G_o)^-1`, with `G_c` the reachability and `G_o` the observability Gramian, is **false** for segments started from `P_0 = 0` -- 223/352 violations on integrator-chain systems. The information-form derivation needs a finite `Y_0 = P_0^-1`, which a zero-start segment does not provide. A correct analytic lower has to come from the closed-loop reachability sum, whose transitions depend on the gains.
4. **Different Lyapunov representation** if the complete-segment covariance map remains unsuitable.

**Selection:** pursue (1) once. It directly attacks the failed dependency mechanism, and the point map is smooth with generalized source-137 ratio `P(x)` versus the low-`x` endpoint close to `1 ... 1.08`, so a centered model should be proving a relative statement near unity rather than recovering five orders of lost SPD margin. If this Taylor-model attempt fails after one mathematically motivated refinement, freeze it and move to (3), not another interval subdivision variant.

## DEAD_ENDS

- **REJECTED: endpoint-only 800-node P2 ancestry.** It loses staged/committed path memory.
- **REJECTED: independent Cartesian `tau/sigma/R_S` extrema.** They destroy source-history correlation.
- **REJECTED: recursive absolute entrywise Loewner point lower plus deeper subdivision.** Tractable depth is insufficient.
- **REJECTED: one-step congruence normalization as the primary architecture.** It improves conditioning but still solves an unnecessary one-step property.
- **REJECTED: post-hoc complete-segment normalization of natural interval Riccati boxes.** Dependency has already exploded before normalization.
- **REJECTED: natural-interval derivative/monotonicity recursion.** It shares the same dependency loss.
- **REJECTED for #471: additional P4 micro-certificates before complete-word feasibility.**
- **PARKED: rigorous H=18/A=21 complete-word dissipation producer.** Written and
  unit-tested at commit `4d68493` (`tools/ou3_p4_complete_word_dissipation.py`),
  then removed: a rigorous certificate before the non-promoting `rho_w`
  diagnostic is the wrong order. Resurrect only after the diagnostic reports
  `rho_w` clearly below 1.

## Open contract inconsistency

`tools/ou3_p4_canonical_gate.py` requires a candidate's `outer_angle_rad` to
equal `0.8` exactly, bound to the Cayley and remainder artifacts, while the
operating domain declares a P4 certificate search over
`p4_complete_word_full_attitude_candidate_deg = [30, 25, 20, 15]`. If the
0.8-rad formulation is later abandoned for a narrower cell, the gate rejects
every candidate on that list. Resolve before relying on the search list.

## Verified theorem steps

Randomised verification against `ou3_interval.symmetric_positive_definite_ldlt`.
All hold, so none of these are suspect as sources of the P3 blockage.

| Step | Result |
| --- | --- |
| `\|\|(R-I)v\|\|^2 = 4/(4+q^2) \|\|[c]x v\|\|^2` | exact, max err 4.7e-13 |
| `\|\|eta\|\| = sin(theta/2) \|\|h\|\|` | exact, max err 6.3e-13 |
| `\|\|R-I-[c]x\|\| <= (3/4) q^2` | 0 violations |
| `J <= n blockdiag(J_ii)` for PSD J | 0 violations |
| `K_theta S K_theta' <= P_theta_theta` | 0/1200 |
| `Phi_s' Sigma_s^-1 Phi_s <= Sigma_0^-1` | 0/1200 |
| `(1-3d/8)^2 <= 1-d/2` | holds for all `d <= 16/9` |

Prefix nonexpansiveness needs `Phi_s` invertible. It is, via the Joseph form and
invertible OU prediction, but the source-path document does not state it.

## Available relaxations, not yet warranted

These enlarge the certified P4 funnel `W_*` only; none changes
`rho = 1 - delta/2`, which is set entirely by the P3 margin. Do not spend
effort on them before a real `delta` exists.

- `lambda_max(Sigma) <= sum_g U_g` instead of `n_g max_g U_g`. With `U_S ~ 9e4`
  dominating this is 9.04e4 against 6.3e5, about 7x.
- Attitude corrections need only the attitude marginal,
  `|dtheta| <= sqrt(U_theta/lambda_min(R)) |y|`, so `sqrt(0.25)` replaces
  `sqrt(9e4)` -- about 600x on the correction gain.
- The sharp Cayley remainder `q^2(q+2)/(4+q^2)` is uniformly `0.805x` the
  `(3/4)q^2` bound.
- The `4` in `B_m = 4 N_op sqrt(m_+) C / m_-` comes from a `W_s <= 4 W_0`
  bootstrap; prefix nonexpansiveness gives `W_s <= (1+d/8)^2 W_0 <= 1.27 W_0`.

A rejected route may be resurrected only after recording the new mathematical fact that invalidates its rejection.

## Retained facts

- The scalar covariance-form measurement update `P - Pe(e'Pe+R)^-1e'P` requires
  only `e'Pe+R > 0`, never strict SPD, and is Loewner monotone on **all**
  symmetric arguments with that denominator positive. For `H >= 0`,
  `x'(dU)x = x'Hx - 2(x'He)(x'Pe)/d + (x'Pe)^2(e'He)/d^2` is a quadratic in
  `x'Pe` whose discriminant is `4((x'He)^2 - (x'Hx)(e'He))/d^2 <= 0` by
  Cauchy-Schwarz in the `H` semi-inner product, so `dU >= 0`; `d` is affine in
  `P` so it stays positive between ordered arguments. Any backend may therefore
  carry singular or indefinite intermediates without invalidating the order.
- Canonical P3 useful gate remains exactly `1e-18`.
- `OU3_P2_CORRELATED_STAGE_TRANSFER_V1` and same-history source correlation are retained.
- Both H=18 and A=21 remain required.
- Zero/disabled lever arm and dormant/transparent vibration-guard branch remain the current proof scope.
- No replay fitting, operating-domain shrink for PASS, or deployed-filter change is permitted.
- Existing exact Joseph, co-rotated accelerometer, reset, and finite-angle identities remain structural facts only; they do not authorize P4 promotion.
- P4 is blocked until canonical P3 passes; P5 is blocked until canonical P4 strictly contracts.

## Next falsifiable experiment

Build a **source-137 / gap-13 diagnostic-only univariate centered Taylor model** for the complete segment in `x=h/tau`.

Requirements:

- carry the single normalized scalar source coordinate symbolically through all 13 predictions and both scalar measurement updates per sample;
- use the existing validated transition/process series and outward remainder bounds;
- use a verified reciprocal expansion for innovation denominators rather than natural interval division;
- compare the final matrix against a high-precision point reference by congruence/generalized eigenvalue;
- report polynomial order, remainder norm, certified relative `alpha`, limiting direction, and number of structural branch splits.

**Predicted success criterion:** point diagnostics suggest the exact generalized relative floor is near 1 across source 137, so the first rigorous model should certify a clearly positive `alpha` with only the required `x=0.01` analytic-branch split and modest model order. If it cannot do that, one refinement of model order/remainder formulation is allowed; a second failure triggers architecture review and the analytic segment Gramian/Riccati alternative.

Before any new P4 proof producer, first obtain the canonical P3 numeric verdict. If P3 passes, the only P4 work allowed in #471 is the non-promoting high-precision complete-word ratio diagnostic `rho_w = V_after(F_w(x)) / V_before(x)` required by the PR scope.
