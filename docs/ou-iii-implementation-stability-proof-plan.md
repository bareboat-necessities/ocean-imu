# OU-III implementation stability proof-completion plan

## Objective

Prove stability of the **actual deployed OU-III implementation**, including the initialization sequence before the MEKF enters Live mode. The target is not an existential local result of the form “there exists a sufficiently small neighborhood.” The target is a numerical, machine-checked certificate for the source-realizable hybrid implementation.

The proof target is:

> Under an explicit configured-runtime and marine operating envelope, the deployed startup reset/Mahony observer reaches a certified Live handoff set; the exact held-bias Live mode (H), the exact held-to-active transition, the exact active-bias Live mode (A), magnetic lock/refinement, tilt reset/re-lock/cooldown, and periodic covariance synchronization remain inside source-node Lyapunov funnels; the exact nonlinear source-word maps enter a certified invariant inner funnel in finite time and thereafter satisfy the practical-ISS theorem. The stochastic layer must separately satisfy the paper's localized mean-square/concentration theorem.

“Actual implementation” means `SeaStateFusionFilter_OU_III` together with `Kalman3D_Wave_OU_III`, using the shipping update order, gains, adaptation law, clamps, covariance operations, reset maps, and startup state machine.

## Existing evidence retained from `main`

The existing tuning, eight-sea replay, multi-seed validation, exact executed-map replay, RMS gates, and evidence/provenance machinery stay in place. They remain important for performance regression, falsification, and locating active constraints.

They are **not** reimplemented in this PR and are **not** substituted for source-complete theorem bounds. A replay counterexample can invalidate a proposed bound; replay success alone cannot promote a theorem claim.

## Proof scope and unavoidable qualifications

The final claim is conditional only where the physics makes an unconditional theorem impossible:

1. **Configured runtime timing.** The public API accepts arbitrary positive finite `dt`; until code enforces a bounded interval, the numerical theorem is scoped to the source-defined configured scheduler.
2. **Heading persistent excitation.** Full-heading convergence requires recurring accepted, non-collinear gravity/magnetic (or equivalent heading) information. Permanent rejection/collinearity is unobservable. Gravity-only operation is treated on the yaw quotient.
3. **Physical startup/motion/noise envelope.** Specific force, initial gyro bias, heading-reference error, and stochastic noise must be bounded by explicit numerical theorem inputs. These are deployment assumptions, never fitted from the eight trajectories.
4. **SO(3) topology.** The smooth observer cannot have a globally asymptotically stable attitude equilibrium. However, the implemented first-accelerometer reset discards prior attitude and maps every stored prior attitude into a physical hemisphere whenever the non-gravitational specific force is smaller than gravity. The proof is therefore global with respect to the discarded pre-startup estimate and regional/almost-global with respect to the physical measurement ambiguity described in the paper.

## Non-negotiable proof rules

- No proof-specific filter, retuning, Schmidt restriction, fixed gain, disabled cross-gain, or altered covariance.
- Keep the complete `S=0` Kalman gain, including `S -> attitude` coupling.
- Bind each proof primitive to source and invalidate the certificate when source changes.
- Use outward-rounded/validated arithmetic for every continuous-source bound used for promotion.
- Cover accepted/rejected measurement branches that satisfy the theorem hypotheses, not only favorable replay words.
- Keep H (18-state) and A (21-state) fixed-dimensional words separate; H→A is an explicit dimension-changing hybrid jump.
- Include quaternion injection and left-error covariance reset in every applicable source word.
- Initialization is part of the final theorem composition.
- The final result must contain numerical margins and finite capture counts/times. `exists epsilon > 0` is not a completion criterion.

## Execution plan

### P0 — Freeze the exact implemented hybrid system

Create/extend a source-derived manifest containing:

- H/A state dimensions and order;
- prediction and process-noise construction;
- `S=0`, accelerometer, and asynchronous magnetometer update order;
- Joseph covariance updates and complete Kalman gains;
- quaternion injection and left-error reset;
- `tau`, `sigma_aw`, `R_S`, pseudo-update cadence and source clamps;
- startup Mahony gains, gravity gate, aligned-branch test, hold times, rate veto and timeout;
- `goLive` state/covariance initialization and bias-held entry;
- H→A bias release;
- magnetic lock/refinement/re-gauge;
- tilt reset/re-lock/cooldown;
- periodic `a_w` covariance synchronization;
- configured runtime sampling contract.

**PASS:** the proof manifest is regenerated from the source tree and semantic parity tests pin every operation used in the theorem.

### P1 — Numerical startup/reset certificate before Live

Instantiate the paper's Mahony/reset theorem with the deployed gains `2 k_P = 0.2`, `2 k_I = 0.02` and validated arithmetic.

Prove numerically:

1. the first accepted accelerometer reset is independent of stored prior attitude and maps into the desired gravity hemisphere for the declared `|a_f| < g` envelope;
2. positive startup Lyapunov/chart margins for the declared initial gyro-bias and disturbance bounds;
3. the normal gravity-gate tilt bound from the implemented `0.075` sine gate plus the declared measured-gravity error;
4. the finite normal-gate comparison time `T_Q` when its robust disturbance floor permits it;
5. the implemented timeout path at 150 s, including the aligned-branch constraint;
6. the exact `goLive` map into the H-mode source coordinates, including seeded covariance, linear-block enable, current OU/R_S state, and held accelerometer bias;
7. full-heading and yaw-quotient startup branches separately.

**PASS:** every source-reachable accepted startup/handoff branch has an explicit finite numerical handoff family in H-mode coordinates. Containment of that family in a Live capture funnel is a separate P5 obligation and is not assumed by P1.

### P2 — Complete source-word language

Make the normal-Live source language finite-window complete relative to explicit theorem hypotheses.

- Supply an explicit finite vector-PE recurrence window as a deployment theorem input.
- Preserve arbitrary accepted/rejected branches between required PE packets.
- Combine it with the rigorous pseudo-measurement firing-gap bound.
- Tile every fixed-mode normal-Live execution by bounded H/A words.

**PASS:** every normal-Live execution satisfying the declared PE/timing hypotheses belongs to the machine-defined language; no favorable-subset selection is possible.

### P3 — Validated H/A Riccati/information-word certificate — COMPLETE

P3 is now closed by two machine-checked layers.

First, the outward-rounded direct matrix backend certifies the source-uniform H/A generalized endpoint inequality with the full coupled `[v,p,S,a_w]` translation block. Four `S=0` firings provide the complete translation observability qualification, with the validated integer spread search choosing the strongest rigorous admissible four-firing information bound. The three-firing `[v,p,S]` construction plus stable `a_w` is used only to sharpen a covariance upper bound; it cannot by itself satisfy the P3 observability qualification. The H and A endpoint margins are both strictly positive, and the arbitrary former `1e-18` threshold is not a theorem condition.

Second, `tools/ou3_p3_word_algebra.py` binds the endpoint comparison to every fixed-dimensional normal-Live covariance operation in the shipping source. For each prefix it proves the exact decomposition

`P_s = Phi_s P_0 Phi_s^T + Omega_s`, with `Omega_s >= 0`.

Prediction, accepted Joseph corrections, rejected/not-due identity branches, the left-error covariance reset, and the periodic PSD `a_w` covariance increment all have the common affine-PSD form

`P+ = A P A^T + B`, with `B >= 0`.

Thus for every `0 < delta < 1`, an established information margin propagates exactly through any subsequent source branch:

`Omega+ - delta P+ = A (Omega - delta P) A^T + (1-delta) B >= 0`.

The implemented left-error reset has attitude block `G = I + 0.5 [dtheta]_x` and

`det(G) = 1 + ||dtheta||^2 / 4 >= 1`,

so reset congruence is nonsingular for every finite injected correction and does not require an additional small-angle hypothesis at P3. The same covariance decomposition gives, by Schur complement,

`Phi_s^T P_s^-1 Phi_s <= P_0^-1`,

hence the exact source-uniform worst prefix information gain is `1.0` rather than a replay estimate.

The final producer `tools/ou3_explicit_information_word_certificate.py` emits `P3_IMPLEMENTATION_WORD_CERTIFICATE=PASS` only when the direct H/A interval inequalities, optimized four-S qualification, current source-word language, exact Live-operation algebra, reset nonsingularity, branch coverage, and unit prefix-information bound all validate together.

**PASS:** strict source-uniform H and A information contraction with explicit positive endpoint margins and an exact finite prefix information-gain bound. P4 consumes this information-word geometry directly through the exact Cayley lift below.

### P4 — Exact nonlinear SO(3) word certificate — COMPLETE

P4 is closed by `tools/ou3_p4_nonlinear_word_certificate.py` and the dependent `exact Cayley nonlinear H A source-word P4` proof-fast job. The sole quantitative nonlinear metric is the exact Cayley lift of the P3 source information geometry,

`c(R) = 2 tan(theta/2) u`,

`W_g(R,xi) = s_m [c(R);xi]^T Sigma_KF(g)^-1 [c(R);xi]`,

where one positive normalization `s_m` is shared by every source node in a fixed-dimensional mode. This preserves all attitude-linear information cross terms and does not change generalized contraction ratios or physical level sets.

The validated word map follows the shipping order exactly: prediction; a due `S=0` correction inside `time_update()` before accelerometer correction; asynchronous magnetometer handling; and immediate quaternion injection plus left-error covariance reset after each accepted `S`, accelerometer, or magnetometer correction. The complete `S -> attitude` gain remains present. The deployed normalized polynomial quaternion branch is enclosed directly, rather than replaced by a linearized/exponential correction.

P3's exact prefix information gain upper bound of `1.0` transports source-uniform quadratic nonlinear defects through arbitrary admissible accepted/rejected/not-due placements, so P4 covers the complete source branch family without enumerating an exponential list of rejection strings. The prefix bootstrap proves the certified Cayley norm is below one, hence `theta < 1 < pi`, and proves every accepted correction remains below `1e-2`, fixing the exact deployed quaternion branch throughout the inner funnel. In A mode the certified prefix also stays strictly inside the shipping accelerometer-bias projection ball, so the exact projection branch is the smooth identity-interior branch there; the nonsmooth projection surface is not silently linearized.

#### Metric-consistent structured defect transport

The nonlinear layer supplies one number: an upper bound `B` on `||r_word||_M / W_0`. The original route obtained it by leaving the metric -- `||z||_2 <= sqrt(W/m_-)` and `||r||_M <= sqrt(m_+) ||r||_2` -- and by bounding the gain isotropically as `||K|| <= sqrt(Sigma_max/R_min)`. That pays the full `sqrt(cond(Sigma))` of a covariance whose diagonal upper bounds span thirty-four decades, and it charges attitude-driven defects against the translation block that carries the whole spread.

`tools/ou3_p4_metric_defect_transport.py` supplies the structured replacement. Three exact facts, each stated against the same source-derived P3/process bounds the certificate already consumes:

1. **Gain transport.** For the shipping gain `K = P H^T S^-1` write `A = R^-1/2 H P^1/2`. Then `P^-1/2 K = A^T (A A^T + I)^-1 R^-1/2` and `||A^T (A A^T + I)^-1|| = max_i sigma_i/(1+sigma_i^2) <= 1/2`, so every residual defect satisfies `||K q||_M <= sqrt(s) ||q|| / (2 sqrt(lambda_min R))`. The node metric uses `Sigma_KF(g) >= P`, so the same bound holds there.
2. **Chart transport.** Marginalising the exact quadratic form gives `min_xi [c;xi]^T Sigma^-1 [c;xi] = c^T (Sigma_cc)^-1 c`, hence `||c||^2 <= lambda_max(Sigma_cc) W / s`. The P3 covariance upper is a Loewner diagonal dominator, so the marginal block maximum bounds it. The exact word defects are quadratic in the attitude, gyro-bias and `a_w` coordinates only; the translation block never enters them.
3. **Attitude-injection cost.** The quaternion injection remainder is supported on the attitude coordinates, so it is charged on the conditional attitude covariance. Every source node is post-prediction, where `Sigma >= Q` gives `(Sigma^-1)_tt <= I/(rho_att q_theta)`, or is reached from one by at most three in-sample corrections, each adding at most `||H_theta||^2/lambda_min(R)` in the exact information form `P+^-1 = P^-1 + H^T R^-1 H`.

Both gains bound the same quantity, so the producer keeps `min(B_isotropic, B_metric)` and the refinement can never widen the certificate. The structured route binds, taking `B` from `8.26237725542113e+34` to `2.19808143397921e+09`.

The final validated CI certificate reports:

- H: `W_* = 4.65099776131798868e-90`, `mu_W >= 1.89616809385829038e-35`, endpoint relative decrease `>= 1.89616809385829092e-35`, prefix canonical norm `<= 4.31323440648337294e-45`, certified attitude Cayley radius `<= 5.58190376455273201e-47`;
- A: `W_* = 4.65099773925429715e-90`, `mu_W >= 1.89616808936071053e-35`, endpoint relative decrease `>= 1.89616808936071107e-35`, prefix canonical norm `<= 4.31323439625267857e-45`, certified attitude Cayley radius `<= 5.58190375131284112e-47`.

That is fifty-one decades of `W_*` above the retired isotropic envelope (`3.29172575174270652e-141`), and it is still not a practical basin. The reason is no longer the defect constant: with `B` fixed, `sqrt(W_*) = delta/(8B)` is capped by the P3 word endpoint margin `delta = 3.79233618771658e-35`. Strict positivity is stored directly, so the proof never forms `1 - delta/2` when binary64 would round that quantity back to one.

#### The route, not its constants, is now the binding obstruction

`tools/ou3_p4_p5_route_ceiling_certificate.py` closes the question of how far this accounting can go. The certified attitude radius of a level is exactly `theta(W) = a_t sqrt(W)` with `a_t = sqrt(Sigma_tt_upper/s)`. Under four hypotheses the shipping producer satisfies -- `delta <= 1` because `Omega_word <= Sigma`; the word must cover the all-accepted branch, so the injecting-operation count is at least the sample count; the exact Cayley cross term `0.5 d x c` is not slack; and the uniform accepted-injection bound is `a_t sqrt(W)`, which is the exact supremum of the attitude part of `K H z` on the metric ball -- an attitude-supported defect costs at least `1/a_t` in the metric, so

`kappa >= 0.5 a_t`, `B >= F N 0.5 a_t`, `theta_capture <= delta / (F N)`.

`a_t` cancels. The ceiling does not depend on the covariance bound, on the metric normalization, or on any constant the search has been sharpening. At `delta = 1` and with no prefix overshoot at all it is `4.95049504950495e-03` rad; at the shipping prefix factor it is `1.23762376237624e-03` rad. The bounded P1 handoff nodes are `0.2721648148683776` and `0.5947333355555983` rad, so the route is short by at least a factor of `120` even at its own theoretical maximum, and the word would have to inject fewer than `1.69` accepted attitude corrections per second for it ever to reach the handoff.

This is a ceiling on the proof route, not on the filter. It retires further sharpening of `kappa`, of the covariance enclosure, or of `delta` as a P5 direction: the ceiling is independent of all three.

**PASS:** `P4_EXACT_NONLINEAR_WORD_CERTIFICATE=PASS` for both H and A with explicit positive `W_*`, positive `mu_W`, exact source-operation semantics, complete branch coverage, and prefix/chart safety. The level remains a theorem seed, and the route ceiling proves it cannot be promoted to the P1 handoff radius without the structural change named in P5.

### P5 — Initialization-to-inner-funnel finite capture — UNIFORM-TRANSPORT ROUTE PROVED INSUFFICIENT, MATCHED ROUTE REQUIRED

The original composition audit correctly found that the useful P1 handoff family lies far outside the microscopic P4 inner seed. That result remains a guard against illegally extrapolating the P4 local recurrence, but it is no longer the current P5 obstruction. The outer bridge has since progressed through the source-staged startup covariance, first-`S`, large-angle geometry, quotient correction, and exact transport layers.

The currently closed P5 prerequisites are:

- the goLive H covariance seed and scheduler phase are source bound;
- the first-due `S=0` gain is bounded with source-staged theta/`S` covariance structure rather than the translation-dominated global covariance box;
- the conditional first-`S` state prefix is finite with the complete `S -> attitude` gain retained;
- the first-due `S` quaternion/Cayley prefix is certified as `PASS_WIDENED_CHART`; the convenient `||c||<1` test is diagnostic only and is not a theorem gate;
- exact finite-angle vector information is positive on both gauged H handoff nodes;
- the isotropic raw `V_R` outer-sector route has a validated source counterexample and is retired rather than repaired by a larger adverse remainder;
- the false yaw-only/full-gyro-bias gravity route is retired; the detectable gravity quotient carries the gravity-parallel gyro bias as an explicit bounded neutral input while retaining the complete translation word;
- every shipping normal-Live operation class is bound to the exact Joseph/quaternion/reset transport calculus, including sequential immediate resets and PSD `a_w` covariance synchronization;
- the exact Cayley vector-defect geometry and the signed Cayley correction primitive are validated, so the backend does not replace `1-a^T c/4` by `1-|a||c|/4`.

A further exact reduction now removes the remaining standalone vector-`eta` penalty from the active P5 numerical route. For the configured isotropic magnetometer,

`H_m = -[v]_x`, `R_m = r_m I`, `H_m^T v = 0`, `S_m v = r_m v`, and therefore `K_m v = 0`.

If `y_m` is the exact finite-angle magnetic residual, define

`d_m = H_m^T y_m / ||v||^2`.

Then the implemented correction satisfies the exact identity

`K_m y_m = K_m H_m d_m`.

The radial finite-angle residual therefore produces no state correction, quaternion injection, or reset correction. In Cayley coordinates the effective tangent coordinate is nonexpansive and is enclosed cellwise by `tools/ou3_p5_effective_vector_input.py`.

For the accelerometer, the shipping Jacobian contains the orthogonal full-rank block `J_aw = R_wb`. For `y_a = H_a z + eta_a`, define `e_eta = J_aw^T eta_a` and insert it in the `a_w` coordinate. Then

`H_a E_aw e_eta = eta_a`,

so exactly

`K_a(H_a z + eta_a) = K_a H_a(z + E_aw e_eta)`.

This does **not** declare `eta_a=0`; it represents the same source residual as a source-correlated effective `a_w` tangent input. The Joseph identity remains valid. The proof backend therefore propagates the joint `P,H,R,K,r,d_eff`/effective-`a_w` cells instead of subtracting an unrelated `eta^T R^-1 eta` norm budget. This reduction is scoped to the configured theorem domain with optional IMU lever-arm compensation disabled.

The bridge remains fail-closed. Its current first unclosed numerical obligations are:

- gauged H: `COMPLETE_WORD_EFFECTIVE_VECTOR_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED`;
- gravity quotient H: `GRAVITY_QUOTIENT_EFFECTIVE_ACCEL_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED`.

Thus the next work is the outward source-correlated subdivision over the **later** prefixes of the complete 1 s word: prediction, accepted/rejected vector corrections, due/not-due `S`, covariance/source-schedule evolution, and each immediate reset. Each cell must keep `P,H,R,K,r,d_eff` (or the effective accelerometer `a_w` input) jointly reachable and accumulate the exact reset/prediction budget. Only after those prefix cells prove an outer recurrence into the existing P4 inner seed may the certificate set a finite integer `N_H_words` and capture time. `N_H_words` remains intentionally unset today.

#### Sample-1 first-accelerometer `q<8` chart line: perturbation route closed out

A separate machine-checked line (`tools/ou3_p5_sample1_*_v29..v50.py`) works the
sample-1 first-accelerometer Cayley chart, whose `q<8` target keeps the deployed
normalized-quaternion branch inside its validated composition range. V41's
complete source-cell-0 cover leaves one authoritative first survivor at
`(p,t,a)=(0,0,23)` with `q = 8.344528951460543`. V42 through V49 all attacked the
same object: the correction-perturbation caps that widen the nominal V10
directional correction.

`tools/ou3_p5_sample1_zero_perturbation_barrier_v50.py` closes that line out.

First it supplies the strongest remaining refinement on it. In the sample-1 body
gauge the shipping accelerometer Jacobian is exactly `H = [-[f]_x | I]`, because
the source block `J_aw = R_wb()` is orthogonal and the gauge places it on the
identity while `J_att = -skew(f_cog_b)` carries the whole modelled dependence.
The Jacobian perturbation is therefore supported on the attitude columns,
`Delta H = [E_theta | 0]`. That is exactly the fact which lets the certified
V12C/V12D `Delta C` parent be the four-term sum
`dP ||H_theta|| + ||P_theta|| dH + dP dH + dP` with no `P_theta,aw Delta H_aw`
term, so V50 reproduces that parent bit-for-bit before reusing the support
anywhere else. Three of V34's seven first-row `Delta S` terms then collapse onto
attitude-restricted nominal factors,

`e_i^T H P E^T -> ||(h_i P)_theta|| dH`, `e_i^T E P H^T -> dH ||(P H^T)_theta||`,
`e_i^T E P E^T -> dH ||P_theta,theta|| dH`,

each taken as a minimum against its own V34 parent term and still intersected
with V12D's full `||Delta S||`. The certified row bound improves from
`1.5212753710318652e-07` to `1.3328773849908888e-07`, which unpins the theta-y
gain row from the V12D parent but leaves theta-z pinned.

Second, V50 decomposes the certified sample-1 reduced covariance perturbation
`dP = 2.76914978691018e-10`, reproducing the V40/V12C parent exactly. Its
largest constituent is the raw sample-1 prediction attitude-covariance epsilon,
`1.4828249737599615e-10`, at 53.5 percent, ahead of the reset-gauge direction
term and the transported first posterior.

Third, and decisively, V50 reruns the authoritative V48 composition with every
componentwise correction-perturbation cap forced to zero. The composed q is
`8.344528951460543` — bit-identical to the parent. The perturbation caps do not
move q at all, so no refinement of `Delta S`, `Delta C`, `dP`, or the
componentwise split can close the authoritative survivor. The reason is
quantitative and is reported alongside: on the binding geodesic branch the
principal angle must fall by `0.019476337434169544` rad to reach `q<8`, while the
entire certified perturbation budget is `0.003654031568545583` rad, 18.8 percent
of it. Removing that whole budget from the parent radius as well still leaves the
geodesic branch at `8.277775079246444`.

The next work on this line is therefore the nominal geometry, not the
perturbation: sharpen V10's exact first-accelerometer directional correction
magnitude, or the sample-0 current chart `q1`, on the authoritative V45 parent.
This is a redirection of the proof search only. No filter setting, source
domain, six-radian correction limit, `q<8` target, source language, whole-word
criterion, or `N_H` state changes, and `N_H_words` remains unset.

`tools/ou3_p5_sample1_exact_monotone_source_gain_v51.py` discharges that
obligation at the witness. The first-accelerometer block is described by eight
rational functions of the same three source intervals - the attitude variance
`t`, the `a_w` variance `p`, and the accelerometer noise variance `r`:

`a = t(p+r)/D`, `b = p(g^2 t + r)/D`, `c0 = -g t p/D`, `bz = p r/(p+r)`,
`det = t p r/D`, `k_theta = g t/D`, `k_aw,t = p/D`, `k_z = p/(p+r)`,
with `D = g^2 t + p + r`.

The parent backend evaluates each as a straight interval expression, losing the
dependency between numerator and denominator. The loss is not cosmetic: at the
witness the parent encloses `k_z = p/(p+r)` in
`[0.5594923342554586, 1.0537323143362434]`, even though `p/(p+r) < 1` holds
identically for positive `p, r`.

Each expression is monotone in `t`, `p` and `r` separately, so its exact range
over the parameter box is attained at a corner. `a` and `b` reduce to two
variables (`u = p+r` and `w = g^2 t + r`) that each range over an exact interval
independent of the remaining variable, so their corners are genuine corners of
the original box. V51 evaluates every expression at its extremal corner with the
same outward-rounded backend and intersects the result with the parent, so a
refinement can never widen a bound or leave its parent, and fails closed if it
would.

At the authoritative witness, reproduced from source rather than copied:

| quantity | parent | exact monotone |
| --- | --- | --- |
| `k_z` upper | 1.0537323143362434 | 0.8001468619320714 |
| sample-1 `\|f\|` | 21.395742136954993 | 18.606777069495593 |
| post-first `a_w` axial | 4.911252706804307 | 2.9343396015749246 |
| sample-1 `rho` | 17.922551201967796 | 15.229738748335985 |
| `k_perp` | 0.9753682347137846 | 0.8468904975139163 |
| `k_parallel` | 0.09899544770387604 | 0.09772949400653325 |
| V10 correction | 2.0466720610769817 rad | 1.7313776836494923 rad |

V50 measured that the geodesic branch needed the correction principal angle to
fall by `0.019476337434169544` rad. This removes `0.3152943774274895` rad,
sixteen times that. Composing the refined correction with V41's archived
sample-0 chart - an upper bound computed with the unrefined gains, hence a
conservative partner - gives `q = 4.8010333986449245` against the archived
parent `q = 8.344528951460543`. The authoritative first survivor closes.

V51 evaluates the single authoritative witness cell.
`tools/ou3_p5_sample1_exact_monotone_cover_lift_v52.py` lifts the same
enclosure over the complete cover. `ou3_p5_sample1_structured_full_gain_v8`
gained a module-level `_first_block_quantities` helper holding the eight
expressions unchanged, so a refinement can be installed without editing the
producer, exactly as V12D/V40 do for the PSD perturbation; a full 12816-cell
V10 rebuild before and after that refactor is bit-identical in every summary
field. With V51's exact path installed there:

| cover quantity | parent | exact monotone |
| --- | --- | --- |
| cells evaluated | 12816 | 12816 |
| cells narrowed / widened | - | 12816 / 0 |
| per-cell narrowing ratio | - | 1.0106991820188718 .. 1.495825177084859 |
| max sample-1 residual | 46.01460061569009 | 42.67399431981034 |
| max `k_perp` | 1.3986770467171177 | 1.2187023949259914 |
| max `k_parallel` | 0.19289137244367335 | 0.19082905441803827 |
| max V10 correction | 7.016940736774492 rad | 6.7576910288356276 rad |

Containment is checked cell by cell, not only in aggregate: a single widened
cell fails the producer closed.

`tools/ou3_p5_sample1_exact_monotone_q8_cover_v53.py` then re-runs V41's
signed-chart `q<8` composition itself, regenerating both the shipping and the
refined cover in one process rather than quoting an archived count, and
requiring the parent run to reproduce V45's archived first-survivor
`q = 8.344528951460543`. The refinement closes the archived first survivor and
strictly improves the cover, but does not close it:

| `q<8` cover | parent | exact monotone |
| --- | --- | --- |
| signed-Cayley cells | 461376 | 461376 |
| open cells | 235738 | 219574 |
| first open cell | `(0,0,23)` at `q = 8.344528951460543` | `(0,2,23)` at `q = 8.475205389989586` |
| worst composed `q` | 525593.677323337 at `(20,16,20)` | 499303.8238549043 at `(23,18,4)` |

So the exact monotone enclosure is a real and uniform tightening of the whole
first-accelerometer block - it closes 16164 further signed-Cayley cells and
retires the survivor that V42 through V50 could not move - but roughly 48
percent of the cover is still open, and the worst cell remains five orders of
magnitude above the target. The remaining obstruction is again nominal
geometry, now at `(0,2,23)` and at the far worse `(23,18,4)`, whose correction
radial upper bound of `3.1593913566884573` rad is past `pi` and so wraps rather
than composing.

`tools/ou3_p5_sample1_open_cell_correction_budget_v54.py` measures how far the
remaining cells are. For each it reconstructs the refined chain from source,
decomposes `corr^2 = k_perp^2 rho_x^2 + k_par^2 (rho^2 - rho_x^2)`, and compares
against the largest correction the SO(3) triangle admits for that cell's
recorded chart, `phi_d <= 2 atan(4) - 2 atan(q_current/2)`:

| cell | correction | admissible | gap | `rho` reduction that would suffice |
| --- | --- | --- | --- | --- |
| `(0,2,23)` first open | 1.9595039241752017 | 1.9389861206167838 | +1.05% | 1.25% |
| `(23,18,4)` worst | 3.1985263367953887 | 1.257263970049669 | +60.7% | 65.0% |
| `(0,0,23)` retired | 1.7313776836494923 | 2.0308497552113565 | -17.3% | already inside |

So the two remaining named cells are in completely different regimes. The
first-open cell misses by about one percent - the same order as the gap V50
measured at the retired witness, which the exact monotone enclosure then beat by
sixteen times - while the worst cell needs the residual cut by roughly two
thirds, and there zeroing `rho_x` or `k_perp` outright cannot reach the target
because `k_par rho` alone already exceeds it.

The two cells also fail for different reasons. At `(0,2,23)` the residual is 58
percent predicted force and 39 percent post-update `a_w` error; at `(23,18,4)`
that inverts to 28 and 69 percent, with the post-update `a_w` bound driven by
the tangent term `|(1 - k_aw,t) r_t| = 9.2` m/s^2 that the large `p` cell
produces. Sharpening the post-update `a_w` bound needs signed structure on the
first tangent residual, in the way V10's equations (1) and (2) already supply
for the `x` row.

V54's reduction figures are arithmetic on outward-rounded enclosures, not
enclosures themselves, and they cover the geodesic branch only: V41 closes a
cell on `min(geodesic, product)`, so reaching the geodesic target closes a cell
but missing it does not prove one open. V54 therefore reports a distance and
never a verdict.

None of V51, V52, V53 or V54 composes `q<8`, promotes sample 1 or P5, or sets
`N_H_words`.

The independent implementation-stability composition gate continues to consume P5 directly. The older generic affine deployment-capture arithmetic is not accepted as a substitute.

#### The uniform transported-defect route cannot close P5

Every P5 attempt so far, including the whole V1..V54 sample-1 line, sharpens a
constant inside one accounting:

`||r_word||_M <= B W_0`, `B = F N kappa`, strict decrease requires `sqrt(W_0) <= delta/(2B)`.

`tools/ou3_p4_p5_route_ceiling_certificate.py` bounds what that accounting can
report, before any constant in it is chosen. Because the certified attitude
radius of a level is exactly `theta(W) = a_t sqrt(W)` and an attitude-supported
defect costs at least `1/a_t` in the metric, the chart scale `a_t` cancels and

`theta_capture <= delta / (F N)`.

With `delta <= 1` and `N` at least the word's sample count, the ceiling is
`4.95049504950495e-03` rad with no prefix overshoot and `1.23762376237624e-03`
rad at the shipping prefix factor. The bounded P1 handoff nodes are
`0.2721648148683776` rad (normal gauged) and `0.5947333355555983` rad (timeout
gauged). The route is therefore short by at least a factor of `120`, and would
need the source word to inject fewer than `1.69` accepted attitude corrections
per second before it could reach the handoff at all.

The ceiling is independent of `kappa`, of the covariance enclosure, of the metric
normalization and of `delta`. Retightening any of them -- including the
fifty-one-decade metric-consistent widening of `B` recorded in P4 above -- moves
the reported radius but not the ceiling. The `Delta S`, `Delta C`, `dP` and
componentwise-split search directions of V42..V54 are constants of exactly this
kind, so that line is retired here rather than continued.

#### What a closing P5 route must change

The failure is structural and localised: the accounting compares a defect that is
*fast* (attitude, refreshed every accepted vector correction) against a margin
that is *slow* (`delta`, rate-limited by the gyro-bias channel of the same word),
and then accumulates the defect `N` times with no credit for the contraction
happening in between. Three changes follow, and all three are needed:

1. **Match each defect to its own operation.** For an accepted correction the
   exact information identity gives `W - W+ = ||H z||^2_{s S^-1}` on the same
   step that injects the defect. Charging the injected defect against that step's
   own decrease removes both the `N`-fold accumulation and the dependence on the
   whole-word `delta`.
2. **Carry a directional or block margin.** The attitude channel must not be
   rate-limited by the slowest source channel of the word. This means a block
   Lyapunov function with a small-gain coupling constant, not one scalar `delta`
   over the full 18/21-state box.
3. **Replace the Cayley quadratic remainder on the attitude channel with an
   exact finite-angle sector contraction** of the deployed correction. The
   quadratic remainder is what forces the basin to shrink like `delta/(F N)`;
   the finite-angle statement is what makes a radius of order `0.3` rad
   expressible at all. `tools/ou3_p5_large_angle_sector_certificate.py` is the
   existing entry point for that geometry.

Until those land, `N_H_words` stays unset -- now for a proved reason rather than
a pending computation.

**PASS criterion for closing P5:** every certified normal or timeout handoff lies in a validated outer H capture funnel; every source-word prefix stays safe; the outer recurrence reaches `W_*` in a finite machine-certified number of H words/time. The current certificate does **not** yet satisfy this criterion, and the route ceiling above proves the present accounting never can.

### P6 — Prove every implemented hybrid jump

For each source-reachable hard event, validate the exact jump image and a strict destination inequality:

1. startup handoff;
2. H→A accelerometer-bias release (18→21 dimensions, including new-coordinate energy);
3. magnetic lock;
4. magnetic refinement/re-gauge;
5. tilt reset;
6. tilt re-lock;
7. cooldown re-entry;
8. periodic `a_w` covariance synchronization (nonexpansive PSD/Loewner proof, with strict decrease supplied by surrounding words).

After any jump that lands in an outer funnel, recompute finite recapture to the inner funnel.

**PASS:** every source-reachable jump lands strictly inside a certified destination funnel and the hybrid execution cannot escape through an unproved branch.

### P7 — Gravity-only quotient route

Machine-certify the no-magnetometer route on the yaw quotient using the same startup, translation, nonlinear, prefix and hybrid machinery, but without making an absolute-yaw convergence claim.

**PASS:** regional practical stability modulo yaw for the source-complete gravity-only implementation.

### P8 — Stochastic localized theorem

Using the same certified source words/funnels from P4-P6, validate the paper's Gaussian/localized drift and finite-horizon concentration constants from primitive source/noise bounds.

**PASS:** positive stochastic drift margin, localization radius inside the deterministic funnel, and recomputed finite-horizon failure probability below the declared theorem budget.

### P9 — Independent final implementation-stability composition gate

The final gate must regenerate source contracts and recompute all composition arithmetic. It must not trust an upstream `PASS` field.

The final status may be `PASS_IMPLEMENTATION_STABLE` only if all of the following are true simultaneously:

- source/implementation parity P0;
- startup/reset/handoff P1;
- source-word completeness P2;
- H linear certificate P3;
- A linear certificate P3;
- H nonlinear/prefix certificate P4;
- A nonlinear/prefix certificate P4;
- finite initialization capture P5;
- all hybrid jumps P6;
- gravity-only quotient result P7 when that mode is claimed;
- stochastic theorem P8 when the stochastic claim is requested;
- existing `main` performance/regression/evidence gates remain green.

Anything less is reported by its first failed proof obligation, never as “stable by existence.”

## PR completion criterion

This PR is complete only when it either:

1. produces `PASS_IMPLEMENTATION_STABLE` with the numerical margins/capture bounds above; or
2. identifies a mathematically genuine obstruction in the implemented filter/theorem assumptions, with the first failing validated inequality and no weakening of the proof gate.

A tiny unspecified local basin, replay-only success, or an analytical existence theorem is explicitly not an acceptable completion result.