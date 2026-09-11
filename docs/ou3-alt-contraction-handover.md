# OU-III parallel ALT proof handover

## Where to resume

This is the canonical handover for the independent whole-word joint24 hybrid dissipativity proof introduced by PR #517. After PR #517 is merged, start the next ALT continuation as a NEW PR from the latest `main`. Read, in order:

1. `AGENTS.md`
2. `docs/ou3-alt-contraction-handover.md` (this file)
3. `docs/ou3-alt-contraction.md`
4. the ALT section of `docs/ou3-proof-research-state.md`
5. `docs/ou3-brmm-main-handover.md` for the shared physical/runtime contracts

The existing P2/P3/P4/P5 proof is a separate, independently continuable route. Do not delete, rewrite, bypass, or silently make it depend on ALT. Likewise, ALT must not consume old PASS labels as proof of its own master inequality. Shared source/runtime facts may be reused only with their actual hypotheses and same-history ancestry.

## Immutable shared constraints

Preserve the unchanged shipping implementation and corrected COMPLETE-BRMM physics from main: padded `||a_wave|| <= 8.8 m/s^2`, `D_S <= 1100 m*s`, one-time Live S origin, zero/disabled lever arm in current proof scope, actual anisotropic `R_S`, same-signal WPE -> sigma -> tau -> T_S -> R_S relation, staged tuner/scheduler state, actual Joseph/reset/projection path, BIAS0/BIAS1/BIAS2 ancestry, H18 and A21, and deployment finite precision as a theorem obligation. P3 remains frozen at `1e-18`. No replay fitting, finite-seed qualification, wordwise S re-zeroing, covariance consistency as hard entry membership, or filter changes for proof convenience.

## ALT theorem target

Use the full true-minus-estimate augmented state

`z = (e_theta,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta_true) in R^24`.

Do NOT demand strict homogeneous contraction of every coordinate. H18 held bias and BIAS2 constant physical bias provide neutral/persistent directions. The selected target is coercive joint storage with a justified supply:

`m||z||^2 <= V(z,xi) <= Mbar||z||^2`

`V(z_next,xi_next) <= rho V(z,xi) + w^T Gamma w + c^T Beta c`, with `0 < rho < 1`.

Only independently bounded physical/bias/reference quantities may enter `c`; unknown motion errors may not be relabeled as bounded inputs. Keep motion/bias cross terms in storage. The ultimate bound must be quantitatively useful, not merely finite.

## Implemented proof primitives

`tools/stability/ou3_alt_contraction/core.py` contains non-promoting algebra for:

- inverse-free measurement graph `S q = r`, correction `N q`, where `N` is the ACTUAL shipping numerator after masks;
- exact finite endogenous innovation increments
  `S1*dq + dS*q0 = dr`, `dcorrection = N1*dq + dN*q0`;
- joint24 bias prediction retaining `(phi_true-phi_hat)*beta` inside the state map and the SAME physical driver in bias error and true bias;
- a radial projection IQC retaining the same true bias on both sides;
- sign-correct descriptor/IQC dissipativity assembly;
- exact structural regressions showing why strict full-state homogeneous contraction and unqualified H18 information-form replacement are invalid shortcuts;
- an exact rational two-state masked-update analogue showing that joint storage plus bounded neutral supply can close despite a neutral coordinate. This analogue is not a shipping certificate.

`tests/ou3_alt_contraction/test_core.py` has 15 algebra/anti-shortcut tests. The independent `ou3-alt-contraction` workflow passed on PR #517.

## Executed diagnostic evidence

The unchanged shipping observer produced 29 retained H18 and 359 retained A21 approximately three-second same-mode words. Worst stored-map ratios were approximately:

- H18: `0.999513639919454`, margin about `4.8636e-4`
- A21: `0.995985717839491`, margin about `4.0143e-3`

The selected-direction signed ledgers also contracted. Treat these only as falsification/exploratory evidence. The observer reconstructs H from covariance/PCt, freezes endogenous gain dependence, accumulates a binary32 homogeneous product, and H18's 18-state ratio fixes held-bias error at zero. Re-evaluating that stored matrix at 80/120 digits does NOT produce the missing actual nonlinear joint24 word or a source-uniform certificate. No metric was fitted to replay.

## Current blockers

1. **Primary ALT blocker — actual source-uniform joint24 word attachment.** The exact finite solve-increment algebra exists, but `dN`, `dS`, residual, covariance, frontend/tuner, guards, true-bias recurrence, source forcing and all event branches are not yet bound together as one same-history source-uniform word.
2. **Fresh-Live hard entry is not closed.** The shared aggregate bias-supply builder currently stops with `qualified fresh-entry source failed` and reports hard-coordinate membership / source-uniform other-coordinate / fresh-Live-cover failures. Do not bypass it. The individual BIAS0/BIAS1/BIAS2 definition validators pass, but that is not startup/hardware admission.
3. **Uniform storage is not solved.** No common, parameter-dependent, or piecewise joint24 storage has yet been certified for the actual family. Search common joint storage with bounded neutral supply first; only then move down the hierarchy.
4. **Nonlinear/hybrid closure is open.** Hard finite-horizon IQCs or equivalent exact graph constraints for finite-angle attitude/reset/chord/floor/clamps/acceptance branches, compatible H18/A21 and H18->A21 edges, and literal every-prefix chart/domain retention remain open.
5. **Startup is separate and open.** The actual Mahony/proxy/learning path must be proved to land in the ALT Live basin in finite time for both measured-period takeover and prior-frequency timeout under the padded physical family.
6. **Finite precision is open.** Real-arithmetic solve identities do not enclose LDLT, Eigen, reset, projection, covariance-floor, and other deployment roundoff defects.
7. The existing proof's PE/vector-domain consistency problem remains separate; ALT should not inherit it unless a shared lemma is actually required.

## Next falsifiable sequence

1. Build an **actual source-uniform finite-increment joint24 shipping word**. Preserve nonzero nominal residual, endogenous `N/S`, true bias, physical source, frontend/tuner and guard ancestry. Instrument or derive missing actual finite increments rather than inferring them from the old frozen observer.
2. Before rigorous interval/outward work, run a high-precision complete-word feasibility diagnostic of THAT graph. It must report worst H18/A21 directions, bias/source supply ports, operation-by-operation margin consumption, hybrid edges, and distance to `rho=1`. If rho is above 1, abandon that formulation rather than refining unrelated bounds.
3. Search the common joint-storage + bounded-neutral-supply master first. If it fails twice by the same mechanism, obey the AGENTS two-strike rule: record failure analysis, compare at least three qualitatively different alternatives, and perform an architecture review before more refinement.
4. If a full master is feasible, make it rigorous over the analytic source family: source cover, metric coercivity/compatibility, hard nonlinear graph constraints, every-prefix reachability, first-exit/domain retention, and finite precision.
5. Determine the largest physics-compatible ALT Live basin and its quantitative ultimate bound.
6. Separately prove actual Mahony/proxy finite-time capture and H18->A21 landing into that basin.
7. Only after all obligations close may ALT claim an end-to-end theorem. Do not promote existing P4/P5 merely because ALT closes; reconcile theorem/gate policy explicitly in a later PR.

## Fail-closed status at PR #517 handover

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Existing `P4_PASS=false` / `P5_MAY_START=false` remain untouched by this route.

The purpose of merging PR #517 is to preserve a tested second proof architecture and its falsifiable continuation path on main, not to claim the stability theorem is finished.
