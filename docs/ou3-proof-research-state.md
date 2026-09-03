# OU-III proof research state

This file is the current research ledger required by the root `AGENTS.md`. Keep it short and replace stale state rather than accumulating PR history.

## State

**FAILURE ANALYSIS / REPLAN**

Current PR research question: make canonical P3 reach one rigorous numerical PASS/FAIL verdict at the unchanged `1e-18` gate; only if P3 passes, run one non-promoting structure-exact P4 complete-word feasibility diagnostic. P5 is out of scope.

## Current hypothesis

The present P3 failure is an enclosure/conditioning failure in the translation covariance lower, not yet a numerical P3 theorem failure. The absolute entrywise enclosure of `F L F^T + Q` followed by `C - eps I` can erase the small integrated-`S` one-step SPD direction before the complete 13--26-sample segment floor is formed.

The leading alternative is a relative/congruence-normalized Loewner enclosure applied to the prediction or complete segment, because P3 needs a rigorous complete-segment floor, not an artificially strong absolute point lower at every 5-ms intermediate step.

## Evidence

- The deterministic source-node 137 / gap-13 case reproduces `Loewner prediction lower lost strict SPD; split x cell`.
- Repeated adaptive subdivision reaches the configured depth limit without closing the current absolute point-lower construction.
- Canonical P3 therefore does not yet reach the translation, H, or A margin calculation; this run is not evidence that the `1e-18` theorem gate fails.
- The frozen P2 contract retains same-history source correlation and rejects independent Cartesian `tau/sigma/R_S` extrema.

## Current limiter

The limiting object is the rigorous Loewner lower used to construct `common_boundary_floor` for a complete legal 13--26-sample translation segment. The immediate experiment must measure conditioning of the first failed prediction and how the relative enclosure error scales under subdivision.

## DEAD_ENDS

- **REJECTED: endpoint-only 800-node P2 ancestry as a P3 source model.** It loses the staged/committed path memory required by the frozen P2-V1 interface.
- **REJECTED: independent Cartesian `tau/sigma/R_S` extrema.** They destroy source-history correlation and are forbidden by the P2 consumer contract.
- **REJECTED: recursive absolute entrywise Loewner point lower plus blind deeper subdivision.** The same strict-SPD mechanism has already failed repeatedly through the configured adaptive depth. Do not increase depth without a quantitative scaling argument showing the theorem threshold can be crossed.
- **REJECTED for current #471 scope: accumulating additional P4 micro-certificates before complete-word feasibility is known.** Local reset/sector/remainder bounds are subordinate to the complete-word ratio.

A rejected route may be resurrected only after recording the new mathematical fact that invalidates its rejection.

## Retained facts

- Canonical P3 useful gate remains exactly `1e-18`.
- `OU3_P2_CORRELATED_STAGE_TRANSFER_V1` and same-history source correlation are retained.
- Both H=18 and A=21 remain required.
- Zero/disabled lever arm and dormant/transparent vibration-guard branch remain the current proof scope.
- No replay fitting, operating-domain shrink for PASS, or deployed-filter change is permitted.
- Existing exact Joseph, co-rotated accelerometer, reset, and finite-angle identities remain useful structural facts, but they do not authorize P4 promotion.
- P4 is blocked until canonical P3 passes; P5 is blocked until canonical P4 strictly contracts.

## Alternatives under review

1. **Relative/congruence-normalized Loewner enclosure.** Normalize the covariance family near identity with a verified factor/solve, certify `B(x) >= alpha I`, then map the result back to the required segment floor.
2. **Complete-segment Gramian/Riccati lower.** Avoid demanding strict intermediate point floors and certify the complete 13--26-sample segment directly.
3. **Verified generalized-eigenvalue formulation.** Compare the interval covariance family directly with a chosen positive reference matrix instead of subtracting an absolute `eps I` in poorly scaled physical coordinates.
4. **Different Lyapunov/covariance representation.** Use only if the first three show that the present covariance representation itself is the limiting architecture.

## Next falsifiable experiment

For source node 137 / gap 13, record at the first failed prediction and at successive existing subdivision depths:

- `x` interval and depth;
- center matrix and radius matrix;
- absolute row radius `eps`;
- diagnostic `lambda_min(center)` and `eps/lambda_min(center)`;
- high-precision midpoint/end-point references for diagnosis only;
- limiting direction/coordinate.

**Decision rule:** if interval width shrinks but the relative SPD margin does not improve enough to plausibly close the required bound, the absolute-collapse tactic remains frozen and implementation moves to one of the qualitatively different alternatives above. Do not increase `MAX_ADAPTIVE_X_DEPTH`.

Before any new P4 proof producer, first obtain the canonical P3 numeric verdict. If P3 passes, the only P4 work allowed in #471 is the non-promoting high-precision complete-word ratio diagnostic `rho_w = V_after(F_w(x)) / V_before(x)` required by the PR scope.
