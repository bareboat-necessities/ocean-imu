# OU-III proof research state

This file is the current research ledger required by the root `AGENTS.md`. Keep it short and replace stale state rather than accumulating PR history.

## State

**REPLAN COMPLETE — NEXT EXECUTION AUTHORIZED**

Current PR research question: make canonical P3 reach one rigorous numerical PASS/FAIL verdict at the unchanged `1e-18` gate; only if P3 passes, run one non-promoting structure-exact P4 complete-word feasibility diagnostic. P5 is out of scope.

## Failure classification

The current canonical run is **not a numerical P3 theorem FAIL**. It is blocked earlier by a translation-covariance enclosure/conditioning failure:

`Loewner prediction lower lost strict SPD; split x cell`

There is also a separate implementation defect in the initial `x` partition helper: subcells adjacent to the `x=0.01` process-series branch are emitted as nested/overlapping intervals because branch-endpoint clamping is applied to every subcell instead of only the boundary subcell. That creates redundant work but does not explain the adaptive recursion failure by itself.

## Evidence

For the deterministic source-node 137 / gap-13 failure, diagnostic evaluation of the first prediction on the widest initial small-`x` cell gives:

- `x = [0.0073253899652235745, 0.009999999999999998]`;
- center `lambda_min ~= 3.405e-10`;
- absolute row-radius `eps ~= 1.142e-5`;
- `eps/lambda_min ~= 3.35e4`;
- limiting center eigenvector approximately `[0.090, -0.448, 0.890, -0.0077]`, dominated by the poorly scaled `p/S` directions.

Following the existing geometric subdivision down one branch reduces `eps/lambda_min` approximately by a factor of two per split, but at depth 12 it is still about `8.26`. Thus approximately four additional binary levels would be needed merely to make this first 5-ms prediction plausible. A full binary refinement to that scale is on the order of `2^16` leaves for the widest initial cell before considering the complete 800-source x 14-gap family. Deeper subdivision is therefore rejected as a tractable proof architecture even though its local error scaling is understood.

A non-rigorous point diagnostic over the theorem-relevant **complete 13-sample segment** is qualitatively better: across the source-137 tau endpoints the final covariance has `lambda_min ~= 2.02e-5 ... 2.83e-5`, roughly five orders of magnitude above the one-sample floor. This is diagnostic evidence only, but it quantitatively supports moving the enclosure target to the complete segment.

A one-step full congruence normalization of the same wide interval reduces the relative enclosure mismatch substantially (`eps/lambda_min` from roughly `3.35e4` to roughly `1.17e3`) but still leaves a large one-step deficit. Therefore normalizing each 5-ms step is not the selected architecture.

Canonical P3 still has not reached translation/H/A margin calculation, so none of these diagnostic numbers are P3 theorem margins.

## What the failure invalidates

- Requiring a useful absolute point Loewner lower after every uncertain 5-ms prediction.
- Blindly increasing `MAX_ADAPTIVE_X_DEPTH` as the response to the same strict-SPD failure.
- Treating the current one-step enclosure mechanism as part of the theorem rather than as one proof tactic.

## What the failure does not invalidate

- The `1e-18` canonical P3 usefulness threshold.
- The same-history P2-V1 source language.
- Existence of a rigorous complete 13--26-sample translation floor.
- H=18/A=21 full-state joining, which has not yet received a translation margin from this backend.
- P4 feasibility or infeasibility, because P4 remains blocked by canonical P3.

## Critic pass

Assume the current absolute one-step architecture is wrong. The strongest reason to abandon it is that it spends proof precision on a property the theorem does not require: a well-conditioned point lower after every 5-ms intermediate prediction in coordinates whose integrated `S` direction is extremely small. The complete segment accumulates much more process covariance, so certifying the segment directly attacks the actual `common_boundary_floor` quantity.

Qualitatively different alternatives considered:

1. deeper adaptive subdivision of the same absolute lower;
2. one-step relative/congruence-normalized lower;
3. complete-segment relative/congruence or verified generalized-eigenvalue lower;
4. complete-segment Gramian/Riccati formulation that never asks for intermediate strict point floors;
5. a different Lyapunov/covariance representation if complete-segment covariance comparison remains ill-conditioned.

**Selection:** pursue (3), with (4) as the immediate fallback. It attacks the theorem-relevant segment floor and the point diagnostic shows about five orders of magnitude more smallest-eigenvalue headroom at 13 samples. Options (1) and (2) still spend excessive complexity on the unnecessary one-step floor.

## DEAD_ENDS

- **REJECTED: endpoint-only 800-node P2 ancestry as a P3 source model.** It loses staged/committed path memory required by the frozen P2-V1 interface.
- **REJECTED: independent Cartesian `tau/sigma/R_S` extrema.** They destroy source-history correlation and are forbidden by the P2 consumer contract.
- **REJECTED: recursive absolute entrywise Loewner point lower plus blind deeper subdivision.** Its local error scales with width, but the depth/leaf count needed for the first prediction is not a viable source-complete architecture.
- **REJECTED: one-step congruence normalization as the primary architecture.** It improves conditioning materially but still leaves a large one-step relative deficit while solving a property P3 does not require.
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

## Next falsifiable experiment

1. Fix the independent `x`-partition boundary bug and add a focused non-overlap/branch-coverage regression.
2. Build a **diagnostic-only complete-segment reference matrix** for source 137 / gap 13 from a high-precision midpoint calculation.
3. Using verified/outward operations, compare the full 13-sample covariance family directly against that reference by congruence/generalized eigenvalue and attempt to certify `P_segment(x) >= alpha P_ref`, `alpha > 0`, without requiring intermediate point-SPD floors.
4. Report `alpha`, relative enclosure radius, limiting direction, and subdivision scaling.

**Decision rule:** if the complete-segment relative comparison closes with modest subdivision, generalize it to the source-complete boundary-floor calculation. If it exhibits the same conditioning pathology as the one-step architecture, freeze it after one refinement and move to the segment Gramian/Riccati alternative. Do not increase `MAX_ADAPTIVE_X_DEPTH`.

Before any new P4 proof producer, first obtain the canonical P3 numeric verdict. If P3 passes, the only P4 work allowed in #471 is the non-promoting high-precision complete-word ratio diagnostic `rho_w = V_after(F_w(x)) / V_before(x)` required by the PR scope.
