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

Therefore the rejected mechanism is now broader than the original `C-eps I` step: **recursive natural interval covariance boxes that forget the repeated scalar source parameter at every Riccati operation are frozen as a dead end.**

Canonical P3 still has not reached translation/H/A margin calculation, so none of these diagnostic values are P3 theorem margins.

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
3. **Analytic complete-segment Gramian/Riccati lower in a different information/covariance representation.** Avoid intervalizing the covariance recursion itself.
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

A rejected route may be resurrected only after recording the new mathematical fact that invalidates its rejection.

## Retained facts

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
