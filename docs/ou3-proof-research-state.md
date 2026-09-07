# OU-III proof research state

## Handoff checkpoint

Canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional complete-SEA3 P3 is closed and frozen at `delta=1e-18`. P4 is **OPEN** and P5 is **BLOCKED**. H18 and A21 are both required; H18->A21 remains a separate rectangular hybrid event. Zero lever arm and the dormant-transparent vibration branch remain the certified production branch.

The only production/proof-domain change made on PR #496 is the user-authorized accelerometer-bias projection-radius tightening from `0.5` to `0.4 m/s^2`. The declared startup/handoff accelerometer-bias error envelope is also `0.4 m/s^2`, while the Normal-Live active-bias interior bound is `0.35 m/s^2`, preserving a `0.05 m/s^2` projection margin. No other filter tuning, quality gate, source-language parameter, or P3 mathematics was changed.

The paper target remains the finite full-state source-indexed quadratic storage

`V(e,zeta)=e^T M(zeta)e`

on one complete same-history SEA3 word, with a strict finite endpoint inequality, a finite every-prefix gain, and every-prefix chart/source-domain retention. A replay, point experiment, fitted endpoint metric, reduced-state coercivity test, or alternate source language cannot promote P4.

## Exact structural work retained on this checkpoint

The branch retains the theorem-facing algebra and execution machinery needed by a continuation PR:

- exact finite accelerometer coordinate shift with original shipping H/P/K/S and `H0 != H_u` at finite attitude error;
- exact prediction transport with the literal full `F E_aw`, retaining v/p/S/a_w rows;
- exact Joseph/reset signed-information identity
  `Delta V = -I_y + E_eta + X_reset + E_reset`;
- exact deployed Cayley reset transport and A21 `0.4 m/s^2` bias projection with generalized Jacobian;
- every due S=0 event with actual applied anisotropic SpectralMSE `R_S` retained; S events have `eta=0` exactly and therefore contribute favorable information;
- every valid accelerometer update, applicable vector update, full Q, covariance floors and immediate resets retained;
- separate H18->A21 hybrid lift;
- outward interval/differential AD and generalized mean-value machinery for the exact nonlinear physical map;
- complete-word accelerometer covariance channel and exact homogeneous Cayley residual-sector factorization;
- the existing finite-`tau_b` A21 detectability module tied to complete SEA3, which closes the paper-level finite-bias detectability/UES hypothesis but explicitly does **not** close the canonical full 21x21 implementation-word P4 inequality.

Do not replace these with selected-S words, independent tuner/R_S boxes, independent per-sample source boxes, packet-count nonlinear budgets, scalar correction radii, inverse-metric-floor arguments, state elimination, or replay-derived source families.

## Authoritative single-observer linear evidence

The corrected single shipping observer owns one source/tuner/Riccati history and retains every operation. On the genuine PM+Stokes Hs=1.5 m history, the pre-0.4 legal 3-second point words were:

- H18: `rho_linear=0.9998658024147671`, 600 predictions, 600 accelerometer updates, 137 actual-R_S S updates, 75 vector updates;
- A21: `rho_linear=0.9958536807113242`, 600 predictions, 600 accelerometer updates, 108 actual-R_S S updates, 75 vector updates.

Same-observer event ledgers reproduce those ratios within a few `1e-6` and telescope to roundoff. The old duplicate observer that staged tuner-commit data at the wrong boundary is retired; its older A21 rho values are stale.

On the 0.4 head, the 3-second physical reset-normalized diagnostic keeps strict zero-state parity. H18 remains contracting over every retained tested scale. A21 still first crosses one at scale `8.0`; the bias projection is inactive in the problematic cases, and the worst retained-domain endpoint ratio is about `1.0860152320`. Thus the 0.4 clamp is a legitimate production/domain tightening, not a proof fix.

## 6-second and 9-second falsification result: longer-window route stopped

The source-contiguous physical long-window diagnostic completed successfully in GitHub Actions run `34083473117` (`ou3-p4-physical-long-window-feasibility`), artifact `10004621877`, artifact SHA256 `8b3c1e3b83e6774895cf20770764609afd38e044131204c568d49d7f3544aa19`.

It uses the same single shipping observer, the canonical source label, actual applied R_S, all accelerometer/S/vector events, strict zero-state parity, and no source/domain/filter substitution. It is explicitly non-promoting.

6-second result:

- H18 linear rho `0.9917547078907433`; worst retained physical rho `0.993720034085038`;
- A21 linear rho `0.991109058400818`;
- A21 first finite-scale crossing `7.5`;
- A21 worst retained physical rho `1.1357228916109403`;
- A21 worst prefix ratio `1.1358285446129752`.

9-second result:

- H18 linear rho `0.9786791982272148`; worst retained physical rho `0.978951454184692`;
- A21 linear rho `0.9868275426222095`;
- A21 first finite-scale crossing `16.0`;
- A21 worst retained physical rho `1.143708327664902`;
- A21 worst prefix ratio `1.1653466626592186`.

Therefore the proposed 3->6->9 second extension does **not** repair the finite A21 mechanism. Per the research protocol, stop the longer-window route here. Do not spend another iteration optimizing window length or point storage around this replay.

## A21 mechanism and dead ends

The finite A21 obstruction is primarily finite accelerometer curvature in a coupled attitude / latent-acceleration / accelerometer-bias cancellation direction. The exact lever-arm-off residual is

`y = (E-I) f_hat + E R_hat delta_a_w + delta_b_a`,

with first-order row

`H e = [c]_x f_hat + R_hat delta_a_w + delta_b_a`.

These first-order terms can nearly cancel while second-order attitude curvature remains. The A21 bias projection is not active at the problematic scales, and reset terms are small; tightening reset bounds or eliminating a_w attacks the wrong mechanism.

Retired / forbidden rescue routes include:

- estimator-pair shadow promoted as theorem map;
- raw or full-Phi endpoint optimization after their A21 finite-scale failures;
- arbitrary single-map converse-metric fitting;
- extending the failed physical replay to still longer windows merely to seek a green point;
- reduced-state/Schur certificates or a_w elimination;
- selected-S words or independent R_S/tuner schedules;
- scalar Lipschitz, correction-radius, inverse-metric-floor, or packet-count-times-worst-remainder bounds;
- further proof-driven filter/domain tightening beyond the authorized 0.4 change.

## Correct next theorem route

The next PR should start from current main after this handoff and pursue the paper-permitted **finite-`tau_b`, full-21-state A21 cascade/detectability route**. This is qualitatively different from the failed endpoint/window experiments.

The existing `tools/stability/ou3_sea3_a21_detectability_completion.py` already establishes the paper-level finite-bias detectability/UES hypothesis from:

- complete-SEA3 H18 contraction;
- exact finite residual-bias Gauss-Markov decay;
- bounded full-state H18<->b_a coupling on the compact word;
- full A21 process UCC;
- no alternate estimator and no state elimination.

However, that module deliberately remains fail-closed for the stronger canonical implementation-word bridge: `full_21x21_Omega_minus_delta_P_LDLT_closed_here = False` and `P4_MAY_CONSUME_P3 = False`.

The next theorem-facing master inequality must therefore be a **full-rank 21-state** Lyapunov/cascade inequality that consumes the same complete-SEA3 source history and retains every actual shipping event/coupling. It must show how the finite-`tau_b` bias decay/detectability estimate combines with the H18 complete-word dissipation and the nonlinear accelerometer residual sector without eliminating `b_a`, `a_w`, or any cross terms. Only after that full-state bridge is quantitative and outward-certified should the universal nonlinear complete-word P4 endpoint/prefix enclosure resume.

If the full-state finite-`tau_b` construction itself exposes a genuine admissible expanding direction that violates the paper-equivalent P4 inequality on the declared domain, report that obstruction directly rather than inventing a weaker certificate or shrinking the domain.

## Complete-source obligation still open

The upstream complete-source contract still does not materialize the correlated finite-window SEA3 realization as an executable outward source family. Parameter compactness, RAO/moment envelopes, hard pathwise acceleration/body-rate caps, frontend parity, and adaptive-state rate/jump bounds do not replace a same-history transition for the correlated source state. P4 may not substitute independent per-sample boxes, a replay record, or a finite harmonic/grid surrogate.

The continuation should first formulate the full-state finite-`tau_b` master inequality and determine exactly which source-correlated quantities it needs. Then extend the existing complete-SEA3 execution/source machinery only for those theorem quantities, rather than creating another surrogate proof language.

## CI / evidence status at handoff

The long-window falsification workflow is green. Exact Cayley parity, the physical 3-second finite-map feasibility, A21-word, complete-word feasibility, and lever-arm-study workflows were also green on head `80140b762ab7bfd924d1d5dc279153e425a31612` when this handoff was prepared; several broader workflows were still running.

`ou-validation` is red for a known evidence-provenance reason, not because its numerical unit-test body found a new filter failure: `tools/ou_evidence_contract.py --auto` reports replay dependencies changed relative to committed validation/robustness provenance, including the OU-III filter and WavePeriodEstimator dependencies. Genuine validation/robustness evidence regeneration is therefore still required before a later proof PR is declared final/ready. Do not hand-edit provenance hashes.

No claim is made that P4 or P5 is complete at this merge checkpoint.

## Continuation instruction

Create a new PR from the merged main state. Read `AGENTS.md`, this file, the merged #496 handoff comment, `doc/kalman_ou_iii/w3d-sea3-stability-theorem.tex-part`, and `tools/stability/ou3_sea3_a21_detectability_completion.py` first. Keep P3 frozen at `delta=1e-18`; preserve complete same-history SEA3 and every actual-applied R_S update. Build the full-21-state finite-`tau_b` cascade/detectability bridge before adding further universal P4 enclosure machinery. P5 remains blocked until strict P4 closes.
