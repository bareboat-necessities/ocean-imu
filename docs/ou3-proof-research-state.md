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

## Retained finite-`tau_b` premise

`tools/stability/ou3_sea3_a21_detectability_completion.py` establishes the
paper-level finite-bias detectability/UES hypothesis from:

- complete-SEA3 H18 contraction;
- exact finite residual-bias Gauss-Markov decay;
- bounded full-state H18<->b_a coupling on the compact word;
- full A21 process UCC;
- no alternate estimator and no state elimination.

Its stronger implementation-word flags remain deliberately false. The joint
master consumes this result only through the validated canonical P3 chain; it
does not promote the comparison observer or replace the full A21 matrix test.

## Current hypothesis

Use the exact complete-word endpoint identity and the full 21-state finite-`tau_b`
P3/detectability result to build one correlated nonlinear graph-sector master.
For `z=[x;w_W]`, the controlling matrix is
`L_W=[[-D_W,M_W^T J_N B_W],[B_W^T J_N M_W,B_W^T J_N B_W]]`.
Admissible graph sectors `z^T Pi_j z>=0` enter only through the full
S-procedure test `-(L_W+sum lambda_j Pi_j)>0`; the same construction is
required at every prefix.

## Evidence and current limiter

P3 was recomputed with the deployed `0.4 m/s^2` accelerometer-bias projection
limit: H18 and A21 retain `delta=1e-18`, the first active A21 bias full-matrix
margin is `1.2499987189052501e-9`, and the H18 worst interval LDLT pivot is
`4.987499868870966e-14`. `P3_DEPLOYMENT_PASS` remains false only for the
separate physical-language inclusion obligation.

The canonical 6 s A21 payload has a strict small-error margin
(`rho=0.9911176` at scale `0.125`) but crosses one at scale `7.5` and
reaches `1.1357229` at the retained-domain boundary. The finite-`tau_b`
detectability rerun passes with bias energy gap `1.1992803e-3` and
asymptotic A21 gap `1e-18`. Thus linear bias decay is not the limiter;
the limiter is the correlated attitude/`a_w`/`b_a` accelerometer curvature.

## Failed approaches / DEAD_ENDS

On the same canonical payload, diagonal bias precision, fixed-frame and
source-framed `a_w`/`b_a` cross penalties, and separate latent diagonal
energies all failed after their allowed refinement. The closest result was
`rho_linear=0.996797290`, `rho_finite=1.00009035805` at
`(beta_aw,beta_b)=(500,310000)`. Do not resume metric-grid tuning.

## Retained facts, alternatives, and next experiment

Retain the complete same-history source, frozen P3 `delta=1e-18`, full H18/A21
states and cross terms, exact Cayley/reset residuals, all valid accelerometer
and vector events, and every due S event with actual applied `R_S` inside its
suffix. Packetwise radii, state elimination, replay fitting, and further
domain/filter changes remain forbidden.

The remaining alternatives are a dense source-structured storage LMI, a
path-dependent joint storage, or direct falsification of a source-uniform
inner funnel. The next falsifiable experiment is to materialize only the
same-history state/residual selectors needed by the master, lift the exact
Cayley and A21 projection bounds into dense graph sectors, and evaluate the
full augmented matrix first on the canonical point word. Proceed to outward
source-uniform LDLT only if that non-promoting matrix diagnostic is below one.

## Complete-source obligation still open

The upstream complete-source contract still does not materialize the correlated finite-window SEA3 realization as an executable outward source family. Parameter compactness, RAO/moment envelopes, hard pathwise acceleration/body-rate caps, frontend parity, and adaptive-state rate/jump bounds do not replace a same-history transition for the correlated source state. P4 may not substitute independent per-sample boxes, a replay record, or a finite harmonic/grid surrogate.

The full-state master now identifies the needed quantities: same-history prefix
state selectors, stacked nonlinear residual selectors, reset/boundary terms,
and the A21 projection graph. Extend the existing complete-SEA3 execution only
to materialize those quantities, rather than creating another source language.

## Current CI / evidence status

The exact `source-foundation` command passes 73 tests. The 0.4 P3 recomputation,
finite-`tau_b` detectability rerun, canonical P4 integration, and all 91 P4
discovery tests pass. `make all` is locally blocked before compilation because
this environment has neither `/usr/include/eigen3/Eigen/Dense` nor a vendored
Eigen tree; the build command and include policy were not changed.

`ou-validation` is red for a known evidence-provenance reason, not because its numerical unit-test body found a new filter failure: `tools/ou_evidence_contract.py --auto` reports replay dependencies changed relative to committed validation/robustness provenance, including the OU-III filter and WavePeriodEstimator dependencies. Genuine validation/robustness evidence regeneration is therefore still required before a later proof PR is declared final/ready. Do not hand-edit provenance hashes.

No claim is made that P4 or P5 is complete.

## Promotion boundary

The master and graph-sector machinery are non-promoting. P4 remains open until
the source-uniform endpoint and every-prefix augmented LDLT plus domain
retention close on the same complete SEA3 history. P5 remains blocked until
strict canonical P4 contraction closes.

## PR #500 experiment 1 failure analysis — arbitrary two-prefix covariance fixture

Hypothesis: the trusted typed complete-SEA3 execution kernel can be instrumented passively to retain branch-correlated every-prefix ancestry and event-local H18/A21 Riccati cells (P-before/P-after plus the exact F/Q, floor increment, or Joseph H/R) without reimplementing shipping transitions or losing actual-applied R_S provenance.

Execution: PR #500 added passive event capture inside `ou3_sea3_complete_window_execution_kernel.advance_branch` and a selector smoke that executed two identical point samples from the synthetic covariance fixture `P0_H=2 I_18`, `P0_A=2 I_21`. The canonical `ou3-proof` run `34132984840` passed source-foundation, complete-SEA3 source, frozen `riccati-p3`, the P4 geometry/metric/reset gates, and all pre-existing P4 tests before reaching the new selector test.

Observed failure: `test_ou3_p4_complete_sea3_same_history_prefix_selectors.P4CompleteSea3SameHistoryPrefixSelectorsTest.setUpClass` failed during the second synthetic prefix in the unchanged outward interval Joseph backend. `matrix_inverse_gauss_jordan` rejected innovation pivot 1 because its enclosure crossed zero: `[-31.51349023865584, 33.33589453040376]`. The failure is an interval-enclosure failure of the arbitrary `2 I` multi-sample smoke fixture; it is not a P4 contraction result and not a shipping filter failure.

Interpretation: using `2 I` as a convenient multi-prefix covariance is unjustified and substantially wider than the source-generated Normal-Live covariance structure. Refining the Gauss-Jordan backend, selecting a favorable interval branch, dropping the second prefix, or weakening the Joseph inversion would hide the fixture problem and is rejected. The typed capture concept itself was not falsified: the first sample and all pre-existing exact differential-event tests reach the same P/H/R semantics successfully.

Replan: replace only the arbitrary smoke covariance with the canonical source-generated Normal-Live H18/A21 seed structure (including the actual attitude handoff seed, v/p/S seeds, committed stationary a_w covariance, and the A21 bias-release seed where applicable), while keeping the same trusted transition, every front-end successor, exact event order, actual applied R_S, and unchanged outward inverse backend. If the source seed still loses a pivot on the second prefix, stop and treat that as a separate interval-representation obstacle rather than widening/regularizing the inverse by hand.
