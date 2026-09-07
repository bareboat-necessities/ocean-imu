# OU-III proof research state

## Current checkpoint

Canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional complete-SEA3 P3 is **closed and frozen** at `delta=1e-18`. P4 is **OPEN** and P5 is **BLOCKED**. H18 and A21 are both required; H18->A21 remains a separate rectangular hybrid event. Zero lever arm and the dormant-transparent vibration branch remain the certified production branch.

PR #496 is merged at `75f8ecc58bec654bb994ebcb7d74921355514544`. Continuation PR #498 is the active research branch for the paper-permitted finite-`tau_b`, full-21-state A21 route. No merge is authorized.

The only production/proof-domain change inherited from #496 is the user-authorized accelerometer-bias projection-radius tightening from `0.5` to `0.4 m/s^2`. The declared startup/handoff accelerometer-bias error envelope is `0.4 m/s^2`; the Normal-Live physical active-bias bound is `0.35 m/s^2`, leaving a strict `0.05 m/s^2` interior margin. No further filter tuning, source-language parameter change, quality-gate relaxation, or P3 change is permitted by this research route.

The P4 target remains the finite full-state source-indexed quadratic storage

`V(e,zeta)=e^T M(zeta)e`

on one complete same-history SEA3 word, with a strict finite endpoint inequality, a finite every-prefix gain, and every-prefix chart/source-domain retention. Point experiments, replay-fitted metrics, reduced-state coercivity, alternate source languages, selected-S words, independent tuner/R_S boxes, state elimination, packet-count nonlinear budgets, scalar correction radii, or marginal inverse-metric-floor arguments cannot promote P4.

## Exact structural work retained

The retained theorem/execution machinery includes:

- exact finite accelerometer coordinate shift with original shipping H/P/K/S and `H0 != H_u` at finite attitude error;
- exact prediction transport with literal full `F E_aw`, including v/p/S/a_w rows;
- exact Joseph/reset signed-information identity
  `Delta V = -I_y + E_eta + X_reset + E_reset`;
- every due S=0 event with the **actual applied anisotropic SpectralMSE R_S**; S has `eta=0` exactly and therefore contributes favorable information with no nonlinear charge;
- every valid Normal-Live accelerometer update, applicable vector update, full Q, covariance floors and immediate resets;
- separate H18->A21 hybrid lift;
- exact deployed Cayley reset transport and generalized A21 bias-projection branch;
- outward interval/differential AD and generalized mean-value machinery for the exact nonlinear physical map;
- complete-word accelerometer covariance channel and homogeneous Cayley residual-sector factorization;
- finite-`tau_b` A21 detectability/UES support tied to complete SEA3.

Do not replace these with replay families, independent sample boxes, finite harmonic/grid surrogates, selected-S words, packet-count-times-worst-remainder bounds, correction balls, inverse-metric floors, `a_w` elimination, or any further proof-driven filter/domain tightening.

## Authoritative point evidence from #496

The corrected single shipping observer owns one source/tuner/Riccati history and retains every operation. The genuine PM+Stokes Hs=1.5 m point diagnostics are falsification/debugging evidence only.

Pre-0.4 legal 3-second linear words:

- H18: `rho_linear=0.9998658024147671`, 600 predictions, 600 accelerometer updates, 137 actual-R_S S updates, 75 vector updates;
- A21: `rho_linear=0.9958536807113242`, 600 predictions, 600 accelerometer updates, 108 actual-R_S S updates, 75 vector updates.

On the 0.4 head the 3-second physical reset-normalized diagnostic keeps zero-state parity. H18 contracts on every retained tested scale. A21 first crosses one at scale `8.0`; the worst retained-domain endpoint ratio is about `1.0860152320`. The 0.4 clamp is therefore a legitimate production/domain tightening, not a proof fix.

The source-contiguous 6/9-second physical diagnostic is GitHub Actions run `34083473117`, artifact `10004621877`, SHA256 `8b3c1e3b83e6774895cf20770764609afd38e044131204c568d49d7f3544aa19`. It uses one shipping observer, complete event order, actual R_S and zero-state parity and is explicitly non-promoting.

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

The 3->6->9 longer-window tactic is stopped. Do not optimize window length further around this replay.

## A21 mechanism

The finite A21 obstruction is concentrated in a coupled attitude / latent-acceleration / accelerometer-bias cancellation direction. With lever arm disabled,

`y = (E-I) f_hat + E R_hat delta_a_w + delta_b_a`,

while the first-order row is

`H e = [c]_x f_hat + R_hat delta_a_w + delta_b_a`.

The first-order pieces can cancel while second-order attitude curvature remains. The genuine #496 long-window maximizing directions are strongly accelerometer-bias dominated in Euclidean direction energy: approximately **84.5% b_a at 6 s** and **89.0% b_a at 9 s**. This supports the finite-`tau_b` route, but does not justify a bias-only metric.

The nonlinear accelerometer residual has **no nonlinear b_a term**: b_a enters the physical residual exactly linearly. Thus finite-`tau_b` should be used to break the temporal `b_a`/`a_w`/attitude cancellation through the full same-history word, not charged as an artificial nonlinear bias remainder.

## Finite projection lemma closed on #498

The shipping 0.4 m/s^2 accelerometer-bias projection is the Euclidean projection onto the closed ball `B_2(0,0.4)`. The theorem-domain true active bias satisfies `||b_true||<=0.35`, hence lies strictly inside the ball. Orthogonal projection onto a closed convex set is nonexpansive relative to every point in the set:

`||Pi_B(x)-b_true||_2 <= ||x-b_true||_2`.

Therefore projection cannot increase finite physical bias-error norm. P4 does **not** need an adverse projection penalty, a correction-radius bound, or an assumption that the sampled projection branch is inactive. This is finite-error, not tangent-only.

## Fixed additive bias metric: falsified and retired

A natural diagnostic candidate was

`M_alpha(P) = P^{-1} + alpha E_b^T E_b`,

with one fixed positive bias weight `alpha` and no state elimination. The genuine #496 A21 data falsify this tactic.

For the selected worst nonlinear physical words, endpoint contraction requires approximately:

- 6 s: `alpha >= 115346.17`;
- 9 s: `alpha >= 86656.74`.

But scanning **all legal A21 linear words already present in the same genuine artifacts** shows the largest fixed `alpha` preserving linear `rho<1` is only approximately:

- 6 s: `alpha < 46392.08` (limiting legal word near `t0=1176.9927 s`);
- 9 s: `alpha < 59331.27` (limiting legal word near `t0=125.9704 s`).

At `alpha=100000`, legal linear words reach roughly `rho=1.10785` (6 s) and `rho=1.11385` (9 s). The nonlinear-repair lower bound and linear-contraction upper bound do not overlap. Per the research two-strike rule, **stop fixed additive bias weighting**. Do not retune alpha or replace it by a fixed diagonal bias congruence.

This is a proof-method failure, not a counterexample to the paper-permitted source-indexed full-cross-term metric family.

## Quantitative finite-tau_b comparison cascade on #498

The existing `ou3_sea3_a21_detectability_completion.py` closes the paper-level finite-bias detectability/UES hypothesis using a triangular comparison observer

`E_A = [[E_H,C_Hb],[0,Phi_b]]`

but its old `N*K*L^N` C_Hb number is only a finiteness diagnostic and may not be consumed by P4.

#498 now formulates the full-rank comparison storage

`M_A(zeta)=diag(M_H(zeta),mu I3)`

with target comparison gap `delta_A=5e-19` while leaving frozen P3 at `delta_H=1e-18`. If `c^2 >= ||M_H,+^(1/2) C_Hb||_2^2`, the exact two-block Schur condition is

`c^2/mu < ((delta_H-delta_A)(beta_b-delta_A))/(1-delta_A)`.

Strict gaps are stored directly; no binary64 `1-delta` comparison is used.

The comparison C_Hb bound is now expressed without a suffix-norm product or packet-count nonlinear budget. Choose the H18 comparison observer with zero b_a correction row. Starting from `h_0=0`, the only direct linear b_a input to H18 is the accepted accelerometer residual. The exact Joseph identity gives

`Delta V_H <= s_H b_j^T R_acc^{-1} b_j`,

and the comparison bias follows exact homogeneous GM decay. Therefore

`||C_Hb b_0||^2_{M_H,+} <= s_H lambda_max(R_acc^{-1}) sum_{j=1}^N exp(-2 j h_-/tau_b) ||b_0||^2`.

The geometric sum is evaluated in closed form with outward-safe `expm1` arithmetic. Every valid accelerometer update remains present; S/vector/Q/floor/reset operations remain inside the complete H18 information word, including actual applied R_S. The resulting comparison-cascade Schur weight is represented logarithmically if very large rather than overflowed.

This closes only the **triangular comparison-observer** quantitative cascade if its tests pass. It does **not** close the actual A21 Kalman word because the comparison observer deliberately has a zero bias correction row.

## Actual A21 shipping bridge: current controlling obligation

The next theorem object must retain the actual A21 Kalman bias correction row and all H18<->b_a covariance/metric cross terms. The fixed additive bias metric is ruled out, and the comparison observer cannot be relabeled as the shipping estimator.

The preferred next route is a source-indexed **full-cross-term 21x21** word inequality / conditional-information construction that uses finite `tau_b` as temporal separation inside the complete same-history word:

- b_a has exact finite GM decay;
- b_a enters accelerometer residual linearly;
- a_w/translation are tied to the complete S-chain, with every actual-R_S S event supplying favorable information and zero nonlinear charge;
- attitude/gyro-bias are tied to recurring vector information;
- all A21 Joseph cross terms and bias-row corrections remain in the joint matrix;
- nonlinear accelerometer/vector/reset costs are compared against the joint signed information, not an event-count scalarization.

The open bridge must either produce a source-uniform full 21x21 metric/information inequality for the **actual** A21 homogeneous word or expose a genuine admissible expanding direction. It may not shrink the domain or invent a weaker certificate.

Only after this actual full-cross-term A21 linear bridge is quantitative and outward-certified should the universal nonlinear complete-word endpoint/prefix enclosure resume.

## Complete-source obligation still open

The upstream complete-source contract still does not materialize the correlated finite-window SEA3 realization as an executable outward source family. Parameter compactness, RAO/moment envelopes, hard pathwise acceleration/body-rate caps, frontend parity, and adaptive-state rate/jump bounds do not replace a same-history transition for the correlated source state. P4 may not substitute independent per-sample boxes, a replay record, or a finite harmonic/grid surrogate.

The continuation should request from the complete-SEA3 machinery only the source-correlated quantities needed by the full-cross-term A21 inequality, rather than creating another proof language.

## Retired / forbidden routes

Do not revisit:

- estimator-pair shadow promoted as theorem map;
- raw/full-Phi endpoint optimization after their A21 finite-scale failures;
- arbitrary single-map converse-metric fitting;
- longer replay windows merely to seek a green point;
- fixed additive/diagonal bias metric after the no-overlap falsification above;
- reduced-state/Schur certificates or a_w elimination;
- selected-S words or independent R_S/tuner schedules;
- scalar Lipschitz, correction-radius, inverse-metric-floor, or packet-count-times-worst-remainder bounds;
- further proof-driven filter/domain tightening beyond the authorized 0.4 change.

## CI / evidence

At #496 handoff the genuine long-window diagnostic was green. Genuine validation/robustness replay evidence still requires regeneration before a future final/ready proof PR because committed provenance is stale after proof/filter dependency changes. Do not hand-edit provenance hashes.

On #498, normal build passed on early heads. The new finite-`tau_b` bridge/projection/GM-energy modules require canonical proof/quality validation on the final head. `ou-validation` may remain red for the known stale evidence-provenance reason until genuine evidence is regenerated; numerical failures must still be investigated separately.

No claim is made that P4 or P5 is complete.

## Continuation instruction

Continue on PR #498. Keep P3 frozen at `delta=1e-18`; preserve complete same-history SEA3 and every actual-applied R_S update. Validate the new comparison-cascade modules, but do not promote their comparison-observer metric to P4. Build the **actual A21 full-cross-term finite-tau_b shipping bridge** next. If it closes, attach the signed nonlinear complete-word sector and certify endpoint/prefix/domain retention. If it fails on a genuine admissible direction, record the obstruction directly. P5 remains blocked until strict P4 closes. No merge without explicit authorization.
