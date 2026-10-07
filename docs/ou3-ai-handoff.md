# OU-III proof handoff after PR #653

This records the mathematical state being published through PR #653,
`proof/ou3-planar-service-phase`. Publication is not a stability certificate.
Fetch current main, PR state, head and exact-head CI before continuing; after
this PR is merged, start follow-up work from current main.

## Verified starting point

- Preserved automated evidence head: `ca3bd78d255d1e1b1e5eaceadf4c3f405aa50946`.
- Its validation parent: `b92c5e3078a06ae3c8802b309f0172e10601ab03`.
- Incorporated main: `8cecdf39e60e0802182f0c1b5bc9bab74d4cb558`.
- Merge/source-binding commit: `59939745e4770c4677dd4d112033be1258d4a357`.
- The whole 20-file `ca3bd78d` evidence regeneration was preserved. Main's
  BMM150 compensation and calibration-dependent magnetic source selection
  were retained, including the actual OU-III sketch initialization.

The native probe's exact compiler include closure has 34 shipping files,
all byte-identical to immutable `108ce145623f98d4fb6269172ef26d173b9bcb41`.
All eight changed source paths are outside that closure. The existing finite
diagnostic's source binding was refreshed with its earlier audit retained;
its measured results and qualification flags were not changed or replayed.
This establishes source identity for that diagnostic, not a contraction
margin or qualification of the upstream magnetic residuals.

## Mathematical result and one remaining comparison

No new sign result was proved after CR43. Receipt corner decided: **NO**.
Full Schur cost paid: **NO**. Common positive reserve proved: **NO**.
No existing full-word OPEN dependency was discharged. The manifest remains
**25 PROVED, 4 CONDITIONAL, 1 OPEN**; `theorem_closed=false` and
`regional_practical_stability_claimed=false`.

Use the existing appendix's exact CR41 definitions:

\[
K_c=\lambda\bar H-cJ_{r,OO},\quad L=EL_0,\quad
Z_c=L_0K_c^{-1}L_0^T,\quad Y_c=K_c^{-1}L_0^TZ_c^{-1},
\]
\[
S_c=Z_c^{-1}-\lambda E^T\bar D E,\quad
B_c=Y_c^TX_c,\quad
T_c=A_c-X_c^TK_c^\sharp X_c.
\]

The single quantitative obligation is the **feasible**
`R_W(c)=[[S_c,B_c],[B_c',T_c]] >= 0` for one common `c>0` on the
justified activated shipping family. Start by deciding the receipt corner
using the SAME linked action already in `Hbar`; if `S_c>0`, immediately
pay `T_c-B_c'*S_c^-1*B_c` with that same reserve.

CR42 is lossless on justified linear slices. On directional cones retain
the actual inherited image and every prefix guard; an unconstrained
minimizing kernel component is not automatically feasible. Fixed-word
`K_0>0` does not establish a uniform positive-reserve inverse. Do not
silently pseudoinvert unsupported directions or assume receipt rank fixed.

A feasible CR30 mode with `chi>=1` has nonpositive signed odd action.
At threshold one it has `Hbar*y=L'*Dbar*mu`, `mu=L*y` and **positive**
action `y'*Hbar*y=mu'*Dbar*mu`. Zero signed gap therefore does not imply
zero process or correction losses. Neither existence nor exclusion of an
actual such mode is proved. A feasible `S_0*t=0` additionally requires
`B_0'*t=0` for full nonnegativity on a justified linear slice.

The smaller `Dbar` of CR36--CR40 changes `Hbar` on the same reader directions;
it is not independently an improved margin. Do not use the restricted
`10^-37` kernel reserve to price all cross terms. No compact forward domain,
closed normalized inherited image or physical gauge deletion is assumed.

## Focused continuation

Read `AGENTS.md` and `docs/ou3-shipping-faithfulness-protocol.md` first,
then the current research ledger, relevant status/manifest entries and
`doc/kalman_ou_iii/kalman_ou-w3d-stability-proofs.tex-part`, especially
CR29--CR43 and LR8/LR9. The status record is
`reports/results/ou3_stability/theorem-status.json`.

Use the literal integrated OU `g/Q_ao`, queued-target lag, previous AA
deletions, same-cycle acc and transported actual S/magnetic covariance
columns on ONE inherited tangent. Those losses are already counted in
`Hbar`. Physical magnetic service does not supply an unproved generic
covariance Gram floor. Keep metric cross entries, generated suffix scores,
`E_L+S_q-2C_W`, generated process cross/square work, reset/projection,
inherited auxiliary/reference dependence and both endpoint gauge terms.

The next result must cross the existing feasible receipt threshold and
the full Schur threshold, or quantitatively isolate their unpaid linked
comparison. Another equivalent formula, zero-action theorem, storage,
observer, coordinate reset, independent tuner box or empirical search is
not that result. Do not begin source completion or radius work beforehand.

## Validation and shipping-faithfulness record

Exact-head GitHub results and the final merge commit are publication
metadata and must be fetched again. At `59939745`, both push and PR
theorem/planar/moving checks and the PR evidence contract passed. The local
source-bound proof validator passed with no failures and no theorem promotion.
The primary required local validation is `make all`. Its final result and
the remaining exact-head CI belong to the PR's publication metadata.

1. **Structures preserved:** all 21 states, actual covariance/gain/Joseph,
   integrated OU/S chronology, AW floors/queued targets, applied acc/magnetic
   action, resets/projections, causal tuner/reference/auxiliary histories and
   endpoint gauge; upstream main behavior is retained without proof edits.
2. **Relaxations introduced:** none by merge/handoff preparation. Existing
   real planar partial-fibre, inverse, inherited-image and cone qualifications
   remain; no reachability or forward retention is assumed.
3. **Failures:** E for automated-head workflow admission and stale full-src
   bindings; repaired by preserving evidence and auditing exact source identity.
   An incorrect status-document lookup was corrected to the JSON path above.
   Earlier D-class sufficient-proof failures remain classified in the ledger.
4. **Genuine shipping counterexample:** none. No actual negative feasible
   eigenmode or lossless shipping impossibility was established.
5. **Established results retained:** all 25 proved entries, including the
   integrated process allocation, nonzero receipt accounting and CR36--CR43
   identities under their stated qualifications. No margin or flag promoted.
6. **Next calculation:** decide the SAME feasible `R_W(c)` using its linked
   shipping action, first the receipt corner and immediately its full Schur
   cost. Endogenous inherited work remains endogenous.
