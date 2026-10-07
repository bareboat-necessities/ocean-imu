# OU-III receipt continuation from merged main

## Baseline and separate CI repair

PR #653 is merged. Verified starting main: `a47da994b85f2ec844b199d27b619a7fcc553adb`.
Follow-up branch: `proof/ou3-receipt-schur-continuation`; never reuse the old branch.
Mechanical repair commit: `9b951b6ac`.
Concurrent main `59171c1a4aef654bb9808629e8856e4a56edfcc8` is incorporated,
including PR #663, version 2.4.0 and its automatic validation/robustness bundle.

Exact main reproduction rejected frontend binding and the finite secant whole-src
binding. The 34-file native compiler closure is unchanged in membership: 33 files
are byte-identical; Mahony changes only a constant declaration's whitespace.
Exact-profile preprocessing differs only by that same line-wrap replacement.
The two refreshed bindings retain original payloads, qualifications and all
previous audits. Validators and fail-closed checks are unchanged. Upstream
AtomS3R calibration, compensation and source selection remain. The separate
warnings job's uninitialized test outputs are value-initialized; its compensation
regression passes with `-O2 -Wall -Wextra -Wpedantic -Werror`.
Concurrent main restores the earlier declaration wrapping and initializes the
same failing buffers. The final 34-file closure and complete preprocessing
match the earlier valid baseline exactly. Rebind only the affected source
fields and retain both repair audits; no upstream behavior is reverted.

## New inequality, with its limitation

The existing appendix now contains CR44--CR45, implemented in the existing
`measurement_frame.py`. The target is still CR39 followed by the complete
CR41--CR43 remainder; no new architecture or storage is introduced.

On the actual preceding fixed-target covariance prefix use
`X=F*P*F'`, `D=F*dP*F'`, `Qr=Q_ship-q_aw*e_a*e_a'`,
`v=X_aa`, `w=X*e_a`, `d=w'*Qr^-1*w`, and `mu=D_aa`.
The SAME retained full residual process has the quantitative payment

    L_rem >= alpha*mu^2,
    alpha=(4*d+v)/(v*(2*d+v)^2)>0.

Its exact remainder retains a nonnegative transverse Fisher square and the
full process addition `Qr-w*w'/(2*d) >= Qr/2 > 0`. This prices all nonzero
receipts simultaneously on their common inherited tangent, with full integrated
Q_ao, AG/BG/BA, polynomial/PSD qualifications and every earlier AA deletion.
The process is allocated once. Paired acc squares and actual S/magnetic action
stay in Hbar; they are not counted again. Queued and committed targets remain
distinct. Scalar zero-gap guards and inherited-image constraints remain.

**This does not decide the receipt sign.** Alpha has no certified uniform lower
bound or proved comparison with dbar on the activated shipping family. For this
particular allocation, with `k=dbar*v^2`, complete individual receipt payment
requires `0<k<1` and

    d/v <= sqrt(1-k)/(2*(1-sqrt(1-k))).

For `k>=1` that allocation alone cannot pay it; the retained joint word action
may still do so. This is not a necessary word premise or a negative shipping
mode. The exact source-domain quantity is one inherited covariance-column energy:

    d/v=(P*e_a)'*F'*Qr^-1*F*(P*e_a)/P_aa.

It must be controlled from the SAME coupled tuner/process/acc/S/magnetic
history, not independently boxed. No actual feasible chi>=1 mode has been
exhibited or excluded. Positive process action at zero signed gap remains
possible; no zero-action argument or compactness shortcut is used.

## Full endogenous cost and status

Keep original CR41--CR43:
`K_c=lambda*Hbar-c*J_r_OO`, `L=E*L0`,
`S_c=(L0*K_c^-1*L0')^-1-lambda*E'*Dbar*E`,
`B_c=Y_c'*X_c`, `T_c=A_c-X_c'*K_c_sharp*X_c`.
The required feasible full remainder is `R_W(c)>=0` with ONE common `c>0`.
Where its ordinary inverses and `S_c>0` are justified, pay
`T_c-B_c'*S_c^-1*B_c` immediately. On cones retain every actual guard and
minimizing-kernel feasibility. No receipt rank or inherited origin graph is assumed.

Metric cross entries, generated covariance suffixes, dK*r,
E_L+S_q-2C_W, reset/projection, inherited auxiliary/reference correlations and
both endpoint gauge terms remain endogenous. No cross cost is priced with the
restricted 10^-37 reserve. Source completion must wait for a homogeneous margin.

- Receipt blocker resolved: **NO**.
- Full endogenous Schur cost paid: **NO**.
- Complete homogeneous gap proved: **NO**.
- Existing OPEN dependency discharged: **none**.
- Manifest: **25 PROVED / 4 CONDITIONAL / 1 OPEN**.
- Uniform reserve: **not established**; theorem/regional/capture flags unchanged.

The AtomS3R OU-III sketch still disables continuous hard-iron learning after
refinement. Its inherited refined-reference error and source correlation remain;
this does not qualify magnetic residual envelopes. Its distinct noise/S-cadence
profile is not inferred from the planar facade/native profile.

## Validation and six-item shipping record

Final make/evidence/appendix results, published SHA and exact-head/merge-tree CI
belong in the new PR metadata and must be fetched again. Parent passes do not
certify a new head. Only required source/evidence regressions were run; no
exploratory replay, sampled contraction search or secant sweep was added.

1. Preserved: all 21 states/P/K/Joseph/masks, full integrated OU/BG/BA/S,
   floors/queued lag, scheduler, actual corrections, tuner/reference/source
   inheritance, reset/projection and endpoint gauges; no proof-boundary reseed.
2. Relaxations: none in CR44--CR45; existing qualified real planar partial
   fibre and actual scalar cones remain. Regression slot operands do not
   certify reachability, admission or uniformity.
3. Failures: E for stale bindings, warning buffers and local dependency/path
   handling; repaired without weakening checks. Parallel make all started
   tests before binaries existed. The first serial command exited 2 with no
   final diagnostic in its saved log (last line: passing TFG Jacobians);
   repeat with live capture on the integrated tree. D if pointwise alpha>0 is
   promoted to uniform receipt domination: no threshold-crossing bound exists.
4. Genuine shipping counterexample: none; no A/B conclusion is asserted.
5. Retained: all 25 existing proved entries and qualifications, CR25--CR43;
   plus the nonzero-receipt lower inequality, without margin/count promotion.
6. Next: quantitatively control the actual covariance-column energy and the
   remaining joint process/acc/S/magnetic receipt action, excluding ALL feasible
   CR30 modes with chi>=1. Then pay the SAME full CR43 cost. Another equivalent
   completion or positive-substep assertion is not a sign certificate.
