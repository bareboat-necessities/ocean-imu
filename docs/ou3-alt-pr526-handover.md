# OU-III ALT handoff after PR #526

Current continuation: read `docs/ou3-alt-main-handover.md` and
`docs/ou3-alt-startup-pre-rho.md` for the exact entry obstruction and corrected
finite horizon before using this retained implementation handoff.

Read `AGENTS.md` first. This file is the primary continuation handoff after PR #526. Then read `docs/ou3-alt-proof-plan.md`, `docs/ou3-alt-contraction-handover.md`, `docs/ou3-alt-contraction.md`, the ALT section of `docs/ou3-proof-research-state.md`, and `docs/ou3-brmm-main-handover.md`.

## Non-negotiable architecture

ALT is an independent proof track. Preserve the original P2/P3/P4/P5 proof independently; do not delete, weaken, rewrite, or make it depend on ALT. Do not change the shipping filter, quality gates, frozen P3 delta, corrected COMPLETE-BRMM semantics, actual-applied R_S requirements, or zero-wind-heel ALT scope merely to obtain a certificate.

Do **not** begin common/compatible storage or rho search until `proof_plan.assert_finite_storage_master` accepts the actual finite shipping master. `storage_search_allowed=false`, `ALT_STARTUP_PASS=false`, `ALT_LIVE_PASS=false`, and `ALT_END_TO_END_PASS=false` are intentional at this handoff.

## What PR #526 materially closed

PR #526 roots the startup machine at literal reset and carries the executed persistent chain

`guard -> private Mahony -> tracker LPF/stillness -> band/statistics -> WPE -> TuneState -> Racc`

through startup and the same goLive edge instead of seating a synthetic frontend at Live. Guard/source sample counts, physical/bias history identity, machine Mahony vertical ancestry, band input, sigma/tuner state, pending state, runtime configuration, Racc bootstrap identity, scheduler ancestry and a_w-sync goLive re-anchor are checked by construction.

The full dual-compiler binary32 WPE moment/raw-period/log state is now rooted at startup reset, driven by the same private-Mahony vertical signal, preserved through goLive, substituted into the admitted COMPLETE-BRMM/BIAS Live constructor without a new WPE snapshot, and required on every one of the 600 Live IMU edges. MAG and HOLD preserve that WPE state by identity.

The padded COMPLETE-BRMM private-Mahony prerequisite is attached through `admitted_startup_mahony_invariant.py`. It consumes the shared two-phase certificate and requires both `initial_set_inside_invariant` and `continuous_all_live_PI_invariant_closed`, plus the source-order binary32 discrete invariant. This closes the regional admitted-source Mahony startup-to-Live invariant used by ALT. It does **not** close magnetic accumulation-frame accuracy or universal startup reachability.

Both gauged and timeout/ungauged H18 handoff state/covariance maps are represented. The ungauged branch normalizes the persistent proxy quaternion from the same predecessor and installs the free-yaw covariance about its physical down axis. This is exact-real handoff composition; Eigen/binary32 normalization correspondence and proof that every admitted history reaches the branch remain open.

Bounded default scheduler-nextafter and binary64 a_w-sync clock relations are attached. The finite guard also records the raw WPE period branch topology and same-moment ancestry. These are finite-horizon supplies, not indefinite wrapper-clock/counter lifetime proofs.

## Strongest current finite objects

Start with these files rather than older partial products:

- `tools/stability/ou3_alt_contraction/finite_startup_joined_machine_history.py`
- `tools/stability/ou3_alt_contraction/finite_startup_wpe_machine_history.py`
- `tools/stability/ou3_alt_contraction/finite_admitted_startup_wpe_machine_word.py`
- `tools/stability/ou3_alt_contraction/finite_admitted_machine_joined_frontend_racc_interleaved_prefix.py`
- `tools/stability/ou3_alt_contraction/finite_admitted_machine_clock_qualified_interleaved_prefix.py`
- `tools/stability/ou3_alt_contraction/finite_admitted_wpe_machine_clock_interleaved_prefix.py`
- `tools/stability/ou3_alt_contraction/finite_master_guard.py`
- `tools/stability/ou3_alt_contraction/admitted_startup_mahony_invariant.py`

`finite_master_guard.py` is the authoritative fail-closed readiness aggregator. Do not bypass it by invoking an older storage experiment directly.

## Remaining finite-master blockers, in required order

1. **Universal startup branch/reachability.** Qualify the near-antiparallel Eigen `FromTwoVectors` JacobiSVD branch against the actual target Eigen implementation; close nonfinite outcomes; prove the startup/control clocks and timeout path for every admitted corrected COMPLETE-BRMM + BIAS + disturbance history; finish magnetic accumulation-frame/heading qualification where needed by the gauged path. Representation of a branch is not reachability of that branch.

2. **Source-uniform deployment arithmetic.** Finish target/compiler qualification and outward supplies for WPE raw period/log/exp/sqrt, band/statistics/tuner arithmetic, Q-axis operations, scheduler nextafter beyond already bounded pieces, a_w-sync clock math, trig/quaternion normalization, Eigen LDLT/eigensolver/floor decisions, comparisons, and nonfinite branches. Preserve separate and FMA compiler histories where the shipping toolchain permits both. Do not replace these with sampled replay extrema.

3. **Literal same-history 600-transition master.** Every literal event must be qualified under one COMPLETE-BRMM + BIAS0/1/2 + disturbance history: all IMU edges, magnetic calls and their admission/rejection/not-due branches, hold transitions, reset/watchdog/relock branches that belong to the theorem source, covariance floor/sync events, and H18/A21 dimension/control transitions. The existing interleavers are structural attachment, not yet a source-uniform arithmetic certificate.

4. **Only then storage/rho.** Once the finite status truthfully has finite error identity for every event, all coefficient product graphs, and all configured branches, let `assert_finite_storage_master` unlock the next phase. Search common joint24 storage first; move to compatible/path-dependent storage only if the common search is rigorously ruled out. Then prove every-prefix retention and a disturbance-dependent ultimate bound. No replay contraction and no synthetic representative word may promote the theorem.

5. **Indefinite composition.** After finite-word contraction, still prove no-restart tiling and finite-width wrapper-clock/magnetic-counter lifetime or an equivalent exact deployment argument.

## Original-track PE/source issue

There is an inherited original P2/P3/P4/P5 source issue: the declared PE certificate can fail to refine the vector certificate. PR #525 diagnosed the root cause as the vector PE certificate retaining the unpadded 30 deg/s body-rate constant while corrected COMPLETE-BRMM uses the padded 35 deg/s envelope, and reports positive quantitative margin after rebinding to 35 deg/s. At this handoff #525 is separate work. Reconcile that issue on the original track from current main/its eventual merge; do not use it as an ALT shortcut and do not silently import an unmerged theorem assumption.

## Validation/status warning

This remains a proof-progress checkpoint, not an end-to-end PASS. Startup/WPE attachment, source-audit, master-guard and finite native correspondence regressions are separate from universal stability qualification. The original-track PE/vector refinement issue remains separate work; do not turn inherited red checks into ALT PASS flags. Exact validation commands, results and CI-run references are recorded in the PR discussion.

The source-observation audit is bound to the current shipping header. The local `zhat` const reference aliases an already evaluated const `Vector3`, is read once by the magnetic residual, is never mutated and does not escape. The source-audit regression reconstructs the previous complete-file hash after reversing only that copy-to-reference change. A mismatched pin still raises a re-audit error. Regenerated native H18, H18-to-A21, A21 and sample traces agree byte-for-byte across that cleanup; this is finite implementation evidence, not source-uniform deployment arithmetic.

Published replay/results provenance must be regenerated through its canonical evidence workflow after shipping-source changes. No historical replay fingerprint is made current by manually restamping its source hash. Full build validation still requires Eigen and the pinned simulation archive; consult current-main CI rather than assuming every repository check is green.

## Recommended next-conversation prompt

Continue the independent OU-III ALT proof from latest `main`, which includes merged PR #526. Read `AGENTS.md` and `docs/ou3-alt-pr526-handover.md` first and treat them as authoritative. Preserve the original P2/P3/P4/P5 proof independently. Do not start storage/rho search until `assert_finite_storage_master` accepts. First close universal startup branch/reachability and remaining source-uniform deployment arithmetic, then qualify every literal event of the same-history 600-transition COMPLETE-BRMM + BIAS0/1/2 + disturbance word. Keep all ALT PASS gates false unless their actual theorem prerequisites close.
