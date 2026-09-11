# OU-III ALT post-merge handover

Resume from latest `main` after PR #522. Read, in order:

1. `AGENTS.md`
2. `docs/ou3-alt-proof-plan.md`
3. `docs/ou3-alt-finite-measurement-proof.md`
4. `docs/ou3-alt-contraction.md`
5. the ALT section of `docs/ou3-proof-research-state.md`
6. `docs/ou3-brmm-main-handover.md`

The original P2/P3/P4/P5 proof remains an independently continuable route. ALT must not rewrite it, consume its PASS labels as ALT proof, or import its PE/vector mismatch unless ALT genuinely needs that lemma.

## What is genuinely proved on ALT

- Exact finite accepted-measurement algebra on joint24, including finite Cayley reset, exact secant residuals, physical `r_S=e_S-S_phys`, actual H18/A21 gain masking, radial bias projection, and thin rank-3 storage contribution without state reduction.
- Exact H18 held-covariance invariant and conditional H18->A21 covariance map: BA cross blocks stay zero, the hidden BA marginal stays at the configured seed variance, and the enable-floor is a fixed point on that invariant.
- Separate analytic BIAS0/1/2 contracts with one persistent physical history/root token and one shared driver entering both `e_ba` and `beta`.
- Exact finite prediction algebra for attitude, translation, and physical accelerometer bias. Translation uses the existing correlated q15 `(a0,a1,J0,J1,J2)` forcing relation.
- The deployed quaternion predictor is source-uniformly inside the polynomial branch on the declared regional Live domain: the worst shadow increment is about `0.00311 rad < 0.01 rad`.
- Regional Live/source lineage machinery is separated from startup capture. Startup remains a separate theorem.
- Anti-dead-end guards are executable: replay, unreachable perturbations, Jacobian cocycles, incomplete source tokens, and premature common-metric/rho searches are blocked from promotion.

## What is NOT proved

- The complete finite 600-step physical word is still OPEN. `physical_word.py` currently composes pointwise Jacobians/source responses and is explicitly not a finite endpoint identity.
- The exact finite measurement and prediction descriptors have not yet been composed through all literal covariance/frontend successors, due-S events, asynchronous magnetometer events, H18/A21 branches, and the H18->A21 splice into one source-uniform finite word.
- No common joint24 storage, source-uniform rho, useful ultimate bound, every-prefix chart-retention certificate, deployment finite-precision enclosure, startup capture, or end-to-end theorem exists yet.
- `ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
- Original P4/P5 remain unpromoted by ALT.

## Next theorem obligation — do this first

Build a **finite physical word composer**, replacing the remaining derivative/Jacobian word representation. It must compose the exact finite prediction and measurement descriptors along the already-admitted same-history source relation while retaining:

1. one BRMM generator and one Live S origin over the full word;
2. the correlated q15 forcing relation, never independent per-step boxes;
3. one BIAS-family root/history token and the same driver in truth/error recurrences;
4. actual covariance/frontend/tuner successors and actual anisotropic `R_S`;
5. H18 held covariance with BA uncertainty retained in innovation, then the exact H18->A21 covariance splice when the literal release guard fires;
6. every configured due-S and accelerometer event and all admitted asynchronous magnetometer branches;
7. finite Cayley reset/projection branch graphs and chart denominators;
8. physical reference forcing, especially `S_phys` and continuous-physical-vs-sampled-gyro attitude forcing.

Only when the finite-word guard closes may the first common coercive joint24 storage search run. Search one common `M` first. Parameter-dependent or piecewise storage is allowed only after a genuine falsified common-M attempt on the complete finite master. Rank-3 structure is an exact arithmetic optimization only.

## Dead-end guards to preserve

Do not:

- use replay or finite seeds as theorem evidence;
- perturb unreachable covariance/frontend roots and call that physical admission;
- freeze `K`, `P`, frontend, tuner, or guard state;
- replace the finite map by products of Jacobians;
- Cartesianize correlated source variables or bias history;
- reset `S` per word or re-anchor physical displacement;
- marginalize A21 to motion-only storage;
- import the old PE/vector mismatch unless ALT actually requires that lemma;
- shrink the declared source/error domain or lower the frozen P3 delta;
- run common-metric/rho/high-precision refinement before the complete finite physical word exists.

## Handover validation state

PR #522 was mergeable at handover preparation. Its head immediately before adding this file was `8c5f460c170ae92dbebac87ea42dda117046f026`; this handover commit advances that head. CI for the preceding head had started but was still pending/queued, so merging must not be interpreted as an all-CI-green claim. Earlier focused finite-map tests passed; broader-suite/environment limitations are documented in PR #522. The shipping filter itself is unchanged by ALT proof tooling.

For the next conversation, start a new PR from latest `main` and continue only the finite-word obligation above. Do not resurrect superseded PRs #519/#520/#521 or their replay/finite-perturbation diagnostics as theorem evidence.
