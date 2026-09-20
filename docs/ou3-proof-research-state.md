# OU-III proof research state

## Current hypothesis

One path in PR #560: construction -> finite capture -> finite H18 bridge ->
reference refinement/release -> recurring informed A21 -> regional practical
stability. All stages carry the same physical execution and estimator state.
The controlling identity is
`M_end' P_end^-1 M_end + D_word = P_root^-1`.
The required uniform inequality is `D_word >= delta P_root^-1`, delta>0,
followed by nonlinear supply and every-prefix retention in the same storage.

## Evidence

- Exact complete-word energy algebra and low-rank factors retain all 21
  columns, actual correction/reset transport and singular nuisance ranges.
  Dense rational adversarial checks verify the algebra, not a uniform margin.
- The corrected Hermite endpoint path gives the singular full-state lower
  comparison `P >= E_LIN A^-1 E_LIN'`, with the full rational 4x4 action,
  interleaved corrections and source polynomial defects retained.
- At recurring regular post-prediction roots after 16 s, combine that full
  comparison with fresh AG/BA process injection using equal convex weights.
  This closes the joint 21-coordinate lower covariance without adding
  unrelated marginal floors (`root_covariance_certificate.py`).
- `ou3-nuisance-upper-proof.md` proves the recurring 15-coordinate LIN/BA
  upper comparison after 17 s in the same default regular A21 profile.
  Three actual S observations eliminate the arbitrary neutral root exactly;
  OU forcing, source Q defects and the actual PSD-sync increments have uniform
  upper bounds. A comparison omitting acc/mag corrections uses the same
  realized sync increments, not a recomputed positive-part map. Block Cauchy
  retains all nuisance cross covariance. Rational constants are reproduced by
  `nuisance_upper_certificate.py`. No full-state decay rate follows yet.
- `ou3-stationary-detectability.md` proves that the stationary A21 lifted
  linear pair has no unobservable nondecaying eigenmode with recurring S,
  acc/mag observations and nonparallel gravity/field. Stable AW/BA components
  vanish at eigenvalue 1; the remaining eigenvector is killed by S and the
  two vector observations. This is not nonlinear quiet-water stability.

## Current limiter

The remaining six-coordinate AG upper comparison must control varying
realized attitude, specific force, gyro transport, correction/reset and
accepted magnetic events. Restricted heading service cannot be embedded as
independent full-state information. The nuisance upper comparison removes
one obstruction but does not close this AG step or the full transported loss.
Finite nonlinear radius, float32 supply, capture, H18/release retention and
every-prefix finite-error retention remain open. The theorem is not claimed.

## Failed approaches / DEAD_ENDS

- Entrywise midpoint-radius Riccati refinement failed twice by dependency:
  prediction enclosure [-18.7907040,2502.41442], failed acc/mag inverse boxes.
  Do not resume subdivision of that mechanism.
- Absolute cross-ceiling/Gershgorin reduction gave gamma=-5.1223e8. An LDL
  pivot is not an eigenvalue floor. Preserve matrix action and factor geometry.
- Endpoint accumulated process noise plus per-event information ceilings
  fails under interleaved corrections; retain the exact 2-D counterexample.
  The corrected endpoint path action charges the actual intervening sequence.
- Restricted heading/bias loss I_2 can coexist with full loss
  `[[1,0,1],[0,1,0],[1,0,1]]`, whose nuisance cancellation vector (1,0,-1)
  preserves energy. This falsifies information lifting, not a shipping marine
  execution. Longer words or tighter scalar constants do not fix the inference.
- A finite quiet-core replay, a QR residual alone, and a literature map did
  not discharge the uniform theorem. Those additions are removed. The exact
  factor identity and nuisance range elimination remain because they support
  the full loss calculation.
- Previous infrastructure failures: literal backslash-n serialization broke
  Python; incorrect Fraction string parsing broke a test; wrong unittest cwd
  broke imports. Fixed at source/invocation. Local TeX lacks IEEEtran/luaotfload
  and package installation lacked setgroups/seteuid permission; use the CI
  renderer. HTTPS push lacks credentials; publish identical Git trees through
  the connector. Oversized tool JSON and urllib artifact HTTP403 were avoided
  by per-file reads and curl. None invalidated a mathematical premise.
- Main refresh 43d41dcc caused 14 generated-evidence merge conflicts. The
  merged branch retains the coherent PR set from 1e18e6ad, not mixed rows.
  At 9fcc76b9 the citation gate rejected inline bibitems invisible to its
  BibTeX audit. The now-removed literature discussion has no remaining cites.
  Recheck the unchanged publication and evidence gates after this cleanup.

- The first cleanup test run failed on the retired motion-domain acronym in
  the new note and on a status JSON not yet regenerated. Use the current
  MARINE MOTION name and regenerate the artifact; neither changes the proof.
  The unchanged gates must pass after those fixes.

- This checkout's required `make all` stopped in `tests/ahrs` while compiling
  `ahrs-qmekf-sim.cpp`: `fatal error: Eigen/Dense: No such file or directory`.
  Neither vendored Eigen nor the system fallback is installed. This is a
  build dependency failure; no shipping source was changed. CI installs
  Eigen and must validate the literal shipping regression builds.

- Publication tree verification caught a dropped executable bit on
  `build_evidence.py`; preserving mode 100755 restored exact local/remote
  tree equality before updating the branch. A short-name fetch selected a
  same-named tag; an explicit `refs/heads/` fetch synchronized the actual PR
  branch. Neither failure changed proof contents or the published branch.

- The local Python quality command initially lacked `ruff`. Installing it
  in a temporary tool directory resolved that dependency; the unchanged
  repository ruff and Python-compilation gates then passed.

## Retained facts

The same-history physical contracts, bounded displacement primitive, exact
bias mismatch, accepted-event service, projection sector, finite-bridge
algebra and Joseph identity remain. No shipping instability is inferred from
a failed sufficient certificate. Quiet water remains admissible: its exclusion
is not justified by the stationary linear analysis. If excitation is needed
for varying histories, require a quantitative window bound; excluding exactly
zero amplitude cannot itself supply a uniform margin.

## Alternatives

Within the same theorem, use full transported loss/Schur elimination or a
joint covariance comparison. The new nuisance upper bound permits treating
nuisance contributions as bounded covariance in an AG trial estimator, but
its response to the actual varying measurement/transport geometry must still
be proved. Do not replace that task with independent principal-block claims.

## Next falsifiable experiment

Derive an AG trial estimator whose root cancellation is uniform for the
admitted realized event geometry. Test its complete loss construction before
rigorous enclosure. If quiet motion is an obstruction, exhibit the full
unobservable mode with the active BA predictor and actual applied events;
then derive a quantitative physical excitation condition that removes it.
The stationary detectability proof rules out claiming that every resting
A21 trajectory lacks information. After AG/full loss closure, compose
nonlinear supply and capture/release-to-retention.
