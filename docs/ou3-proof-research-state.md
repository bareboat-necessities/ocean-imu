# OU-III proof research state

## Current hypothesis

The current point-sample domain **does not imply the declared six-degree
physical capture**, even with all-time actual magnetic service. The requested
construction -> capture -> H18 -> release -> A21 theorem cannot be completed
as stated. No motion/sensor assumption, estimator behavior or threshold has
been changed. Conditional local A21 stability is not refuted.

## Evidence

`ou3-sampled-capture-obstruction.md` proves an all-time counterexample, with
exact rational envelope and service constants in
`sampled-capture-obstruction.json`. Two nonzero sinusoidal translations have
constant true rolls +/-acos(3/5), zero bias/noise, and identical stationary
200 Hz IMU / 25 Hz magnetic samples. Position, velocity, acceleration and the
global displacement primitive satisfy all declared bounds. The actual wrapper
completes startup, refinement and A21 release and stays level. Physical tilt
error is always 53.13 degrees, so no finite six-degree capture time exists.

The invariant source heading/bias block admits an exact all-prefix covariance
upper comparison. Root observation information, with all process correlations
and actual sequential innovations retained, exceeds 3.72037 per one-second
window; even injection about the true down axis exceeds 1.33933 > mu_M=1.
This is a real-arithmetic all-time service proof, not inference from packet
cadence or a finite replay. The native `sampled_capture-test` separately checks
the literal float startup/release path, applied corrections, level invariance
and axis separation. Exact sequential-vs-batch root-information tests include
interleaved process noise; they do not reuse the invalid endpoint comparison.

Validation: the clean-checkout evidence/publication suite passes 449 tests
with one existing data-dependent skip; exact artifact reproduction, ruff,
Python compilation and all three literal shipping regressions pass. The
sequential full `make all` compiled every target and passed the native suites
through `tests/spike_filter`, then was interrupted in the repository-wide
architecture text scan of generated simulation files. The identical validation
suite passes in the clean checkout; its skipped record-convention test passes
separately against the fetched data. The remaining wave-direction and wave-sim
suites also pass separately. The interrupted aggregate command is not reported
as passed.

## Current limiter

Classification: mathematical information obstruction to universal physical
capture. Invalidated hypothesis: the current three physical contracts and
point-sample sensor bounds force entry into the declared six-degree region.
The witness defeats every finite deadline, including a history-dependent one.
It does not demonstrate covariance divergence or failure of conditional local
stability. Removing quiet water cannot repair it: both histories have motion.

The missing connection is quantitative sampling fidelity between physical
acceleration and velocity increments. The proof derives a concrete cell-error
condition and its implied sampled force diversity, without adopting it. Real
sensor integration or qualified physical bandwidth/jerk could provide such a
condition. The current formal premises provide neither.

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

- The continuous velocity-balance force-diversity identity is valid, but
  transferring it to point samples without a quantitative cell-error bound
  is invalid. The exact nonzero histories above falsify physical capture
  even when the nominal sampled A21 linearization is detectable. Further
  covariance tightening cannot repair that implication.
- Audit of the service induction found that a handoff can occur inside a
  magnetic cell. A post-mag invariant alone does not initialize that induction.
  The exact seed covariance's first partial-cell prediction is now certified
  below U before the first correction, for every phase; the subsequent
  periodic induction and loss floor are unchanged. No clock coincidence is
  assumed and no service hypothesis is weakened.
- Local dependency recovery: `apt-get download libeigen3-dev` failed with
  `Unable to locate package`; Eigen 3.3.7 then failed C++20 rewritten equality
  candidates (return type is not bool). Official Eigen 3.4.0 in a temporary
  dependency directory compiled the new native regression. The unchanged
  required `make all -j2 EIGEN_DIR=/tmp/ou3-deps/eigen34` was then run.
  These failures concern build dependencies, not a mathematical premise.
- That parallel top-level invocation exposed the existing unordered
  `all: build test` prerequisites: `tests/kalman_ou_ii/run_tests.sh` ran before
  its executables existed (`kalman_ou_ii-sim` and `tuner_schedule-test`:
  `No such file or directory`, exit 127). Build subsequently completed.
  Rerun the required `make all EIGEN_DIR=/tmp/ou3-deps/eigen34` sequentially;
  do not change the unrelated top-level Makefile or relax a gate.
- The sequential aggregate and two redundant local validation invocations
  remained in `ArchitectureCleanupTests.test_no_retired_architecture_survives_repository`
  while scanning generated data. They were interrupted (exit 130); no failing
  mathematical assertion was observed. The same unchanged gate passes in a
  clean checkout. Remaining native suites and the data-dependent convention
  check pass separately. No input or gate was removed to manufacture a pass.
- An article patch failed exact-context matching before writing; applying
  the correct source context fixed the authoring operation.
- The first clean-checkout evidence run ran 449 tests with one failure and
  one existing data-dependent skip: the publication contract required the
  literal phrase `End-to-end regional practical stability remains`.
  Restored it with the truthful predicate `unproved`, immediately followed
  by the six-degree capture obstruction. The unchanged rerun passes.

## Retained facts

The exact full-word identity
`M_end' P_end^-1 M_end + D_word = P_root^-1`, low-rank loss factors, singular
nuisance elimination, corrected 4x4 LIN endpoint action, joint 21-coordinate
lower covariance after 16 s and 15-coordinate nuisance upper covariance after
17 s remain valid in their stated regular default A21 scopes. Stationary A21
detectability remains valid; nominal detectability is not physical capture.
Same-history bias mismatch, projection-sector and finite-bridge identities
are retained. Full AG upper covariance/loss, nonlinear supply, prefix retention
and whole-implementation arithmetic remain unproved conditional obligations.

## Alternatives

For a physically justified revised domain, require and qualify the cell bound
`||h_k a(t_k)-integral_cell a|| <= eta h_k`. Summing it gives sampled diversity
`g B_h,min-(2 V_max/T+eta) B_max`. A quantitative jerk limit J suffices with
`eta=J h_max/2`; J=100 and T=8 give 21.47475 > 0 as an **unadopted example**.
This alone does not prove capture or the nonlinear retained radius. Enlarging
the claimed region is a different target and also is not silently adopted.

## Next falsifiable experiment

Qualify a quantitative physical/sensor sampling condition and test it against
both the exact witness and the complete corrected shipping transport before
claiming it repairs capture. This premise narrows the current physical-history
domain and is not implied by its rules. Preserve the proof of
impossibility for the present six-degree target; do not return to a nominal
covariance calculation as if it could establish the missing physical premise.
