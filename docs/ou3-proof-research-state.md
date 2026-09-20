# OU-III proof research state

## Current hypothesis

The user-authorized MARINE MOTION revision requires locally absolutely
continuous acceleration with jerk <=100 m/s^3. The old sampled ambiguity is
excluded. The same-history construction -> capture -> H18 -> release -> A21
architecture, quiet-water admission, all other limits and shipping behavior
remain. General physical capture and the complete stability theorem are open.

## Evidence

`ou3-corrected-word-proof.md` proves the following complete implications:

- the existing embedded nuisance floor passes through actual acc/S corrections
  and arbitrary attitude resets to roots immediately before prediction;
- with that floor and the existing nuisance upper bound, a uniform SPD bound J
  on the six AG columns of the actual complete loss gives a full covariance
  upper bound by a Schur complement; the first prediction then gives strict
  21-state contraction, retaining all cross covariance;
- actual-gain variation of constants composes finite-error supplies without
  equating independently evolving filters; correction supply is bounded by
  the R-whitened residual without a separate gain norm; the exact covariance
  recursion bounds all prediction/measurement inputs jointly by the square
  root of their summed actions, retaining process correlations;
- an explicit finite-angle reset remainder includes both the injection-squared
  times error term and the real normalized quaternion polynomial defect.

The six-column uniform J premise remains OPEN. The exact nuisance ratio is
about 1.02976e-18 and the scalar full process floor about 2.49159e-26. These
are valid existence ingredients, not a useful practical-radius margin; do not
substitute this coarse prediction-only bound for a well-conditioned complete
loss calculation. The 80-digit cross-coupled diagnostic gives rho about
0.995037970315 and a conditional certified decrement about 2.64851e-5;
its independent rational verification is a supplied-word check, not shipping
history evidence. The new certificate remains fail-closed.

`ou3-sampling-fidelity.md` proves:

- sharp nonuniform trapezoidal acceleration mean error <=J h_max/4=.15;
- exclusion of every constant-attitude stationary-accelerometer alias of
  >=6 degrees over 32 s, including full accelerometer residuals;
- a joint 3-D measured-vector Gram floor >6.29714e-5 over 64 s for constant
  world field, charging both 5-uT hard-iron and 2-uT measurement residuals;
- finite-angle rotation coercivity of that same full 3-D Gram matrix;
- physical LIN mismatch action and prediction supply in the full carried
  covariance metric, with source transition defects and process correlations.

`sampling-fidelity.json` reproduces exact rational constants. The all-time
continuation schema fails closed without both a jerk bound and acceleration
absolute continuity. Finite traces still do not certify all-time membership.
Reference-generator coefficient screening gives maximum harmonic jerk envelope
about 29.1126 m/s^3 across its 20 scenarios; it supports the selected 100 limit
but is not measured-vessel qualification or a numerical step in the proof.

The prior counterexample remains exact on the domain without the new bound:
its jerk exceeds 10500, its physical tilt is 53.13 degrees, and the source
all-time actual magnetic-service floor exceeds one. Its native regression
checks literal finite float startup/release; all-time float totality was not
claimed. Current status explicitly records that this witness is inadmissible.

Validation for this continuation: 463 evidence/publication tests pass with
one existing data-dependent skip. Exact artifact reproduction, ruff and
Python compilation pass. The normalized quaternion polynomial defect was
also checked at 80 decimal digits. No C++ or shipping behavior changed.

## Current limiter

The six AG columns must be bounded for the actual nominal corrected transport,
including realized gains and resets. The physical three-dimensional Gramian
and two-dimensional magnetic-service restriction do not provide that bound.
Once a useful full loss is certified, the actual-gain supplies and projection
sector must close a retained region and startup/bridge/release must enter it.
A finite source supply and a positive but extremely small comparison margin
are not evidence that the declared six-degree region is invariant.

## Failed approaches / DEAD_ENDS

- Entrywise midpoint-radius Riccati refinement failed twice by dependency:
  prediction [-18.7907040,2502.41442] and failed acc/mag inverse boxes.
- Cross-ceiling/Gershgorin reduction gave gamma=-5.1223e8. An LDL pivot is
  not an eigenvalue floor; keep matrix action and factors.
- Accumulated process noise plus endpoint information ceilings fails with
  interleaved corrections. The corrected endpoint path action is retained.
- Restricted heading/bias loss I_2 can coexist with full loss
  `[[1,0,1],[0,1,0],[1,0,1]]`; nuisance cancellation (1,0,-1) defeats lifting.
- Finite replay, QR positivity alone, and a literature map did not discharge
  uniform stability. They are not substitutes for the full loss implication.
- Continuous velocity balance could not be transferred to point samples
  without sampling fidelity. The old all-time witness refuted six-degree
  capture. The newly adopted jerk assumption now closes that specific gap.
- Service induction must initialize the first partial magnetic cell, not
  assume handoff/clock coincidence. Its exact seed/prefix comparison is fixed.
- Previous local infrastructure failures: missing Eigen; Eigen 3.3.7 C++20
  equality incompatibility; unordered `make all -j2` build/test prerequisites;
  generated-data architecture scan interrupted (exit 130). Eigen 3.4 and
  clean-checkout validation resolve the dependencies/scanning issue. The last
  sequential `make all` compiled all targets and passed through spike_filter
  before interruption; remaining native suites passed separately. Do not
  report that interrupted aggregate as passed.
- Local TeX lacks IEEEtran/luaotfload; use the unchanged CI renderer. HTTPS
  push lacks credentials; use identical Git trees through the connector,
  preserving executable bits and fetching the explicit branch ref.

- The first broad test invocation in the generated-data checkout again
  stalled after its first test in the architecture scan and was interrupted
  (exit 130). Run the unchanged evidence suite in the clean detached validation
  worktree; the targeted new tests and exact artifact reproduction pass.

- CI article rendering at 5ba1e038 failed with `Environment proposition
  undefined` at line 420. Reuse the declared lemma environment and split the
  longer equations to fit the column. Mathematical certificates and literal
  shipping CI passed; this was a publication-source failure. The provenance
  gate then correctly rejected the changed article hash; refresh its binding
  and rerun the unchanged gate.

- Visual review of the rendered bde2d251 article caught an unqualified
  old-domain capture-refutation sentence in the evidence section and a missing
  jerk row in the constants table. Restrict that sentence to omission of the
  jerk premise and include J_max=100 in the table; proof constants are unchanged.

- The old smooth-remainder helper inferred eta(r)->0 from smooth operations,
  a projection-sector flag and frozen schedules. That implication is invalid
  for the shipping reset comparison G=I+[d]/2 at a nonzero injection: its
  derivative is J_l(d), and J_l(d)-G generally is nonzero. Removed that helper,
  its dependent Boolean small-gain promotion and the accepting test. The
  explicit injection-dependent finite-error bound replaces them. Full theorem
  status was already false; the controlling next check is the actual word
  margin after this residual is charged.
- The new high-precision diagnostic initially failed because mpmath was absent
  (`ModuleNotFoundError`). An isolated environment with mpmath resolved it;
  no repository dependency or quality threshold was changed.

- The first broad evidence run in the new isolated environment failed with
  missing matplotlib imports, which also prevented publication-test adapters
  from loading. Installing matplotlib restored the unchanged suite: all 462
  tests pass with one existing skip. No publication files or gates were altered
  to accommodate that environment failure.

- The first article render for the corrected-word result passed but reported
  a 9.1585-pt overfull inline nuisance-block tuple. Expand the nonzero blocks
  into two displayed rows; the covariance constants are unchanged.

- Visual review of the joint-input-action article found a nearly empty ninth
  page containing only the disclosure. Remove the conclusion that repeats
  the abstract and evidence boundary; retain all mathematical arguments and
  the disclosure. This is a publication-layout correction only.

## Retained facts

The exact full-word identity
`M_end' P_end^-1 M_end + D_word = P_root^-1`, low-rank loss factors, singular
nuisance elimination, corrected 4x4 LIN endpoint action, joint 21-coordinate
lower covariance after 16 s and nuisance upper covariance after 17 s remain
valid in their stated regular default A21 scopes. Same-history physical bias
mismatch, projection-sector and finite-bridge identities remain. Stationary
A21 detectability does not imply capture of the nonlinear physical observer.

## Alternatives

Use the sharp sampled mean and finite-angle vector identity with actual
transport, preserving nuisance columns and the realized schedule. The exact
physical OU action retains the full process matrix; its scalar box bound is
only a finite supply ceiling, not evidence of a small retained radius.

## Next falsifiable experiment

Certify or falsify a uniform SPD bound on the six AG columns of the complete
corrected loss under the admitted varying nominal coefficients. The coupled
covariance reduction now handles the other fifteen coordinates analytically.
For a useful nonlinear margin, retain full process/measurement factors rather
than the coarse scalar prediction-only floor. Charge the actual injection
term before asserting a vanishing nonlinear gain; then check the invariant
root-radius and every-prefix inequalities on the same history.
