# OU-III proof research state

## Current hypothesis

The user-authorized MARINE MOTION revision requires locally absolutely
continuous acceleration with jerk <=100 m/s^3. The old sampled ambiguity is
excluded. The same-history construction -> capture -> H18 -> release -> A21
architecture, quiet-water admission, all other limits and shipping behavior
remain. General physical capture and the complete stability theorem are open.

## Evidence

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

Validation: the clean-checkout evidence/publication suite passes 454 tests
with one existing data-dependent skip. Exact artifact reproduction, targeted
sampling/metric checks, ruff and Python compilation pass. No C++ source,
estimator behavior or quality threshold changes are included.

## Current limiter

The positive physical Gramian uses true world transport and has dimension
three. It does not establish the nominal corrected 21-state loss, eliminate
nuisance adjustment, bound the remaining AG covariance, or prove general
capture. The new physical prediction supply is finite; a finite supply does
not imply that the declared six-degree region is invariant. It must be
combined with strict full loss and quantitative nonlinear/prefix bounds.

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

Check the full transported loss after nuisance elimination for the realized
shipping dynamics on jerk-bounded histories, then prove a uniform positive
margin or identify a genuine admissible cancellation. Bound the actual
finite-error variation between the physical vector rows and those corrected
shipping rows. Do not equate the physical 3-D Gram floor with full-state loss.
