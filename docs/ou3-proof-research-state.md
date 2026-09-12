# OU-III proof research state

## Status after PR #516

PR #516 repairs the BRMM source-definition omission exposed by PR #515 and
continues the OU-III proof without changing the shipping filter or the frozen
`P3 delta = 1e-18` requirement. The next continuation must start from latest
`main`, read `docs/ou3-brmm-main-handover.md`, and open a new PR.

The end-to-end theorem is still open. `P4_PASS=false` and
`P5_MAY_START=false` remain intentional fail-closed outputs.

## Physics-first COMPLETE-BRMM

The theorem now defines the physics before any mathematical source
construction. BRMM wave displacement `p_wave` is bounded oscillatory vessel
response about a local equilibrium, not absolute/global vessel translation.
Current, propulsion, leeway, secular drift, arbitrary global origin offsets and
other slow/global translation are outside `p_wave`.

The hard deterministic wave property is a bounded centered primitive:

`|| integral_u^t p_wave(s) ds || <= D_S`

for every admitted history and relevant `u,t`. Zero mean alone is not used as a
substitute. Harmonic/spectral and bounded shaping-state representations are
only sufficient certificate methods for the physical condition.

The old exact `p=d != 0, v=0, a=0` witness is retained as a regression for the
old finite-window source. Its historical classification is B under the old
formal source surface and E as the intended physical source-specification
omission. Under corrected COMPLETE-BRMM it is excluded by theorem, not by an ad
hoc flag. No shipping-filter instability is inferred from that witness.

## Quantified physical envelope

Using the existing Hs=8.5 m reference family with fixed outward engineering
padding, the current hard proof envelope is:

- `Hs <= 9.35 m`
- `||p_wave|| <= 8.10 m`
- `||v_wave|| <= 5.50 m/s`
- `||a_wave|| <= 8.80 m/s^2`
- `||omega_body|| <= 35 deg/s`
- `f in [0.018, 0.88] Hz`
- `D_S <= 1100 m*s`

For the 28-ft finite-harmonic reference, the deterministic primitive derivation
gives `D_S <= 863.7794 m*s`. Applying the declared amplitude and lower-frequency
padding gives `<1056 m*s`, rounded outward to 1100. The old 300 m*s P4 value is
not a physical source cap and must not be used to define admissibility.

The one-time Live coordinate remains `S_L(t)=S(t)-S(t_L)`. There is no wordwise
S reset and no position reanchor.

## Correlated innovation result

The proof now preserves covariance provenance through one correlated
measurement object:

`(P,H,R) -> PHt -> S=HPH^T+R -> S^-1 -> K -> Joseph`.

This removes the earlier artificial innovation singularity caused by replacing
`S` with an independent entrywise rectangle. Structural `P >= 0` and `R > 0`
are retained in the inversion argument; impossible singular members of the
rectangular hull are no longer interpreted as physical source histories.

This repair is proof machinery only. It is not a filter change.

## H18 eta6/a_w information certificate

After removal of the false singularity, the first quantitative H18 lower was

`7.092471820569811e-19`,

which was strictly positive but only 70.9247% of the frozen 1e-18 gate. The
limiter was the coupled `eta6/a_w` block; other translation directions had large
headroom.

Two same-history proof tightenings are now canonical:

1. Four guaranteed actual S=0 firings are selected from scheduler windows
   `[0,g]`, `[4g,5g]`, `[8g,9g]`, `[12g,13g]`. The last ends by
   `1.9499999564 s` inside the same 3 s word. All other due S updates stay in
   the literal word. OU response and process-noise charges are recomputed over
   the longer selected horizon.
2. The second PE accelerometer occurrence retains the homogeneous OU attenuation
   of the same initial `a_w` coordinate. It is not conservatively re-created
   with unit sensitivity as though it belonged to an independent history.

The canonical post-Live correlated-numerics CI now reports:

- `H18_information_lambda_min_lower = 4.253919518541475e-18`
- frozen gate = `1e-18`
- gate ratio = `4.253919518541474`
- `coupled_eta6_aw_lambda_min_lower = 4.253919518541476e-18`
- `accelerometer_translation_cross_norm_squared_upper = 3767421.6507477993`

Directional translation lower bounds are:

- `S = 1.8434197928164813e-11`
- `g*p = 2.3086404646316095e-11`
- `g^2*v = 4.447787288026003e-10`
- `g^3*a_w = 3.615365438779952e-08`

The former eta6/a_w blocker is therefore closed at the unchanged gate. This is
about a 5.9978x improvement over the preceding H18 lower and does not promote
P4/P5 by itself.

## Current failure analysis

### C/D — PE / metric-memory domain consistency

A current metric-memory path fails before its intended diagnostic at

`RuntimeError: declared PE does not refine vector certificate`

from `ou3_brmm_riccati_tube.py::_declared_vector_alpha6`.

This is presently a proof-domain/representation consistency issue. It is
separate from the now-gate-passing H18 eta6/a_w information certificate. Do not
lower the 1e-18 gate or retune the filter to hide it. Inspect which declared PE
object and vector certificate are being compared, preserve same-history source
ancestry, and determine whether this is C (dependency/enclosure) or D
(entry/working-domain modeling).

### F — padded-family startup/Mahony requalification

The physical acceleration cap was widened from the earlier 8.0 m/s^2 proof
surface to 8.8 m/s^2. The real private Mahony/proxy runtime handoff still needs
a rigorous invariant/capture requalification over that padded family. Preserve
both measured-period takeover and prior-frequency timeout Live entry. Do not
make measured-period availability a hidden theorem prerequisite.

### C/G — downstream same-history joint24 closure

The controlling P4 theorem remains the complete-word, same-history joint24
error/true-bias inequality with compatible consecutive storage. Required
remaining pieces include actual source forcing, H18/A21 prior-free and
finite-bias transport, source-uniform endpoint contraction, every-prefix
augmented LDLT, hybrid metric compatibility, and finite precision.

No source-uniform worst endpoint rho, every-prefix gain, maximum retained basin,
finite H18 capture time or H18->A21 basin landing is yet certified.

## Retained proof architecture

Preserve:

- branch-complete production lineage;
- source-reachable selector family;
- H18 source-indexed prefix transport;
- A21 joint24 source-indexed prefix transport;
- BIAS0/BIAS1/BIAS2 separate recurrences and same true-bias ancestry;
- actual-applied anisotropic `R_S`;
- same-signal `f -> sigma -> tau -> T_S -> R_S` relation;
- candidate/active tuner and staged commit/scheduler state;
- correlated same-history P/H/R/K;
- Joseph/reset/projection splitting;
- one-time Live S-origin handling;
- finite-precision obligations;
- joint24 compatible storage.

Do not resurrect the A21 18-state marginal-motion storage. Prior diagnostics
showed that discarding motion/bias cross-information can make the marginal
storage catastrophically noncontractive while the full joint storage remains
useful.

## Dead ends / forbidden shortcuts

Do not:

- choose `D_S` because a P4 search likes it;
- assert `S <= 300 m*s` as a physical premise;
- use zero asymptotic mean as an all-time bounded-primitive proof;
- use replay, captured trajectories or finite seed sets as universal source
  qualification;
- multiply independent coefficient rectangles that destroy same-history
  dependence;
- wordwise re-zero S or re-anchor p;
- use Gaussian/high-probability source events in a deterministic theorem;
- use covariance consistency as hard entry membership;
- lower P3 delta below `1e-18`;
- change the shipping filter merely to stabilize the excluded DC-position
  history.

## Next falsifiable sequence

1. Fix/classify `declared PE does not refine vector certificate` while retaining
   the now-canonical H18 information construction.
2. Requalify Mahony/proxy startup for the 8.8 m/s^2 padded physical family.
3. Complete continuous same-history physical source -> IMU/frontend -> tuner ->
   P/H/R/K attachment for every allowed word and BIAS family.
4. Close the consecutive compatible joint24 endpoint inequality
   `W_next <= rho W + C`, `rho < 1`, with explicit coercivity and metric-change
   bounds.
5. Close literal every-prefix augmented LDLT for every shipping subevent.
6. Prove first-exit chart/domain retention and determine the largest certifiable
   P4 basin from physics, not convenience.
7. Prove finite H18 capture, actual H18->A21 release/guard transport, and
   indefinite continuation.
8. Add deployment finite-precision enclosure and compose the end-to-end theorem.

Use failure classes A/B/C/D/E/F/G exactly as established. The old DC-position
witness is outside corrected COMPLETE-BRMM and may not be recycled as A or B.

## Parallel ALT contraction/dissipativity track

### Current hypothesis and retained facts

Construct the actual finite source/runtime joint24 word before searching
coercive storage with bounded source supply. Preserve the original route above;
`ou3-alt-proof-plan.md` is normative. Full-21 covariance, exact physical
prediction, inverse-free innovation/Joseph/reset, same-beta projection,
startup-rooted H18 and persistent frontend/calibration/control state remain.
The same-operand finite-real tilt reset is attached; universal deployment
arithmetic correspondence is not.

`ou3-alt-source-continuation.md` proves the moment concatenation identity
`G(a+b)=T(b)G(a)T(b)' + G(b)` and its nonnegative projection loss. Therefore
`E(J_prefix,H) <= sum_i E(J_i,h_i) <= A_max^2 H` on every finite prefix.
The actual physical segments drive this derived budget and the existing finite
predictor. Endpoint vector caps, the coupled nine-moment IQC, a necessary
rotation chord/rate bound, exact BIAS envelope endpoints and one actual selected
bias factor are checked before an IMU event. These are necessary outer
constraints, not complete BRMM/BIAS generating-history membership.

The fresh Live origin now has a separate checked outer endpoint: sample-zero
magnetic calls consume the actual fresh `Reference` without manufacturing a
predecessor transition, and leave source ordinal 1 available to the first IMU.
This removes the sample-zero topology gap but does not establish generator/QO or
full BIAS-history membership.

The magnetic graph retains default continuous calibration, refinement and
coupled offset/reference writes, including rejected-branch bookkeeping.
Real-arithmetic calibration <=28.7 uT and active reference/corrected observation
<=110.7 uT are boundedness results, not calibration accuracy or contraction.
The corrected unlock theorem uses continued locally finite calls and
`first_gap + max((n-1)*gap,1+gap)`; default internal unlock remains within 10 s.
An upper gap bound never forces the n-th call to satisfy the strict >1 s guard,
and arbitrary external hold does not imply eventual A21.

### Failure analysis and independent critic

The source-product fixture referenced nonexistent `PhysicalSegment.delta_theta`
/ `delta_velocity` / `delta_position`, not the actual coupled `J0/J1/J2`.
This test implementation defect is repaired without changing runtime physics.

The earlier source constructor checked labels, clocks and bias envelopes but
not physical vector/moment constraints. Matching generator/history strings
cannot certify O^601_BRMM membership. An algebraically consistent 100 m/s^2
endpoint, or an independent-sign moment corner with energy 193*h*A_max^2,
passed the missing checks. These are source-binding regression inputs, NOT
admitted shipping instability counterexamples. The checks now reject them.

The critic's strongest objection is model drift through qualification metadata.
Alternatives are (1) an executable necessary correlated outer graph with honest
open membership, (2) a verified admitted-generator evaluator, and (3) symbolic
source-function substitution into shipping code. Use (1) to attach the already
declared physical conditions; (2)/(3) remain the full-source binding obligation.
This constructs ports in the controlling storage inequality, not a refinement
of an unproved rho margin. No storage feasibility or contraction strike follows.

The focused guard integration exposed an unchanged legacy serialization defect:
`ou3_p4_brmm_physical_acceleration_witness_sector.build` calls
`float(mom['A_max_mps2'])` on the exact string `'44/5'`. This is an implementation
failure before any storage inequality, not a negative contraction margin.
The pure finite-storage barrier is independently tested without pretending a
source builder passed. Full `phase1_closure.build/validate` still execute every
source prerequisite; the inherited integration test remains, with its obsolete
storage-true expectation corrected to false. Shared P2/P3/P4/P5 code is untouched.

### Evidence and validation boundary

The full focused finite-identity selection passes 511 tests with native checks
required. The native shipping finite-identity check passes with bit-identical
observed/plain sample states. Native preserve-yaw reset regression passes,
including its intermediate covariance axis and near-parallel branch. These are
implementation/algebra regressions, not all-history deployment qualification.
The full inherited ALT integration suite is not claimed passing; the shared
`float('44/5')` failure remains exposed by its source-building path.

Default `make all` stops compiling `tests/ahrs/ahrs-qmekf-sim.cpp` with
`Eigen/Dense: No such file or directory` on `/usr/include/eigen3`. Native checks
use an existing separate Eigen include. No full-build PASS is claimed.

### Current limiter, DEAD_ENDS and next falsifiable work

The complete finite source-uniform word remains open: bounded generator/
potential and Q/O continuation, full BIAS driver/parameter histories, remaining
same-source exp/trig/Eigen/roundoff branches, full sample-zero generator membership (the checked outer endpoint topology is
closed), startup/ungauged capture and continuation beyond the current
600-step container.
The firing tilt/reset operands are no longer free; their deployment boundary
and nonfinite behavior still need proof. Floating clocks and the signed magnetic
counter cannot be replaced indefinitely by Python rationals/integers.

Do not revive derivative cocycles as finite maps, frozen gains/replays,
independent coefficient boxes, more-seed source qualification, wordwise S resets,
18-state A21 marginal storage, convenient entry sets or covariance consistency.
The two-strike rule and `assert_finite_storage_master` remain in force.

Next bind a verified source realization or symbolic same-history source graph
to the checked moments, BIAS driver and every runtime operand; demonstrate
universal inclusion and all branch successors rather than more matching labels.
Only a complete finite master may authorize common joint24 storage search,
followed by every-prefix retention, capture and an ultimate bound.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
