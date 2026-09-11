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

The independent track is specified in `docs/ou3-alt-contraction.md`. The
existing architecture, its next sequence above, and all P2/P3/P4/P5 gates remain
unchanged. ALT does not promote those gates. The current ALT research question
is whether an inverse-free whole-word joint24 dissipation formulation is a
viable alternative, not whether isolated algebra tests establish stability.

### Current hypothesis and exact architecture failures

Keep a coercive joint24 storage, but charge non-decaying held-error/true-bias
coordinates through bounded, same-history supply ports. Seek strict decay of
the motion performance with finite ultimate bounds, not contraction of the
physical source itself.

The stricter common-metric target fails structurally: a held coordinate or
BIAS2 constant truth supplies an identity block. If `A v = v != 0`, then
`v^T(A^T M A - rho M)v = (1-rho)v^T M v > 0` for every SPD common M and
`rho < 1`. This invalidates strict homogeneous contraction of the entire
augmented state; it does not invalidate motion ISS or the existing proof.
Limiter: an exact neutral direction, not interval width. The independent
critic's strongest objection is that merely renaming P4 as contraction cannot
remove that direction. Alternatives are (1) bounded-bias supply with full
joint storage, (2) a quotient/incremental metric at fixed physical source,
and (3) finite-duration held-mode storage followed by active-mode performance
storage. The first is selected without dropping bias cross-information.

A second shortcut fails at the shipping H18 gain mask: the actual `PCt` is not
unconditionally the full residual Jacobian's `P H^T`. The unqualified full-state
information-addition identity is therefore not an allowed rewrite. Retain
`S q = r`, `delta_e = -N q`, with actual `N=PCt`, and the implemented Joseph
polynomial. An effective H18 information model requires its own invariant
bridge. This is a proof-method mismatch, not a filter defect.

### Shared-contract execution failure

`ou3_p4_bias_family_joint_iss_supply.build()` currently stops in
`ou3_p4_hard_entry_set.build()` with:

`qualified fresh-entry source failed: ['all_hard_entry_coordinates_inside_full_scale not true', 'source_uniform_reachable_other_entry_coordinates_closed not true', 'fresh_live_entry_cover_closed not true', 'one or more coordinate memberships failed']`.

Classification: entry/working-domain qualification failure. This invalidates
using the aggregate supply builder as already-qualified ALT entry evidence.
It does not invalidate the separately defined BIAS0/BIAS1/BIAS2 recurrences.
Limiter: fresh-Live hard membership, not covariance consistency. Critic:
consuming the aggregate as a successful theorem would import an unproved
startup assumption. Next falsifiable experiment: invoke each family definition
and validator separately, preserve their exact shared-error/truth recurrence
in the new lift, and report aggregate entry admission as false. No change to
the shared builder or admission gate is authorized by this experiment.

Local `make all` is infrastructure-blocked at `tests/ahrs/ahrs-qmekf-sim.cpp`
because `Eigen/Dense` is absent from `/usr/include/eigen3`. This is not a
mathematical failure. The independent CI installs Eigen before the unchanged
shipping observer is built.

### Current ALT evidence and next experiment

The separate BIAS0/BIAS1/BIAS2 definition validators all pass; this does not
close aggregate hard-entry or hardware admission. Fifteen ALT algebra tests
pass, including exact finite gain/innovation increments and a rational
masked-update analogue with joint cross storage and bounded neutral supply.
The analogue is not a shipping certificate.

The captured observer contains 29 H18 and 359 A21 retained three-second words.
Re-evaluation of its stored products at 80/120 digits gives worst ratios
0.9995136399194538514 and 0.9959857178394911698. The existing H18 ratio fixes
held-bias error at zero, while A21 retains all 21 estimated-error coordinates.
No metric is fitted to this replay. The limiting observed frozen-map margin
is H18's approximately 0.00048636.

A new attachment limitation is explicit: the observer reconstructs H=P^-1 N
and accumulates a binary32 homogeneous map, not the complete endogenous
nonlinear joint24 finite-increment map. Treating it as the latter would be a
proof/implementation correspondence error. Its strict ratios do not justify
outward certification of the unsupplied ALT master. The strongest critic is
that high precision cannot restore omitted physics or gain dependence.
Alternatives: instrument actual finite increments; derive an analytic graph
of actual N/S/residual/tuner operations; or prove a fixed-trajectory nonlinear
storage identity that legitimately avoids incremental gain derivatives.
Selected next experiment: bind the already implemented exact finite solve
increment graph to the analytic same-history shipping source, then run the
full-word high-precision feasibility experiment before any interval refinement.
All ALT and existing final theorem gates remain false/unpromoted.

### ALT continuation experiment — source attachment and rank-three product ports

The next falsifiable attachment experiment was executed without changing the
shipping filter or the original P2/P3/P4/P5 route. The event-local
estimator-owned JOINT/kernel relation executes for both H18 and A21: H/A cells
bind, literal event orders agree, and the authoritative next frontend is the
JOINT image. However, invoking the inherited theorem-facing source-selector
builder pulls in the separate Mahony startup invariant and currently stops at
`continuous_all_live_PI_invariant_closed` / `initial seed angle not closed`.
Classification: proof-architecture coupling, not an ALT Live theorem failure.
It invalidates reuse of that theorem-facing wrapper as the ALT Live attachment;
it does not invalidate the inverse-free finite-increment algebra, the joint24
bounded-supply target, or the independently continuable original proof.

The attachment gap is now sharper. The available relation is
single-execution/differential: it does not supply paired admissible estimator,
covariance, frontend/tuner and guard states with same-history `r0/r1`, `N0/N1`,
`S0/S1`, and `q0/q1`. Therefore `deltaN` and `deltaS` are not yet tied to one
paired COMPLETE-BRMM continuation and
`ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED=false`.

Rank-three structure is retained only as an exact representation optimization.
The new thin kernels use `Psi+ = Psi-K(H Psi)`, the actual-numerator Joseph
polynomial, and the rank-at-most-six storage correction without discarding any
state or cross term. The endogenous increment master can carry
`u_S=deltaS*q0` and `u_N=deltaN*q0`; for joint24 this reduces the local lifted
coordinate count from 87 to 33. This is exact only when hard same-history
product graphs enforce both identities. The ports may not be treated as
independent bounded disturbances.

Critic: the strongest reason to abandon the inherited wrapper is that it makes
a Live-word feasibility experiment depend on an unrelated, still-open startup
theorem. The selected alternative is a Live-only paired JOINT+kernel transition,
not more refinement of that wrapper. Current limiter and next falsifiable
experiment: take two estimator states over one admitted COMPLETE-BRMM source
continuation, emit paired literal H18/A21 event cells, bind both inverse-free
solves, and enforce hard `u_S/u_N` same-history product graphs. Startup remains
a separate later theorem. All ALT and original P4/P5 gates remain fail-closed.
