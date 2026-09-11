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

The independently continuable ALT route targets a finite physical joint24 word
and coercive dissipativity storage with bounded neutral/source supply. The
original P2/P3/P4/P5 route above is unchanged. The normative plan is
`docs/ou3-alt-proof-plan.md`; continuation status is in
`docs/ou3-alt-contraction-handover.md`.

### Current hypothesis and proved algebra

Use exact finite descriptors, not a product of local derivatives. Accepted
measurement/reset/projection identities retain physical `r_S=e_S-S_phys`,
held e_ba, full 21-state covariance and H18's latent BA innovation contribution.
The physical predictor now uses the actual continuous physical rotation
increment, retains angular/model forcing, all correlated q15 translation
moments, gyro-bias drift, and the same physical accelerometer-bias driver in
both e_ba and beta. Quaternion/translation identities are proved by exact
all-coefficient polynomial checks. See `docs/ou3-alt-finite-measurement-proof.md`.

A branch-correct rank-three Joseph simplification follows from the ACTUAL
inverse-free `K Sigma=N` relation, without assuming `N=P H'`. The implementation
checks that equality before cancelling; a nonzero solve defect is not dropped.
This is exact arithmetic on the full covariance, not a state/storage reduction.

### Evidence and its limits

A passive observer runs the real wrapper from startup without setting synthetic
state/covariance/gain/mode values. Three 600-step windows exercise H18, A21 and
an actual H18->A21 release. Core finite identities, consecutive state/covariance
ancestry, nonzero physical S residuals, H18 latent BA, both observed quaternion
branches, floor, due/not-due scheduler and same-beta projection are checked.
The first H18 window retains fresh centered e_S=0. Observation statements erase
to the original source, and the two host builds must have bit-identical recorded
sample states. These are implementation regressions, not source admission,
complete branch coverage, contraction evidence or deployment qualification.

### C/E — controlling limiter

The missing object is still the complete FINITE same-history source-uniform
runtime graph. Exact local algebra and recorded operand ancestry do not supply
all coefficient/product/guard equations over the entire analytic family. The
frontend/tuner/covariance successors, physical angular-defect source relation,
repair/rejection/watchdog branches and all asynchronous/hybrid continuations
must be bound, rather than substituted by snapshots, tokens or independent boxes.

This is a finite representation/source-attachment gap (C/E), NOT instability.
No common-storage formulation has been tested on a complete finite master in
this continuation. There is no worst rho, failed common-M impossibility, metric
strike, certified basin or capture time to report. Interval refinement and
storage search remain blocked.

### Critic, retained facts and alternatives

The strongest objection is that another successful trace only checks another
execution; it cannot close the missing universal quantifier. Increasing seed
count, precision, interval subdivision or metric complexity would not repair
that gap. Retain the finite local identities and passive correspondence checks
as tests of the actual future runtime graph, never as its source cover.

The alternatives remain (1) direct finite descriptors with explicit runtime
successors, (2) a complete anchored mean-value construction retaining F(0), or
(3) incremental storage with runtime memory and a separately proved physical-truth
bridge. Continue (1); do not start another route merely because the source
attachment requires work. The two-strike/architecture-review rule remains in force.

### Validation boundary

The broad inherited ALT suite still reaches a shared startup prerequisite and
can fail with `continuous Mahony invariant invalid: ['continuous_all_live_PI_invariant_closed is not true', 'initial seed angle not closed']`.
The failure is reproduced by the unchanged `test_coarse_endpoint_outer_attempt`
path. Do not weaken that shared original-route premise to make ALT green. The
finite-identity job is independent and must not interpret its own success as
success of the broad suite or source qualification.

### Next falsifiable theorem work

Attach the finite predictor and measurement graphs to the actual endogenous
covariance/frontend/tuner successors and every configured guard, with one
physical segment chain, one Live origin and one BIAS0/1/2 history per word.
Recompute geometry and operands from each finite predecessor; do not reuse the
sampled-shadow predictor or derivative lineage as a physical endpoint.
After the full finite master passes the existing representation guard, and only
then, run its high-precision feasibility diagnostic and first common joint24
storage search with justified supply.

Fresh-Live aggregate admission, useful storage/ultimate bound, every-prefix
retention, actual H18/A21 hybrid closure, startup capture and finite precision
remain OPEN. All ALT final gates remain false; original P4/P5 are unpromoted.
