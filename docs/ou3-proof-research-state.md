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

Continue only PR #522 (`proof/ou3-alt-physical-master`) while it is open. The
original P2/P3/P4/P5 route and the sequence above remain independently
continuable. ALT does not promote their gates or consume them as a replacement
for its own complete-word certificate. The normative ALT execution plan is
`docs/ou3-alt-proof-plan.md`.

### Current hypothesis and exact positive result

Use a same-trajectory finite physical joint24 descriptor and common coercive
storage with bounded neutral/source supply. Strict homogeneous decay of held
bias and BIAS2 constant truth is not the target: a persistent coordinate cannot
contract in a uniformly coercive metric. Do not discard motion/bias cross terms.

The accepted real-arithmetic measurement graph now has an exact finite form,
not a linearized reset. With C=[c]x, U=4I+2C+cc', the deployed quaternion
parameters w,k and correction d give

`W=2w+k*c'd`, `W*(c_after-c)+k*U*d=0`.

The coefficient identity holds for all c,d,w,k, so it covers both quaternion
branches subject to their actual hard graphs and a nonzero Cayley denominator.
Together with exact residual secants and the same-beta radial projection it
gives `z_after=P_alpha*(z-B*q)`, `Sigma*q=Hbar*z+nu` with a 24x3 B. The exact
storage contribution retains every joint cross term and uses a 3x3 core. See
`docs/ou3-alt-finite-measurement-proof.md` for the derivation, hypotheses and
source/numerical limitations. This is finite-map algebra, not a contraction or
source-cover PASS.

H18 covariance must retain the latent BA marginal. On the held zero-cross-block
invariant, the reduced accelerometer innovation is
`H18*P18*H18'+R_acc+B0`. A kernel using bare R_acc has not implemented that
shipping relation; changing only its last gain would leave its earlier P
ancestry wrong. The finite graph requires the full 21-state operands and keeps
held e_ba in the residual. The 24-state storage is not marginalized.

### Limiter and failure classification

**C — finite-map/derivative representation gap.** The current endpoint object
is made from pointwise Jacobian and source-derivative products. Without a finite
identity or an anchored mean-value theorem, those products do not give the
physical endpoint error; F(0) is nonzero at a physical S event. Thus the earlier
Phase-1 assembly labels are not a closed finite-state master. This invalidates
using that object directly in a finite-state common-storage inequality. It
does not prove filter instability or invalidate the separate original route.
The common-metric entry point and candidate helper now reject it before any
coarse endpoint calculation or scale-based candidate is constructed.

The controlling missing quantity is the complete **finite** same-history error
map and its physical reference forcing, not an interval pivot or a near-one
rho. No source-uniform rho or failed common-metric impossibility is claimed.
No interval refinement is authorized merely by naming the endpoint family.

**Implementation defects in retained ALT prerequisites.** The BIAS2 adapter
re-widened an already-certified phi endpoint 1 to greater than 1, then rejected
its own result. It now consumes the certified endpoints directly; the
nonrelaxing limit is unchanged. Quadratic bias supply coefficients enclose the
exact real squares instead of using a rounded point product. The H18 enable
fixed-point check unnecessarily widened max after endpoint selection; max now
preserves its represented fixed point, while the seed square is outward.
These fixes do not change source assumptions or deployment arithmetic.

### Critic and alternatives

The strongest objection to the old candidate is mathematical: high precision,
a common metric or a thinner gain bound cannot turn an unanchored derivative
product into a finite error map. The alternatives are (1) exact finite
measurement/prediction descriptors, (2) a fully anchored radial mean-value
construction with its reference defect and domain proof, or (3) incremental
storage augmented by covariance/frontend memory plus a separate physical-truth
bridge. The first is selected for the measurement part because its exact
rational identity eliminates the finite-reset approximation rather than merely
shrinking a bound. The second and third remain alternatives, not simultaneous
new proof tracks. No replay, unreachable-root probe, metric fitting or
independent source boxes are accepted as universal evidence.

### Next falsifiable theorem work

Compose the finite measurement descriptors with physical prediction forcing and
actual held/active covariance/frontend lineage. Retain physical angular
integration/model defects rather than modeling truth as a second estimator.
Enforce every configured branch, one bias root/driver, and one Live primitive
origin. Then form the first common-joint24 master on that actual finite
relation. Before a numerical attempt can run, the finite-representation guard
must have genuine supplying lemmas, not edited flags.

A useful source-uniform storage bound, every-prefix chart retention, actual
H18->A21 landing, deployment roundoff and Mahony/proxy capture remain OPEN.
Individual BIAS definitions do not imply aggregate fresh-Live hard admission.
The old PE/vector issue remains separate unless ALT uses that lemma. All ALT
final gates and original P4/P5 remain false/unpromoted.
