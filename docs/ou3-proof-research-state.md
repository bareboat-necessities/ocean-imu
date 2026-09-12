# OU-III proof research state

## Status

The non-ALT (BRMM) route now carries the named magnetometer theorem classes and
has a quantified startup obstruction. The next continuation must start from
latest `main`, read `docs/ou3-brmm-main-handover.md`, and open a new PR.

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

## Named magnetometer theorem classes (non-ALT)

Two named classes are now first-class on the non-ALT route, with their own
modules and no ALT dependency:

- `MAG-BMM150-DET-v1` (`ou3_brmm_magnetic_source_envelope.py`): commissioned
  installation class `20 <= ||B_W|| <= 75 uT`, horizontal `>= 15 uT`,
  `||b_HI|| <= 5 uT`, `||n_m|| <= 2 uT` per theorem sample, for the exact
  identity `m_body = R_true B_W + b_HI + n_m`. `Rmag` is never read as a
  deterministic bound.
- `MAG-CALL-SCHEDULE-v1` (`ou3_brmm_magnetic_call_schedule.py`): first post-Live
  `updateMag()` call within 0.04 s and every later gap at most 0.04 s. This is a
  deployment call-cadence assumption, kept separate from magnetic values.

Closed by them at the unchanged gates:

- the declared magnetic PE band is now a consequence, not an extra hypothesis:
  `||m_body||` lies in `[13, 82] uT`, so the declared
  `normal_live.magnetic_vector_norm_lower_uT = 10` and upper `200` follow, and
  the shipping `mag_init_min_mag_norm = 1e-3` guard clears unconditionally;
- the deterministic startup yaw capture: at a supplied tilt-frame error the mean
  perturbation is `E = ||b_HI||+||n_m||+2 Bmax sin(delta/2)`, and
  `|sin(theta_yaw)| <= E/H_min`. At the declared 0.02 rad startup direction
  error `E <= 8.5 uT`, `|sin(theta_yaw)| <= 17/30`, and
  `sin(0.61) >= 0.61 - 0.61^3/6 > 17/30` gives `theta_yaw < 0.61` rad. With the
  SO(3) triangle inequality the full attitude error is `< 0.63` rad
  `= 36.0962 deg`, strictly inside the declared 45 deg
  `initial_filter_entrance.attitude` set. No `1/sqrt(N)` reduction is used;
- finite-time `H18 -> A21` release timing: the shipping counter is incremented
  on every post-delay `updateMag()` call regardless of innovation acceptance, so
  250 counts are reached within `0.04 + 249*0.04 = 10.0 s` of Live and the
  elapsed `9.96 s` strictly exceeds the 1 s guard.

Not closed by them, and explicitly fail-closed:

- `MAG-BMM150-DET-v1` as declared does not force universal per-sample admission
  at the shipping `MagAutoTuner::max_sample_norm_ratio_from_mean = 0.35` gate.
  The sufficient condition is `2P/(Bmin-P) <= 0.35`, i.e.
  `P <= 0.35*Bmin/2.35 = 2.9787 uT`, against a declared combined
  hard-iron + residual budget of `7 uT`;
- the vector sine separation is not derivable from this class: with
  `||a_ng|| <= 8.8` the specific force can lie within 63.8 deg of vertical and
  `B_W` within `asin(15/75) = 11.5` deg of it, so `0.1` stays a declared PE
  hypothesis;
- eventual A21 under an arbitrary external `acc_bias_hold_` is not claimed, and
  the `H18 -> A21` joint24 covariance transport edge is untouched.

## Current failure analysis

### D — PE / metric-memory domain consistency (CLOSED)

`RuntimeError: declared PE does not refine vector certificate` from
`ou3_brmm_riccati_tube.py::_declared_vector_alpha6` was a padding omission, not
an enclosure defect. `ou3_vector_uco_certificate.PE` carried the *unpadded*
`REFERENCE_BODY_RATE_DEG_S = 30.0`, while every declared COMPLETE-BRMM domain
carries the padded `BODY_RATE_NORM_MAX_DEG_S = 35.0`, so the refinement test
`rate > base rate` fired on every call. `PE["body_rate_norm_upper_deg_s"]` is now
bound to `ou3_brmm_complete_physical_envelope.BODY_RATE_NORM_MAX_DEG_S`.

The quantitative permission was checked before the change: the two-packet
bracket moves from `0.98953` to `0.98778`, so `1+2/gamma^2` moves from
`1277.5970` to `1282.1123` and `alpha_6` loses 0.353 per cent. Declared
`alpha6 = 1.0605504334607669e-4 > 0`. `H18_information_lambda_min_lower` stays
`4.253919518541475e-18` and the conditional P3 composition still validates.
Raising a PE ceiling weakens the hypothesis, so the certificate is strictly
stronger, not tuned.

### C/F — padded-family startup/Mahony requalification (CLOSED, two-phase)

The single-level formulation was infeasible at the padded 8.8 m/s^2 envelope.
`INITIAL_TILT_RAD_UPPER` had been frozen at `0.955` rad from the retired
8.0 m/s^2 surface; taken from the source's own
`asin(8.8/9.80665) = 1.1137282529726653` rad `= 63.8119 deg` it left three
mutually unsatisfiable level bounds under the old metric
`P = [[1,-6.5],[-6.5,163.25]]`: seed containment needed `C >= 1.8540056`, the
87 deg chart allowed `C < 1.4912551`, emptiness factor `1.2432518`. The boundary
flow still closed with margin `0.0357940`, so it was a failure of the
*formulation*, not of the enclosure.

The replacement keeps one metric and splits the roles. Writing `u=Rz`, `y=u/||u||`
and `M=R A_s R^-1`,

`Vdot = C y'(M+M')y + 2 sqrt(C) y'Rw = -sqrt(C) [ sqrt(C) q - 2 sup ]`,

and `q`, `sup` are level-*independent*. So inward flow is a **lower** bound on
the level, `sqrt(C) >= W_in := 2 sup_max/q_min`, and it holds on *every* level
above `W_in`, not only on a certified boundary. That gives two nested levels of
the better-conditioned metric `p=3, c=11` (`P=[[1,-3],[-3,130]]`, `det=121`):

- seed level `1.6707881129` <= outer level `1.7689` < chart ceiling
  `1.8726722863`, chart headroom `5.5414%`;
- outer boundary margin `0.0262974`, inner boundary margin (`sqrt(C_in)=1`)
  `0.0134738`;
- `W = sqrt(V)` obeys `Wdot <= -(q_min/2)(W - W_in)` with
  `q_min = 0.0387576`, so the capture rate is `0.0193788 /s`
  (time constant `51.6029 s`) and `W_in = 0.8578833`.

`{V <= C_out}` is therefore forward invariant and contains the seed, so the
observer never leaves the chart, and the state enters the inner level
exponentially. Certified tilt bounds:

- all-time `84.7161 deg` (was a meaningless level radius of `86.2567 deg`);
- at the deployed 150 s startup horizon `58.2102 deg`;
- asymptotic `56.6779 deg`.

The whole dependent chain closes again: the source-order binary32/discrete
charge (`discrete_V_margin_lower = 1.518757e-4`), the frontend transition,
`ou3_startup_timeout_capture` (`timeout_aligned_branch_margin = 4.1379 deg`) and
`ou3_p4_brmm_frontend_predecessor_invariant`.

Two corrections were made while doing this. The discrete charge had used
`2*sqrt(C)*margin` as the guaranteed first-order decrease; the derivation above
carries **one** factor of `sqrt(C)`, so it had claimed twice the decrease it had
proved. With the factor corrected the margin is still comfortably positive. The
charge also hardcoded the old Cholesky row `[-6.5, 11]`; it now reads the metric
from the continuous certificate.

### C/D — startup tilt supply for the magnetic classes (ROUTE RULED OUT)

The shipping startup magnetic accumulation frame is
`tiltOnlyQuatFromBoatQuat_(attitudeReferenceQuat_())`, i.e. the private
observer's own tilt, and the handoff seed is
`boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)`, so the proxy tilt
drives both the gauge error and the seed's tilt error. The magnetic classes
therefore need a private-observer tilt *accuracy* bound, and the route only has
a level-set radius. Exactly over the rationals:

- `min_horizontal_fraction = 0.05` needs `E <= (15 - 0.05*75)/1.05 = 75/7 uT`,
  hence `sin(delta/2) <= 13/525` and `delta <= 0.04952887185726867` rad
  `= 2.8378 deg`. This is the binding requirement;
- non-vanishing north (`E < H_min`) needs `sin(delta/2) < 4/75`, i.e.
  `delta < 0.10671729940461795` rad `= 6.1145 deg`;
- the certified all-time tilt is `84.7161 deg`, a shortfall factor of `29.8528`
  against the binding gate and `13.8551` against north capture; the asymptotic
  `56.6779 deg` still falls short by `19.9725` and `9.2695`;
- the certified tilt alone already exceeds the declared 45 deg entrance set by a
  factor `1.8826`, so no yaw gauge, however accurate, repairs the entrance from
  this supply.

The declared `startup.world_averaged_gravity_direction_error_upper_rad = 0.02`
satisfies the binding requirement with margin (1.1459 deg against 2.8378 deg),
but it is the low-passed world-gravity direction error, not the private
observer's tilt error. It is recorded as an *unsupplied hypothesis*, and
`declared_startup_direction_error_is_the_accumulation_frame_error` is false.

**The route, not just the metric class, is ruled out.** The certificate
quantifies over the sector value `s` as an independent parameter in
`[sinc(87 deg),1]`, so its conclusion must also hold for the member `s=1`. For
that member the deployed PI gains give the tilt transfer function

`theta/r = -(0.01 + 0.1 j w) / (0.01 s - w^2 + 0.1 s j w)`,

and an admitted primitive-bounded sinusoid `r = j w xi` produces
`|theta| = w |theta/r| |xi|`. That response **peaks inside the declared physical
band**, at `w = 0.117 rad/s = 0.0186 Hz` against a declared support of
`[0.018, 0.88] Hz`, giving `0.14678877` rad from the primitive channel alone;
superposing an admitted DC mean chord on the linear `z`-dynamics adds `0.05`
rad. Hence no certificate in this formulation — any metric, any level, any
subdivision depth — can bound the tilt below

`0.19678877` rad `= 11.2752 deg`,

which is `3.9732` times the binding requirement and `1.8440` times north
capture. This supersedes the earlier `10.0855 deg` scaling estimate, which used
the sector value at the chart edge rather than at the small equilibrium and was
not a defensible lower bound. The limiter is the declared forcing pair
`(mean chord 0.05, primitive 1.0 s)` together with the commissioned magnetic
band — not the storage geometry. Because the obstruction is a frequency-response
peak, a static quadratic metric cannot see it; a frequency-dependent
multiplier/IQC on the primitive channel is the instrument that could.

### F — timeout branch has no yaw gauge

`ready_by_timeout` requires only `proxy_ready`, `t_ >= timeout_sec` and
`mag_gravity_aligned_branch_`; it does not require `north_ready`. The handoff
then takes `proxy_handoff_yaw_sigma_free_rad` with `pending_yaw_abs_rad_` NaN.
`timeout_branch_yaw_gauge_forced` is false. Making `north_ready` a consequence
rather than an assumption needs a forced-acquisition argument
(`mag_tilt_fallback_sec = 30 s` plus `mag_min_window_sec = 15 s` against
`timeout_sec = 150 s`), which in turn needs the un-forced `0.35` norm-ratio
admission above.

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
- joint24 compatible storage;
- `MAG-BMM150-DET-v1` and `MAG-CALL-SCHEDULE-v1` as *named* classes with their
  own non-ALT modules, separating magnetic values from magnetic call cadence;
- the two-phase private-Mahony certificate: one metric, an outer level for seed
  containment and chart retention, an inner level for the ultimate bound, and
  the level-independence of `q` and `sup` that makes both work;
- the exponential capture estimate `Wdot <= -(q_min/2)(W - W_in)` and its
  `W_in = 2 sup_max/q_min`.

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
  history;
- refine the private-Mahony metric, level or boundary subdivision toward the
  magnetic gates: the declared-formulation tilt floor `0.19678877` rad is
  `3.9732` times the binding requirement and `1.8440` times north capture, and
  that floor is metric-independent (two-strike frozen). Metric work aimed at
  *chart retention* was permitted and is now done; metric work aimed at the
  magnetic gates is not;
- freeze a startup seed angle as a constant detached from the source's own
  `instantaneous_direction_angle_upper_rad`;
- claim a first-order Lyapunov decrease of `2*sqrt(C)*margin`: the derivation
  carries one factor of `sqrt(C)`;
- read the declared `world_averaged_gravity_direction_error` as the private
  observer's accumulation tilt-frame error;
- use `Rmag`, or any covariance, as a deterministic magnetic source bound;
- claim eventual A21 under an arbitrary external `acc_bias_hold_`.

## Alternatives for the magnetic capture

The private observer cannot supply the declared gates in the current
formulation, and no static quadratic metric can. Two independent gaps must both
close, and they are now quantified:

**Gap A — a tilt certificate near the formulation floor.** The floor is
`11.2752 deg`; the best *certified* member of the quadratic class is
`56.6779 deg`, and a non-promoting scan of the whole `(p,c)` grid found nothing
below about `41 deg`. Since the obstruction is a frequency-response peak at
`0.0186 Hz`, the instrument is a frequency-dependent multiplier/IQC on the
primitive channel, not a tighter static metric. `ou3_p4_affine_hard_tube_iqc.py`
and `ou3_brmm_acceleration_moment_iqc.py` already carry that machinery.

**Gap B — a commissioned band that forces the shipping gates at that tilt.** The
gates depend on `Bmax` and `H_min`, not on tilt alone:
`E <= (H - f Bmax)/(1+f)` with `E = P + 2 Bmax sin(delta/2)`. At the
`11.2752 deg` floor the required horizontal floor is

| `||B_W||` ceiling | required horizontal floor |
| --- | --- |
| 45 uT | 18.8833 uT |
| 55 uT | 21.4462 uT |
| 65 uT | 24.0092 uT |
| 75 uT | 26.5721 uT |

Every row is feasible within its own total field, so **Gap B is viable** — a
commissioned mid-latitude class such as `||B_W|| <= 55 uT` with horizontal
`>= 21.4462 uT` would force the gates, at the cost of excluding polar and
very-strong-field installations. It must also fix the un-forced `0.35`
norm-ratio admission, which needs `P <= 7 Bmin/47`.

Two further routes remain open and do not depend on Gap A:

- **Tighten the gravity-direction forcing qualification.** The primitive channel
  contributes `w |theta/r| |xi|` at the in-band peak, so `(0.05, 1.0 s)` must
  come down together; reaching the binding `0.0495289` rad budget needs roughly
  mean chord `<= 0.0163` and primitive `<= 0.226 s`. This is a claim about
  marine gravity-direction rectification and must be argued from the
  COMPLETE-BRMM spectral support, not chosen to fit. Note the direction of the
  risk: a pathwise primitive bound derived from the band's low-frequency edge
  could be *larger* than 1.0 s, not smaller.
- **Prove the quality handoff strictly precedes the timeout handoff.** With
  `mag_tilt_fallback_sec = 30 s`, `mag_min_window_sec = 15 s` and
  `timeout_sec = 150 s`, forced acquisition would make `north_ready` a theorem
  consequence and remove the ungauged branch entirely. It depends on the
  norm-ratio admission result, not on a new tilt bound.

Do not revisit the static-metric route for the magnetic gates: the
`3.9732` formulation-floor shortfall is metric-independent.

## Next falsifiable sequence

1. Attempt Gap A with a frequency-dependent multiplier/IQC on the primitive
   channel and measure how close to the `11.2752 deg` floor it gets. Do not
   retry a static quadratic metric.
2. Decide Gap B: either narrow `MAG-BMM150-DET-v1` to a commissioned band that
   forces the gates at the achieved tilt, or tighten the forcing qualification
   from the COMPLETE-BRMM spectral support. Then close the ungauged 150 s
   timeout branch.
3. Complete continuous same-history physical source -> IMU/frontend -> tuner ->
   P/H/R/K attachment for every allowed word and BIAS family.
4. Close the consecutive compatible joint24 endpoint inequality
   `W_next <= rho W + C`, `rho < 1`, with explicit coercivity and metric-change
   bounds.
5. Close literal every-prefix augmented LDLT for every shipping subevent.
6. Prove first-exit chart/domain retention and determine the largest certifiable
   P4 basin from physics, not convenience.
7. Prove finite H18 capture and the actual H18->A21 joint24 covariance
   release/guard transport; only its release *timing* is closed
   (`<= 10.0 s` of Live under `MAG-CALL-SCHEDULE-v1`).
8. Add deployment finite-precision enclosure and compose the end-to-end theorem.

Use failure classes A/B/C/D/E/F/G exactly as established. The old DC-position
witness is outside corrected COMPLETE-BRMM and may not be recycled as A or B.

## Parallel ALT contraction/dissipativity track

### Current hypothesis

Construct the finite physical joint24 runtime word before searching coercive
storage with bounded neutral/source supply. Keep the original route above
independently continuable. `docs/ou3-alt-proof-plan.md` is normative; the ALT
handover and supplying finite notes specify the current representation boundary.

### Retained facts and representation

Finite physical prediction retains continuous attitude increments/angular defect,
correlated q15 translation moments, gyro-bias drift, one shared accelerometer-
bias driver and one Live S origin. Finite measurement descriptors retain the
nonhomogeneous physical S residual, full H18/A21 covariance, inverse-free
innovation relation, masked Joseph identity, finite injection and same-beta
projection. Prediction covariance/runtime and frontend/tuner temporal composers
remain conditional where their arithmetic, source or branch inputs are unbound.
None is the completed source-uniform 600-step master.

The private observer now has an optional **profile-specific initialized
binary32 graph**, supplied by `finite_binary32_mahony.py` and
`docs/ou3-alt-mahony-binary32.md`. Exact integer/rational nearest-even rounding,
an independent midpoint-cell checker, the literal bit seed/Newton normalization,
and the initialized Mahony/readout program remove free reciprocal choices on
that path. Every rounding defect belongs to its actual operands. The existing
vertical-result type feeds the same downstream consumers. This does not qualify
the actual target compiler, initial seed or nonfinite branches; the default
conditional real helper is not silently relabeled as a deployed float proof.

### Evidence and validation boundary

All 13 new local tests and the existing private-vertical/core/BIAS tests pass
(39 combined). Native correspondence observes the actual headers from public
startup/update calls, never installs a synthetic runtime root, and checks 767
scalar normalization boundary cases and 51 initialized observer steps. These
are implementation checks, not source admission, state-domain invariance or
storage evidence. The all-input profile theorem is the operation-by-operation
identity in the supplying note, not an inference from those cases.

The default local `make all` stops while compiling
`tests/ahrs/ahrs-qmekf-sim.cpp`: `Eigen/Dense: No such file or directory` on the
configured `/usr/include/eigen3` path. The native binding regression uses an
available Eigen include directory explicitly and passes. This is an environment
failure, not a mathematical failure, and no complete local build is claimed.
The next build check is the unchanged primary build with its proper Eigen/data
prerequisites. Focused CI and the inherited full suite remain separate; passing
one does not qualify the other or the theorem.

### Current limiter: C/E, plus unqualified deployment arithmetic

The actual complete finite source/runtime relation is still missing. Remaining
work includes all unconditional frontend/tuner/sensor-model bindings, numerical
factorization and transcendental branches, asynchronous magnetic histories,
every H18/A21 edge and literal prefix, corrected COMPLETE-BRMM/BIAS admission,
and a same-history uniform deployment-roundoff enclosure. Initial seeding and
compiler/FMA/reduction-profile qualification are explicit on the new Mahony path.
A returned finite map does not establish indefinite retention on its domain.

There is no certified rho, retained basin, ultimate bound or capture time. No
storage attempt, common-metric impossibility or new metric-failure strike is
implied by this representation work.

### Critic and alternatives

The strongest objection is model drift: a Python evaluator plus successful host
traces could miss a different deployed expression/reduction profile. Keep exact
rounding-cell checks and source correspondence, but do not promote host equality
to deployment qualification. Remaining normalization/roundoff alternatives are
an exact profile-bound program (current choice), the existing outward scalar
contract on the SAME operand graph, or a verified bit-vector/compiler extraction.
These are attachment methods, not permission for another frozen observer or a
new proof track. No route is selected merely because a tighter box looks better.

### DEAD_ENDS and next falsifiable work

Do not revive derivative cocycles as finite maps, independently boxed operands,
frozen/replayed words, more-seed source qualification, wordwise S resets, position
reanchoring, an 18-state A21 marginal, or covariance-consistency entry assumptions.
The two-strike rule and `assert_finite_storage_master` remain in force.

Complete the remaining same-predecessor runtime bindings and actual target
arithmetic profile; verify the generated finite event/prefix graph against its
supplying source without freely chosen numerical outputs. Attach one physical
and bias history through the complete word. Only after that relation satisfies
the finite-master guard run the high-precision feasibility diagnostic and common
joint24 storage search. Then close useful supply, every-prefix retention,
H18/A21 transport, fresh entry, startup capture and deployment precision.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Original P4/P5 are not promoted by ALT.
