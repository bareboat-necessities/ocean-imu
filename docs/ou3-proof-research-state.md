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
- the startup yaw *algebra*: at a supplied accumulation-frame excursion the mean
  perturbation is `E = ||b_HI||+||n_m||+2 Bmax sin(delta/2)`, and
  `|sin(theta_yaw)| <= E/H_min`. At `delta <= 0.02` rad, `E <= 8.5 uT`,
  `|sin(theta_yaw)| <= 17/30`, and `sin(0.61) >= 0.61 - 0.61^3/6 > 17/30` gives
  `theta_yaw < 0.61` rad; the SO(3) triangle inequality gives a full attitude
  error `< 0.63` rad `= 36.0962 deg`, inside the declared 45 deg
  `initial_filter_entrance.attitude` set. No `1/sqrt(N)` reduction is used. The
  algebra is closed; its *supply* is not -- see the heading finding below;
- the `H18 -> A21` release time **after north lock**: `<= 10.0 s`, derived by a
  case split on whether the 250-count or the strict 1 s guard binds last.

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
  the `H18 -> A21` joint24 covariance transport edge is untouched;
- the release is **not** reachable from Live alone. The only call site of the
  counter-owning inner method is
  `if (mag_ref_set_ && stage_ == Stage::Live) { impl_.updateMag(...) }`, so a
  host call schedule does not advance `mag_updates_applied_` until north lock.
  On the admitted ungauged timeout path `mag_ref_set_` stays false and the
  release is unreachable however fast the host calls. North lock is exactly
  what the tilt floor below denies, so the yaw gauge and the A21 release share
  one unmet prerequisite;
- the accumulation mean is **not** bounded by a tilt error alone.
  `tiltOnlyQuatFromBoatQuat_` strips the *estimator's* yaw, not the vessel's, so
  `q_tilt * m_body` lives in a level frame that turns with the boat -- which is
  why the gauge is north *relative to the boat*. A heading change during the
  window rotates accepted samples in the horizontal plane and smears the mean,
  and a large enough excursion cancels its horizontal component. The parameter
  is therefore the total excursion `delta = delta_tilt + delta_heading`, and
  neither MAG-BMM150-DET-v1 nor any declared operating-domain quantity bounds a
  startup heading excursion. `DECLARED_SUPPLY_ENTRANCE_CLOSED` is false.

### Corrections made to earlier claims in this route

Three claims recorded earlier were unsound and have been retracted in place:

1. the yaw/entrance certificate charged only tilt error into `E`, leaving a
   legal heading excursion uncharged;
2. the release certificate inferred the strict `> 1 s` guard from
   `(n-1)*gap_max`, an *upper* bound on elapsed time -- 250 calls 1 ms apart
   clear the count at `0.249 s` with the shipping predicate still false;
3. the release certificate read parity from the *inner* method only and missed
   the outer north-lock gate on its single call site.

The discrete Mahony charge additionally claimed `2*sqrt(C)*margin` of
first-order decrease where the derivation gives one factor of `sqrt(C)`.

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
- charge a heading excursion as tilt error in the magnetic accumulation, or read
  `world_averaged_gravity_direction_error` as a total accumulation-frame
  excursion: it carries no heading content;
- infer a strict elapsed-time guard from an upper bound on call spacing, or
  claim the accelerometer-bias release from Live without north lock: the
  counter-owning call sits behind `mag_ref_set_`;
- certify shipping parity by loose token presence: compare the complete
  assignment or condition so a branch cannot change semantics while the
  certificate stays green;
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
3. Supply the total accumulation-frame excursion: declare and justify a startup
   heading-excursion bound over the window, or retain the heading history in the
   accumulation so a legal turn cannot smear the mean. Without it the yaw
   algebra has no input and the A21 release has no north lock.
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

### Current hypothesis and retained facts

Construct the complete finite source/runtime joint24 map before any storage
search. The original P2/P3/P4/P5 route above remains independent and unchanged.
`ou3-alt-proof-plan.md` is normative. All ALT theorem gates and
`storage_search_allowed` remain false.

The admitted-history product retains corrected COMPLETE-BRMM, one BIAS0/1/2
history, bounded symbolic IMU ISS forcing, full-21 covariance, the one-time
Live/S origin, frontend/tuner state, continuous magnetic calibration and
asynchronous control. Primary physical-history restriction, not token identity
or a required spectral generator, supplies the finite endpoint/moment constraints.
The exact moment concatenation and source-owned tilt/reset relations remain.

`finite_startup_joined_machine_history.py` now roots the machine guard, private
observer, default 6 Hz tracker LPF, stillness, band/statistics and TuneState at
construction. Each represented startup event joins the machine vertical and
stillness outputs to the SAME lower frontend/sigma event. Pending boundary and
conditional TunerReady goLive preserve upstream memory. The admitted-Live
factory substitutes the same goLive result into the strong joined Racc word;
it accepts no replacement Live guard/observer/frontend/Racc/runtime snapshot.
These are conditional program relations, not proof that every admitted startup
history reaches the represented entry. Startup packet identities and clocks
are checked, but identities alone do not prove BRMM/BIAS membership.

The ordinary binary32 Eigen FromTwoVectors seed is now derived from the same
rounded accelerometer, including both normalizations, followed by the actual
first-sample Mahony update. The nearly antiparallel JacobiSVD branch is still
explicitly unqualified, not excluded from the physical source. The initialized
observer map and the independent original proof were not replaced.

Racc reads the SAMPLE-ENTRY machine TuneState sigma on every sample, before the
candidate update, even when no pending commit fires. Its frequency is the raw
preupdate WPE getter/prior, before statistics and tuning clamps. Active OU
parameters remain separately held by the pending/goLive commit mechanism.
The joined Racc and accelerometer relations consume the same machine-guard
conditioned operand. Source-uniform guard/Racc/libm supplies remain open.

### Failure analysis and independent critic

* **Physical Live input domain and totality:** arbitrary bounded residuals
  alone permitted a finite `2^80` gyro pulse that made the unchanged scalar
  wrapper nonfinite. The user-authorized MEMS API domain now imposes per-axis
  gyro 35 rad/s and accelerometer 160 m/s² bounds before execution; that pulse
  is an outside-contract regression. The new finite-lattice induction proves
  `|integral|<=4096` for every initialized finite prefix, since the largest
  actual increment rounds back to 4096 at that boundary. Corrected gyro is
  below 4132, Euler norm sum below 5000, and vertical acceleration below 322.
  This closes the scalar observer supply under seed, dormant-guard and scalar
  arithmetic premises without assuming a startup deadline or resetting Live
  memory. It does not prove finite MEKF covariance/solver execution or target
  expression correspondence. The next experiment is their actual same-state
  composition; output forcing membership cannot substitute for input totality.

* **Literal startup deadline fails on a physical yaw/wave history:** the
  circular wave `a=(A cos(omega t), A sin(omega t),0)`, with
  `A=3.313664 m/s^2`, `omega=0.64 rad/s`, has position norm 8.09 m,
  velocity norm 5.1776 m/s and centered primitive 25.28125 m s. Zero true
  bias solves every BIAS family. Continuous piecewise-linear yaw uses rates
  0.2, 0.61, 0.41, -0.2 and 0.2 rad/s on intervals ending at 5, 125, 130
  and 135 seconds. A 0.019 rad/s bounded gyro residual fits both commissioned
  profiles. Its actual world direction mean/primitive bounds are below
  0.053 and 0.501 s, with sensor-conversion allowance included. The unchanged
  full host wrapper, constructed with default configuration and no pre-Live
  magnetic calls, is still non-Live at sample 30,002: the initialized proxy's
  literal averaged force has z=+2.84663 m/s^2. At sample 30,602 it remains
  non-Live with z=+4.35638. The actual guard remains transparent, with RMS
  below 0.000131 m/s^2. This is a finite timeout-alignment counterexample
  on the tested scalar host profile, not an enclosure failure or an all-target
  arithmetic claim. It invalidates the universal 30,002-step deadline and
  using 30,602 as a source-uniform startup-plus-word horizon. It does not
  invalidate eventual capture: this same source first reaches Live at sample
  33,447 (167.235 s). The exact source/packet audit and native regression are
  `finite_startup_timeout_alignment_obstruction.py`. The limiter is actual
  body-axis integral/attitude/LPF evolution during legal yaw changes. Critic
  alternatives: prove a later uniform recovery deadline from the unchanged
  source; construct an indefinitely noncapturing admissible history; or seek
  explicit authorization for an additional physical startup premise. Merely
  increasing the old quadratic Mahony level cannot remove this trajectory.
  The exact named scalar observer also has an inverted-axis zero-feedback
  fixed point with canceled yaw integral. Its source reachability is unproved;
  the observed timeout witness still has horizontal integral norm about
  0.123 rad/s, above the available residual cancellation budget. The dedicated
  source proof and scope are in `docs/ou3-alt-startup-timeout-witness.md`.
  Next falsifiable experiment: construct a physical, primitive-preserving
  steering path into an invariant unaligned set, or prove that actual startup
  histories avoid those sets and enter the aligned branch by a later uniform
  deadline. Larger deadline scans alone cannot establish either conclusion.

* **Continuous startup budget is not a sampled budget:** a bounded 0.5 Hz
  physical wave plus continuous 200 Hz sensor residual satisfies the selected
  continuous mean/primitive limits, while the sampled mean exceeds 0.10.
  The rational certificate is in `finite_startup_direction_sampling.py`.
  This is a premise-transfer failure, not proof of failure to capture within
  150 seconds. It invalidates unchanged all-time reuse of the continuous PI
  budget. The limiter is the sampled forcing/primitive relation. Critic
  alternatives: an explicit quadrature/aliasing bound from the declared
  source, direct discrete observer/LPF capture, or an authorized sampled
  observability premise. Next test must bound the actual 30,002-sample
  alignment predicate, not merely extend the old regional Mahony metric.

* **Standalone target math link placement:** a numerical-link probe using
  the SDK ROM script failed with `windowed longcall crosses 1GB boundary`
  because the default linker placed its text outside the MCU instruction
  segment. This is an infrastructure failure, not a libm inequality failure.
  Keep the actual ROM script and place probe text at the ESP32-S3 flash
  instruction address 0x42000000; then inspect symbol resolution. The full
  firmware remains a separate composition obligation.
  The second default-script attempt moved only `.text`; `.text.expf` and
  other function sections stayed outside the instruction segment and repeated
  the relocation failure. Freeze default-script tweaking. Critic alternatives:
  full SDK memory/sections plus runtime; an explicit standalone script placing
  all `.literal*`/`.text*` in the real instruction segment; or a full Arduino
  CI build/map. The selected static namespace probe uses the second, retaining
  the SDK ROM script and explicitly withholding whole-firmware qualification.

* **Target selection is not target arithmetic:** the verified pinned MCU GCC
  uses `-ffp-contract=fast`; its numerical-header object contains ordinary and
  multiply-add instructions. Scalar no-contraction identities and two abstract
  WPE modes do not establish expression-specific target correspondence.
  This is a deployment qualification gap, not instability evidence. Verified
  Eigen/libm/compiler hashes retain program identity. Remaining limiters are
  instruction semantics, final link resolution and all-input library/solver
  error bounds. Alternatives are an instruction-level relation, source-level
  compiler correspondence with retained expression choices, or a qualified
  operation family retaining shared state. No source/firmware change is
  authorized. The Qaxis general envelope has a proved sufficient absolute
  exp error budget `2^-24`; next work must establish the actual target error
  against that budget and WPE's separately required same-argument relations.

* **FCR initialization is not established by the target ISA alone:** the
  pinned ESP32-S3 ISA reference leaves the floating-point control register
  reset value undefined, and the first-FPU-use task path inherits the active
  CPU FCR rather than writing round-to-nearest. The SDK archive scan found the
  lazy save/restore writers and no startup or ROM writer that supplies RM=0.
  The target scalar, libm, and WPE certificates therefore remain conditional
  on `FCR.RM=0` throughout the admitted execution; they cannot be promoted to
  whole-firmware correspondence. This is a documented failed proof method,
  not a claim that the hardware runs in a non-RNE mode. Critic alternatives:
  qualify the complete Arduino/FreeRTOS startup and both-core context path,
  add an authorized explicit FCR initialization in the shipping program, or
  retain the rounding-mode premise as an execution assumption. Do not infer
  RNE from reset, a zeroed task stack, or absence of `fesetround` calls.

* **WPE frequency-ratio contraction:** pinned target compilation emits
  `msub.s` for `ratio_sq-lambda_*lambda_`, while the current moment graph
  rounds the product separately. This is a missing program branch, not a
  stability failure. It invalidates claiming that the scalar moment relation
  represents every compiler expression choice. The limiter is the same-state
  omega-squared operand relation immediately before its positivity test.
  Retain both legal expressions with their exact common operands and test a
  near-threshold case; do not require agreement of the resulting branches.
  The alternative is exact instruction interpretation, which still needs
  hardware semantics. Choosing a no-FMA compiler flag would change deployment
  and is not the selected correction.

* **Machine covariance is not exactly Joseph-nonincreasing:** unchanged
  construction-rooted 600-step H18/A21 diagnostics with admitted raw inputs
  (constant/alternating 35-rad/s gyro or a slow 150-m/s² acceleration ramp)
  stayed finite with a dormant guard. The positive eigenvalue of
  `P_after_Joseph-P_before` reached 6.98e-6 (H18) and 4.69e-7 (A21).
  The source reset is `G=I+0.5[delta_theta]x`, so its norm is not one; a
  measured A21 reset norm ratio reached 1.00001156. This invalidates a proof
  that silently transfers exact-real Joseph nonincrease or orthogonal resets
  to the machine recurrence; it is not a theorem instability result. The
  limiter is source-correlated covariance/solve/reset roundoff. Critic routes:
  the gyro-bias marginal (untouched by attitude reset) with an explicit
  roundoff charge, a coupled gain/innovation energy inequality, or an exact
  machine block invariant. Next derive the marginal bound with a quantified
  charge before using it to bound corrected gyro; do not use sampled maxima
  as universal ceilings. The abrupt acceleration step engaged the guard and
  supplies no dormant-branch evidence.

* **Machine MEKF recurrence:** the strongest represented Live wrapper starts
  each deployment prediction from its exact-shadow MEKF predecessor; the
  previous machine prediction/accelerometer successor is returned as a local
  supply but not retained as the next machine predecessor. This is a
  composition failure, not evidence that the shipping filter is unstable.
  It invalidates inferring a deployed 600-step recurrence from 600 local
  supply checks, even if their individual bounds are established. The limiter
  is machine-state continuity across IMU, MAG and HOLD boundaries. Critic
  alternatives are a persistent product recurrence, an explicit finite-word
  composition theorem with checked state-equality premises, or a machine-only
  evaluator with a derived exact shadow. The next falsifiable check must reject
  a second event rooted at the old exact shadow after a nonzero machine defect.

* **Subnormal arithmetic domain:** `finite_binary32_arithmetic.rn32(2^-149)`
  rejects a finite IEEE binary32 value although the admitted source has no
  positive lower bound on nonzero components. This is a proof-model totality
  defect; normal-input identities and shipping stability are not falsified.
  Excluding small components cannot repair the universal quantifier. The
  limiting quantity is the missing subnormal lattice, not enclosure width.
  Critic alternatives are full gradual-underflow rounding, a separately
  qualified target flush-to-zero relation, or migration to the existing
  Mahony bit kernel. Extend the finite lattice and check signed ties, the
  normal/subnormal boundary and API/EMA/stillness composition against the
  independent bit kernel. Target arithmetic qualification remains separate.

* **Document-build infrastructure:** the full paper's LaTeX build stops before
  reading the theorem text because `IEEEtran.cls` is absent. This is an
  environment failure, not evidence against the new inequalities. Keep the
  paper class unchanged; validate the edited theorem section with an isolated
  standard-class harness and report the full-paper limitation separately.

* **Commissioned startup sensor contract:** the user selected deterministic
  errors throughout startup and authorized datasheet-informed numerical
  assumptions. Use residual vector caps 0.30 m/s^2 and 0.02 rad/s for BMI270,
  0.50 m/s^2 and 0.03 rad/s for MPU6886. These are commissioning requirements,
  not hard guarantees inferred from RMS or production statistics. Preserve
  every BIAS0/1/2 history and its 0.13 m/s^2 component envelope; calibration
  actually applied must justify admission, not a hypothetical correction.
  Keep the physical acceleration cap 8.8 m/s^2. Use explicit total measured
  direction-error mean/primitive budgets 0.10 and 1.5 s, separately from raw error
  magnitudes; neither is inferred from an arbitrary bounded residual.
  The scalar seed audit gives norm floors 0.48065/0.28065 m/s^2 and tilt
  bounds about 71.99/76.26 degrees before arithmetic. The retained old
  Mahony metric needs seed levels 2.05463/2.27142, above its 1.7689 outer
  level and 1.87267 chart ceiling even before enlarging temporal forcing.
  This is a coverage failure of that conditional certificate, not a shipping
  instability or justification to reduce the sensor caps. Do not promote
  that old invariant to the new profile. Critic alternatives: a different
  nonlinear observer storage; a same-history finite-time source/LPF argument;
  or a longer informative-interval recovery theorem. Merely raising the old
  quadratic level cannot fit its chart. The selected execution materializes
  the source contract, proves both seed-norm margins, and bounds the persistent
  WPE recurrence for every |vertical|<=32 sequence. Exhaustive mantissa-cell
  coverage of the literal inverse square root, including zero/subnormal norm
  inputs, gives ||q||²<1.112 after a defined scalar Mahony update. Its actual
  same-packet projection supplies |vertical|<32 without a tilt-accuracy premise.
  WPE moment and log induction margins are positive; raw period lies between
  2^-57 and 2^16 and |log_period|<=48. These close two subordinate supply
  components, leaving 11 master prerequisites. The current limiter is totality
  and capture of the preceding source-driven observer, plus actual target
  arithmetic. Next falsifiable experiment: use the new sensor/temporal budget
  in a nonlinear or finite-time capture argument that covers both physical
  seed cones and produces the literal gravity-alignment predicate. Do not
  retry the rejected level expansion or infer device admission from RMS data.

* **Magnetic frame prerequisite audit:** the small accumulation-to-handoff
  frame premise belongs to the conditional 0.63-rad entry certificate. It is
  not needed to define the current four-chart finite word. The earlier
  tilt-only substitution failed because physical heading remains in the frame.
  Alternatives are (1) prove all-source heading/lag accuracy, (2) delay the
  word until separately proved accurate magnetic capture, or (3) retain the
  full frame, arbitrary nonzero attitude and actual reference discrepancy in
  the finite graph. Choose (3), already implemented by the atlas and magnetic
  product. For proper rotations, ||A_i-G_L||<=2 universally; the existing
  magnetic envelope gives mean perturbation <=157 uT without any small-angle
  premise. This is a finite image bound, not a useful contraction margin.
  Remove the obsolete accuracy gate from the pre-rho master, retain the
  conditional accuracy lemma for later basin/usefulness work, and test gauged
  180-degree fresh entry with its unchanged full covariance and physical
  attitude. Target arithmetic and universal startup remain separate gates.

* **Startup disturbance admission:** an arbitrary finite IMU residual bound
  cannot imply the raw magnitude needed by the first-sample seed. For a level
  stationary boat with zero BIAS0/1/2 bias, let epsilon=2^-11 m/s^2 and
  n_a=(0,0,g-epsilon). Then a_raw=(0,0,-epsilon), its gravity direction is
  correct, and the residual is bounded, but the literal `norm > 1e-3f` seed
  predicate is false forever. The unchanged public wrapper remains unseeded
  and non-Live through 30,602 samples. This falsifies extending universal
  startup to arbitrary bounded IMU residuals; it does not falsify a specified
  small-disturbance startup theorem, the original proof track, or Live ISS.
  `BoundedHistory` currently binds executed post-Live forcing and supplies no
  startup raw-input observability premise. A label or finite W cannot fill
  that gap. The limiting quantity is the source-produced conditioned-accel
  magnitude, before any seed-angle or timeout-alignment proof can apply.
  Critic alternatives are (1) derive a startup magnitude/direction margin
  from an explicit sensor/model residual contract, (2) state startup separately
  from arbitrary-bounded-input Live ISS, or (3) prove eventual informative
  samples from a declared recovery/persistence contract. Neither a longer
  timeout nor more Mahony enclosure refinement addresses this obstruction.
  The exact constant-input guard/observer induction, source identity and
  native public-API correspondence pass. The commissioned startup profiles
  exclude this example through their actual residual caps and prove a uniform
  seed-norm margin. The master still rejects conflating the startup and Live
  disturbance quantifiers; initialization magnitude alone does not prove
  eventual source-produced alignment.
  The exact guard induction exposed a separate implementation defect:
  `_sqrt(0,0)` delegated to a positive-only rounding cell and rejected literal
  zero detector RMS. This invalidates that arithmetic branch's coverage, not
  the source family or shipping. Handle the exact zero-root identity before
  positive cells, reject nonzero witnesses at zero, and re-run the same
  constant-input induction. No norm threshold or source bound is changed.

* **WPE branch attachment:** the preceding machine product rejected
  any exact/machine usable-latch mismatch and the lower ledgers also required
  simultaneous log initialization/production. This is a composition restriction,
  not a proved property of every admitted source. Native boundary checks cannot
  establish a uniform positive comparison margin. Three alternatives are:
  (1) prove such margins for all WPE comparisons; (2) replace branch mismatch by
  an additive perturbation envelope; (3) compose each literal branch separately
  on the same source history. Choose (3): it retains the discontinuous control
  and exact residual, whereas (1) has no source margin and (2) loses the program
  decisions. The controlling prerequisite is equality to the finite shipping
  word before rho, not a local contraction bound. The falsifiable check is that
  different machine/exact initialization, production and takeover decisions
  survive startup/Live composition without accepting detached log operands,
  vertical inputs, frequencies, or successor states. The raw exact getter is
  also bound before either clamp; matching final clamped values cannot justify
  a different exact input to the lower event. Target accuracy and full
  execution totality remain separate obligations.
* **Frequency input lattice:** the getter and statistics-clamp relation used
  the normal-only arithmetic predicate even though these operations only
  select and compare their input. Host `expf(-100)` returns subnormal bits
  `0x1b`; the wrapper selects that positive value after takeover and the
  statistics clamp raises it to its normal floor. This is a component-domain
  omission, not a reachable-source or stability counterexample. Reuse the exact
  full binary32 lattice for getter/input validation and retain normal arithmetic
  after the clamp. Check actual-header selection on NaN/overflow/underflow,
  subnormal and ordinary outputs. No target libm accuracy is inferred.
* **Inherited ALT CI:** the finite-identity job passes. The inherited suite has
  five failures and nine errors, including unqualified phase-1 storage,
  covariance bounds and source admission; it must not be repaired by promoting
  those prerequisites. This does not invalidate the finite component identities.

* **Python quality gate:** unused module-level `dataclasses.replace` imports
  in two WPE regression modules caused F401/F811 findings; the frequency tests
  already import the helper inside the methods that use it. This is test-code
  hygiene, not a proof or shipping defect. Remove the unused imports and run
  the unchanged Python quality gate plus the affected regressions. Do not
  suppress lint rules or change any theorem qualification to clear this gate.

* **WPE source-order and takeover audit:** the earlier moment graph replaced
  `(alpha*v)*v` by `alpha*(v*v)`. With binary32 alpha=8589935/34359738368,
  previous second moment=8589935/8589934592 and v=9369095/8388608, the literal
  uncontracted result is 2816655/2147483648, outside the old graph's singleton
  11266619/8589934592. This is an implementation/correspondence defect, not a
  physical-source or stability counterexample. Separately, the persistent WPE
  machine product previously omitted `usable_period_` and its post-log-update comparison.
  Critic: more bounds on the old graph cannot certify shipping. Retain the
  literal multiplication order and source-produced latch per compiler history;
  do not infer machine takeover from the exact shadow. Next checks: native
  second-moment correspondence, inclusive takeover boundaries, latch retention,
  and disagreement between compiler histories without splicing their logs.
  The corrected graph retains literal products and both machine latches. Native
  header checks cover 200 second-moment updates and 18 latch boundary cases.
  The full startup/Live product now projects its carried machine WPE entry
  into frequency selection and checks each log successor against the same
  moment update. Initialization, production and takeover can differ from the
  exact shadow and between machine histories. The eager frequency getter and
  invalid-result fallback are represented; target correspondence remains open.
  Component regressions exercise different decisions through the full
  frontend/tuner/Racc join and reject a frontend built for the wrong frequency.
  They do not assert that their component predecessor is a reachable startup.
* **Frequency range closure:** the final outer clamp maps every finite input
  into its configured interval and nonfinite inputs to the floor. Interval
  subtraction bounds machine-minus-shadow independently of WPE accuracy or
  branch choices. The accepted statistics-update path has the tighter image of
  its inner interval under the outer clamp. These are range proofs, not a
  complete-execution or small-gain proof; the exact residual remains attached.
  The controlling limiter is still the source-uniform machine word, so this
  result does not authorize a storage or rho search.

* **Authorized counter repair and composition audit:** the user explicitly
  authorizes fixing shipping signed-counter overflow. The inner attempt and tuner accepted/rejected counts saturate
  at INT_MAX; measurements, statistics and unlock checks still execute. The induction
  c_n=min(INT_MAX,c_0+n) proves safety without a maximum event rate and preserves
  every configurable threshold, including later increases. The native saturated
  boundary passes the signed-overflow sanitizer and continues measurement and
  bias release. This removes the old counter contradiction only.
* **Remaining source-composition defects:** equating the mathematical event
  count to the shipping count becomes false after saturation; use the exact
  saturation projection. The source-bound magnetic adapter duplicates the
  lower interleaver and omits its later-north service-clock initialization. Its
  startup source constructor also demands a magnetic call, excluding the
  already represented no-pre-Live-call timeout branch. These are implementation
  and coverage defects, not physical counterexamples. Critic: consume the
  shared event composer and the certified empty prefix instead of maintaining
  a second clock update or inventing a pre-Live service premise. Next check:
  saturated MAG-to-IMU composition and source-owned ungauged IMU-to-north-to-IMU,
  preserving the same source, bias and Live/S origins.
  The first source-qualified test exposed the second forcing adapter reading
  `UngaugedLiveState.active` before north exists. A waiting calibration sample
  has a physical source packet but no inner measurement forcing. Extract that
  case once for both adapters; do not fabricate an active magnetic reference.
* **Source-audit gate after the authorized fix:** the inherited BRMM sensor
  argument rejects the changed wrapper hash. This is an expected provenance
  gate, not a mathematical failure. The full shipping diff consists only of
  `<limits>` and the guarded counter increment. Neither reads physical p/S or
  changes sensor inputs; the saturation comparison equivalence preserves every
  representable unlock threshold. Re-audit that dependency and retain the
  exact-hash gate, with a regression reconstructing the preceding source bytes.
  The independent magnetic schedule parity check also requires the obsolete
  unchecked increment. Its failed boolean invalidates that source audit, not
  the count/time release argument. Critic alternatives for this repeated
  provenance issue are a complete-source hash plus manual re-audit, a normalized
  method-body audit, or native behavior regressions alone. Use the first with
  an explicit saturation/threshold argument; native checks alone cannot prove
  all finite counts. Revalidate the original magnetic assumptions, retaining
  its north-lock premise and independence from ALT.

* **Ungauged Live composition gap:** the old interleaver rejects a timeout
  entry without north. Shipping continues initial acquisition after Live, but
  its `attitudeReferenceQuat_()` now selects the MEKF quaternion. Continuous
  hard-iron accumulation and refinement independently keep the private proxy
  tilt (`startupProxyTiltQuat()`). Reusing the startup frame for all three
  operations would violate shipping's exact operands. This is a missing
  runtime branch, not a failure of the four-chart cover or physical model.
  Critic: preserve each operation's actual frame instead of one shared frame
  chosen for convenience. Next check: distinct MEKF/private tilts during
  delayed initial acquisition, with continuous statistics accumulated once
  and the same packet entering the MEKF only after north is set.

The controlling missing object is the source-uniform finite deployment word,
not a storage margin. No rho/storage search has been run or authorized.

* **Pre-rho entry representation failure:** the timeout does not require
  `isTunerReady()` or magnetic north. A quiet, level, zero-bias source has
  identical IMU history for every constant true heading. With no pre-Live
  magnetic calls the proxy hands off identity, including when true yaw is pi.
  The fresh joint24 definition then divides by the zero scalar component of
  `q_true_WB * conjugate(q_hat_WB)`. This invalidates universal fresh entry in
  one Cayley chart; it does not establish shipping instability or invalidate
  the original gravity-quotient track. The exact limiter is a zero chart
  denominator, not a loose enclosure. Test the unchanged wrapper from reset
  and verify the analytic constant-heading family and its unbounded Cayley
  limit. A TunerReady-only bridge also misses the literal timeout edge.
* **Finite-horizon coverage defect:** exact RN32 accumulation of 0.005 first
  reaches the default 150-second comparison at sample 30,002. The following
  600 IMU transitions end at sample 30,602, outside the 30,600 machine-history
  cap. Even this deadline remains conditional on the gravity-aligned branch;
  merely enumerating the clock does not prove universal handoff reachability.
* **Magnetic premise mismatch:** the conditional `E<=8.5 uT` calculation
  substituted a gravity-direction bound for a full accumulation-to-handoff
  frame bound. Even perfect tilt leaves physical heading in a yaw-stripped
  sample. The exact north/south pair for B=(15,0,20) has mean (0,0,20), so
  tilt alone cannot prove nonzero mean horizontal field. This falsifies the
  premise substitution, not the conditional chord/mean inequality and not a
  complete shipping capture trace. Retain heading excursion and handoff lag
  in the frame relation; sampled means cannot qualify it universally.
The critic alternative to further seed/SVD bounds is to address the startup
representation first. Three distinct choices are (1) carry a gravity quotient
and neutral yaw until source-qualified magnetic regauging, (2) retain a
homogeneous quaternion or overlapping attitude charts with exact transport,
or (3) supply a separately proved all-history magnetic capture theorem before
entering the current chart. An unproved startup magnetic-cadence/heading
restriction is not choice (3). None may shrink COMPLETE-BRMM, change shipping,
or replace unknown yaw by an independently bounded Cayley disturbance. The selected four-chart architecture now covers every nonzero relative
quaternion and carries exact finite event transport; the conditional ungauged
Live composer retains initial acquisition and its later north transition.
Universal source/control/deployment qualification still keeps the guard blocked.

* **ALT implementation/order defect:** holding Racc sigma between pending
  commits contradicts shipping's per-sample pre-candidate TuneState read. The
  limiting quantity is that exact scalar dependency, not a norm bound. The
  wrapper now reads the carried predecessor machine sigma unconditionally;
  a no-pending changed-candidate regression distinguishes it from the former
  held readout. Reading post-clamp frequency was likewise corrected to the
  raw preupdate getter. These failures invalidate those old dependency claims,
  not shipping stability, physical admission or the scalar arithmetic lemmas.
* **Coverage/method gap:** the initialized-only Mahony map could not execute
  literal reset. The ordinary source-bound seed now removes the free quaternion
  port. Unrepresented SVD/nonfinite branches, target sqrt/profile qualification
  and source-uniform supplies still prevent universal startup coverage. A
  30,602-sample bookkeeping cap accommodates the conditional startup+word graph;
  it is neither proof of a 150 s startup deadline nor clock-lifetime closure.
* **ALT regression defects:** old joined fixtures supplied an independent zero
  band input instead of the same Mahony output; runtime/readiness accessors used
  invalid nested paths/keys. Tests now consume the same source outputs and keep
  the rejecting constructors. A new scheduler-splice test initially used an
  invalid predecessor and was corrected to a valid, different service credit,
  so the goLive ancestry check itself is exercised. No guard was weakened.
* **Inherited certificate/conditioning failures:** the broad suite still fails
  at `event_local_covariance_norm_envelope_closed` and dependent finite-master
  prerequisites. The magnetic AD lineage also raises `AD derivative dimensions
  differ`; an endpoint partition exact-equality assertion observes
  `0.49999999999999983` rather than `0.5` under outward arithmetic. These are
  separately classified proof-construction/test failures, not permission to
  narrow outward intervals or assert the theorem false. Original-track modules
  remain untouched.
* **Inherited build-gate failure:** `make all` stops at
  `tests/kalman_ou_iii/live_entry_audit-test`: the quiet pair does not execute
  H18-to-A21 release by its tested endpoint, and two quiet handoffs have
  session-origin S=288.1474364 rather than the asserted S>300. The terminal
  result is `LIVE_ENTRY_AUDIT_PASS=false`, `Makefile:56: run-tests Error 1`,
  then `Makefile:36: test Error 2`. This is an unchanged original-audit assertion
  failure, not evidence against the counter induction or attitude identities.
  Critic: changing thresholds to make this PR green would obscure the original
  entry question. Retain the failure; rederive those assertions against the
  corrected source in that independent track. The limiting quantities and next
  falsifiable check are the actual release time and session-origin S at handoff.
  The same run also reports the pre-existing OU-III simulation quality gate
  failures `accel_z_bias_percent_exceeded` (4.3921% and 4.46846% versus the
  4.3% limit) for the two high-wave cases. Those numerical regressions are
  outside this proof PR; no tolerance or shipping behavior was changed here.

The strongest critic objection is that an identity-preserving startup graph can
still carry unqualified source/branch/arithmetic witnesses; conditional entry
attachment is not universal finite capture. Three materially different routes
remain: (1) complete source-bound seed/SVD and branch relations and compose the
existing regional certificate, (2) generate a single source-order SSA/interpreter
relation from shipping, or (3) independently verify a native interpreter with
explicit compiler/library semantics. The current continuation uses (1) where
actual modules exist. More seeds, scalar enclosure subdivision and storage/rho
search do not address this missing relation and are not authorized.

### Evidence and validation boundary

The new startup tests execute consecutive Cold samples from literal reset;
conditional goLive/admitted-origin tests separately verify memory, pending and
scheduler substitution and reject splices. Their TunerReady fixture tests an
implication, not an asserted capture trajectory. The first-seed native test
compiles unchanged shipping headers with scalar Eigen/no FMA and observes
public first and second updates without installing a quaternion. Host checks
are correspondence regressions, not target ESP32 profile qualification.

The unchanged current `live_mahony_regional` certificate validates with discrete
margin `1.518757454127882e-4`; its startup-entry premise remains open. This is
consistent with the newer shared two-level repair, not the stale claim that the
continuous regional invariant or 35 deg/s PE refinement is still broken.
Private observer accuracy, total startup accumulation-frame heading excursion,
ungauged timeout and source-qualified deployment seed remain distinct issues.

The full inherited suite still has covariance-envelope, endpoint-partition and
magnetic AD failures. It is not claimed green. Latest focused/native,
finite-suite and build results are recorded in the PR validation report; none
promotes source membership, capture, retention or stability.

### Current limiter, DEAD_ENDS and next falsifiable work

The four-chart runtime, conditional ungauged continuation and source-owned
empty-startup/later-north composition close their representation/event gaps.
All magnetic counts now have source-bound saturation safety on every finite
prefix; the master consumes this fact. The current limiter is universal
startup source/control reachability, full magnetic history and target arithmetic.
Every nonzero quaternion is covered, but that does not certify a storage basin
or finite capture. Storage must not identify different chart origins as the
same zero physical error.

Qualify the machine guard displacement into Racc/accelerometer, raw WPE
exp/log/sqrt, band/statistics/tuner, Q-axis, clocks, trig/normalization, Eigen
solves/floors, nonfinite branches and comparisons on every literal same-history
IMU/magnetic/hold/reset edge. Preserve the saturated shipping count separately
from the mathematical event ordinal. The shared 30,602-sample budget covers the
first default timeout comparison plus 600 IMU edges only conditionally on the
actual aligned-branch predicate.

The master reports no falsified counter prerequisite; universal entry and
arithmetic qualifications remain open. Every-prefix retention, joint24
storage/rho, an ultimate bound and no-restart indefinite tiling remain later
obligations. No larger clock horizon establishes missing startup capture.

Frozen shortcuts: derivative cocycles presented as finite maps; frozen gains,
replays or more seeds as universal admission; independent coefficient boxes;
A21 18-state marginal storage; covariance-consistency entry sets; artificial
physical-domain/basin reduction; wordwise S resets; requiring a single generator
for COMPLETE-BRMM. The two-strike rule and `assert_finite_storage_master` remain
in force.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`,
`storage_search_allowed=false`.
