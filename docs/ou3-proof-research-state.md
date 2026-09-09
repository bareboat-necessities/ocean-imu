# OU-III proof research state

## Current hypothesis

BRMM is the primary physical source. The target is bounded accelerometer-bias
error and regional practical ISS of the other 18 errors, with the full
21-state shipping estimator in active mode. Retain one physical history,
actual anisotropic R_S and physical -S_true, all Q/floors/cross covariance,
finite reset/projection, and the closed .4 estimate ball. P3 delta remains
1e-18. Lever arm is disabled and the dormant vibration guard is transparent.
No replay fitting, filter change, reduced operating domain or quality-gate
change is permitted for proof convenience.

The exact sample-entrance motion reduction is proved in
`docs/ou3-p4-bounded-bias-cocycle.md`. Compose the full 21-state subevents
within each sample before selecting motion rows. On the actual joint history,
`x_next=A_i*x_i+G_i*b_i+d_i` holds with corrected active feedback retained.
The bias is charged only at sample entrances, as the theorem specifies.
An exact cascade is not required when the actual internal bias is bounded.
The resulting motion cocycle is generally not the principal motion block
of the full word product. Neither construction changes the shipping filter.

## Evidence

`ou3_p4_bounded_bias_cocycle.py` consumes the same hash-attached PM 1.5 m,
28 ft RAO source capture and both fixed 600-sample windows. It constructs a
Stein metric without an eigenbasis/diagonalizability premise and reports
maximizing directions, 80-digit coefficient checks and signed operation costs.
These are numerical point diagnostics, not outward nonlinear certificates.

| Quantity | H18 | A21 |
| --- | ---: | ---: |
| Sample motion cocycle spectral radius | .997663960671 | .960675313400 |
| Constructive quadratic endpoint ratio | .995550332926 | .953993263418 |
| Worst completed-prefix ratio | 11.2178082213 | 1.63544556380 |
| Endpoint supply rho | .997775166463 | .976996631709 |
| Separate bias-energy gain | 4847.73017438 | 117.303616700 |
| Separate common-template gain | .237245019053 | .025801527123 |

The active principal-word block and sample cocycle differ by operator norm
.516227538638. The lift driven by the actual full-factor bias sequence agrees
to below 1e-15 in both modes. The source remains one correlated template;
its unit amplitude is not a uniform BRMM budget. The numerical metric
condition numbers are 3.33e9 / 5.11e9, so directional high precision and
small Stein residuals do not constitute rigorous enclosure.

### Bias-family driver recurrences

All three mandatory accelerometer-bias families now have their own authoritative
driver recurrence, each from its own declared parameter box in
`ou3_p4_closure_domain.json`; none is read off BIAS1. Every admitted history
satisfies the exact relation `b_i = phi_true*b_{i-1} + w_i`, so the content is
the outward bound on `w` and on the tau-mismatch forcing `m=(phi_true-phi_hat)b`
against the deployed `phi_hat=exp(-.005/5000)=.999999`.

| Family | Driver | Admitted `phi_true` | Driver norm | Mismatch norm | Joint supply |
| --- | --- | --- | ---: | ---: | ---: |
| BIAS0 | composite: turn-on, thermal, strain, non-GM, declared pathwise GM cap | [.9999833335, .9999986111] | 8.7552e-4 | 3.5276e-6 | 8.7553e-4 |
| BIAS1 | one root, one parameter history | [.9999916667, .9999979167] | 4.9266e-6 | 1.6512e-6 | 5.1960e-6 |
| BIAS2 | bounded-variation drift, no relaxation root | [.9999979167, 1] | 2.2011e-6 | 2.4393e-7 | 2.2146e-6 |

All three also reach the same-history projection/Joseph graph, so the
bias-family half of the P4 contract is closed. The prerequisites see a family
only through admission, one retained physical bias history, the shared `w`
column and `|b_true|`: the radial projection map `F_R(e,beta)=beta-Pi_R(beta-e)`
carries no driver term and the Joseph/reset gains come from the reachable
`P/H/R` cell. The three declared boxes share one true-bias envelope to 5.6e-17,
so the fourth dependence is one number and `|e_b| <= R+B_true = .6252` holds
for all three.

A uniform BIAS2 separation constant is NOT a prerequisite of the declared
objective. Bounded bias error comes from the closed radial projection sector,
which holds for a non-relaxing truth exactly as for a relaxing one, and the
motion half is the cocycle contraction against a persistent forcing, i.e. a
finite ultimate bound. `mu_sep` sharpens motion gains only; it stays unproved
and is no longer counted as a blocker.

Each family's factor interval passes through the deployed 24-state event lift
with the shared `w` column retained, and each declared true-bias norm 0.2252
stays inside the .4 hard-entry radius. On the ISS pair (factor interval, driver
bound) BIAS0's class contains BIAS1's; BIAS1's contains neither of the others
and no single family covers the other two, so the three-family quantifier
cannot be discharged by one proof. BIAS0's pathwise GM increment cap is a
declared family hypothesis at 4.33 times the largest admitted one-step
innovation sigma; it is not a consequence of the stationary PSD, which supplies
no pathwise cap at all.

Exact-real projection already preserves the estimate ball in both modes.
The prior .400126 A21 bias-row value extrapolates frozen interior projection
coefficients; it is not an actual nonlinear projection escape. BIAS0/1 plus
projection gives `B_error <= B_true + R`, without assuming bias convergence
or zero active motion-to-bias feedback. The captured binary32 radius is
.4000000059604645. Shipping rounding remains a separate enclosure obligation.

### Which declared entry ball actually limits retention

`entry-block-retention.json` switches each declared entry ball on ALONE and
maximizes over every completed prefix of the same attached capture, so the
subadditive total that `domain-retention.json` reports per coordinate is split
into the terms that produce it. The capture is reproduced bit-for-bit from
`ou3-source-endpoint.cpp` against the pinned generator
e442150682f560384be427df4cc7815956a091c5; all four capture hashes match the
retained record, as do every number the two experiments share.

Two facts follow that the total cannot show.

The accelerometer-bias entry ball limits **no** coordinate, not merely not the
total: its worst single-ball reach is .1544 (H18) and .2335 (A21) of a radius
over both modes and all six motion coordinates. Together with the .2382
same-graph result this closes the bias entry ball as a retention limiter.

The declared 300 m*s integral-displacement ball is the worst single ball for
five of six H18 coordinates and four of six in A21, reaching 32.50 velocity
radii in H18 on its own. With it the H18 attitude coordinate reaches Cayley
norm 4.5788 against the declared chart bound 1.0, so the declared product box
drives the state out of the chart the frozen map is expanded in and that
configuration invalidates its own linearization. Without that one ball the
chart is retained in both modes, at .9411 (H18) and .6537 (A21).

The minimal working radii this word retains at every completed prefix, as
multiples of the declared entry radii and with the independent integral ball
removed, are

| Coordinate | H18 | A21 |
| --- | ---: | ---: |
| attitude | 1.756 | 1.220 |
| gyro bias | 2.315 | 1.428 |
| velocity | 4.033 | 4.921 |
| position | 1.027 | 2.237 |
| integral displacement | 1.000 | 1.000 |
| latent acceleration | 1.935 | 2.199 |

These are a working domain strictly LARGER than the entry set, never a reduced
entry set: the entry radii are unchanged and the retention target is enlarged,
which makes every downstream nonlinear obligation harder. Only attitude carries
a chart constraint, and it is the row that already fits.

These ratios are normalized by the ENTRY radii, which is the producer's
convention and not the theorem's requirement. `thm:brmm-bounded-bias-motion`
already separates the entry level `L` from a chart-valid level `L_chart` and
asks only for `Gamma*L+C_p < L_chart`, so a prefix excursion above an entry
radius is not by itself a failure. What the excursion must not do is leave the
chart on which the nonlinear majorants are declared.

That is where the declared box fails, and it fails the theorem's own
hypothesis rather than a bookkeeping convention. With `L` the declared product
box the H18 prefix excursion reaches Cayley norm 4.5788 against the declared
chart bound 1.0, so no chart-valid `L_chart` satisfies `Gamma*L+C_p<L_chart`
and the hypothesis is unsatisfiable at that entry set. Removing the one
independent integral-displacement ball makes it satisfiable: the excursion is
.9411 (H18) and .6537 (A21), so a chart-valid `L_chart` exists, with 5.9%
margin in H18.

The velocity row is physically forced, and is a statement about how large
`L_chart` must be in the coordinates that carry no chart constraint, not an
obstruction: a 30 degree attitude entry error mis-resolves gravity by
g*sin(30)=4.903 m/s^2, which over the 3 s word is 14.7 m/s against a 5 m/s
entry radius, and the attitude ball alone reaches 2.9178 velocity radii in A21,
which is that number. Velocity carries no chart, so absorbing it is a
declaration, not a lemma.

### The correlated integral entry relation

`ou3_p4_correlated_entry_relation.py` replaces the independent ball with what
the deployed recurrence and the deployed S=0 scheduler actually permit. Four
facts are now materialized and outward.

The deployed translational factor's S row is exactly `(dt^2/2, dt, 1, phi_Sa)`
with `phi_Sa = int_0^dt (dt-s)^2/2*exp(-s/tau) ds` in `(0, dt^3/6]`, and the
same factor drives truth and estimate, so the ERROR obeys
`e_S_next = e_S + dt*e_p + dt^2/2*e_v + phi_Sa*e_aw` with no remainder. The
S=0 error update is exactly `e_S^+ = (I-K_S) e_S^- - K_S S_true` with
`K_S = P_SS (P_SS+R_S)^{-1}`, i.e. the physical `-S_true` forcing.

The declared 300 m*s is not derived and not conservative. The operating domain
records no provenance for it, and at the declared 20 m startup position
envelope the free integral reaches `20*T_handoff`, above 600 m*s over the
declared live-entry timing floor. It is an arbitrary constant, hence a class-D
entry-set modeling defect in either direction.

The scheduler is the missing quantitative ingredient. `T_S` is clamped to at
most .15 s, so consecutive S=0 Joseph events are at most .155 s apart, not one
per handoff interval. Over one such gap at the declared handoff envelope the
exact accumulation is at most **3.2020 m*s**; over one 3 s P4 word it is
60.1376 m*s.

| Integral radius | m*s | Chart threshold 4.5664 | Correction ceiling 10.2458 |
| --- | ---: | :---: | :---: |
| declared independent ball | 300 | fails | fails |
| one P4 word of dwell | 60.1376 | fails | fails |
| one S=0 cadence of dwell | 3.2020 | passes | passes |

The chart threshold is derived from the retained `entry-block-retention.json`:
scaling the integral radius by `s` scales exactly that term of the subadditive
per-prefix sum, so `total <= without_S + s*alone_S` at every prefix, and H18 is
the binding mode at 4.5664 m*s. It is a frozen-map diagnostic threshold, not an
outward certificate.

What is NOT established is an unconditional anchor. The S=0 update is a
non-expansion of `S_hat` in the `R_S`-weighted norm for every cell, because
`I-K_S` is similar to a symmetric operator with norm `1/(1+mu)` and `mu>=0`.
A uniform CONTRACTION needs a uniform positive `mu`, hence a reachable `P_SS`
lower bound. The only source-uniform one available is `P^->=Q` over one
prediction step, which gives `q_SS >= sigma^2*dt^7*(1-x/4)^2/(126*tau)` and
`mu >= 1.3e-26` against `R_S <= 100` -- above 1e25 events per e-fold. So the
deployed S=0 regulation cannot anchor `e_S` uniformly over the admitted cell
family, and an unconditional correlated entry set needs either a proved
reachable `P_SS` lower bound far above the process floor or the P5 capture
argument. That is reported as class E rather than assumed inside P4.

## Current limiter and failure analysis

The controlling unresolved quantity is a useful **uniform, consecutive-word
motion supply and retention bound**, not exact-real bias compactness. The
retention half is now split: against enlarged working radii it holds on this
word for every coordinate once the independent integral ball is replaced, and
the open items are the correlated integral ball, the enlarged working-domain
declaration and its nonlinear majorants, and outward uniformity.
Pointwise `rho(T)<1` and existence of a metric for each T do not imply
compatible contraction along a nonlinear source continuation. A pair of
Schur matrices with an unstable product is included only as a logical
counterexample, not as an admissible OU-III source.

Materializing the three driver recurrences moves the limiter but does not
close it. BIAS0 now dominates the supply axis at 168 times the BIAS1 joint
supply, driven by its declared pathwise GM cap; BIAS2 dominates the retention
axis, because `phi_true=1` leaves the truth with no relaxation at all and the
bias-error mode decays only through `phi_hat` and the corrections. Neither is
open on the same-history graph any more: the bridge's projection/Joseph
prerequisites are family-parametric and each family discharges them from its
own certificate, as the executable gate records in
`bias_family_source_uniform_same_history_closed`. `mu_sep` remains unproved
and is not a prerequisite of the declared objective, so it is not a blocker.
What the three driver recurrences move is the motion half, and that is where
the six remaining gate blockers now sit.

The new separate endpoint supply is feasible at the coefficient point, but
the relaxed L2 bias budget still does not establish coordinate retention.
For zero initial motion deviation, B_true=0, .4 bias radius and unit template,
its sufficient prefix bounds reach 1.3212 chart, 1.1134 velocity and 3.5159
latent-acceleration radii in H18, and 3.3272 latent-acceleration radii in A21.
These are upper bounds over relaxed coefficient inputs, not attained physical
excursions. Initial motion and repeated-word accumulation are not yet added.

Classification: proof-method / loss-of-correlation supply-enclosure failure.
It invalidates using compactness with independently relaxed entrance bias
energy as a complete retained-tube proof. It does not invalidate point Schur
stability, the exact reduction, bias compactness, correlated-source ISS, or
the deployed filter. No further scalar norm or interval refinement of this
independent-energy route is justified by the current margin.

The independent critic's strongest objection is that the bound discards the
temporal bias constraints that control the actual forcing. Compare three
qualitatively different alternatives: (1) a joint bias/source graph retaining
held-bias/physical-driver recurrence and active projection, optionally with a
proved BIAS2 sector; (2) a compatible source-dependent metric over actual
consecutive words; (3) source-centered correlated coordinate tubes with a
qualified entry set. The next falsifiable experiment is a consecutive-word
joint-graph test retaining the actual bias recurrence and source continuation,
with endpoint and every-prefix coordinate budgets before uniform covering. Run
it per family rather than once: BIAS0 tests whether the budgets survive the
largest admitted supply, and BIAS2 whether they survive a non-relaxing truth,
which is the case that decides whether a separation sector is needed at all.

## Retained facts and rejected routes

- `three-storage-routes.json`: source-centered motion endpoint ratios
  .999572658878 / .932948381611 and prefix maxima 1.003245112488 /
  1.000001212298 remain point facts. Separate scalar information minima and
  the fixed SI/gravity diagonal metrics failed; changing units cannot fix
  their ratios. A signed full-matrix information route remains point-feasible.
- `domain-retention.json`: the initial bias ball plus template retains all
  motion coordinates in the frozen full-state map, at most .2382 of a radius.
  Scalar conversion of velocity/displacement energy into worst-direction
  attitude manufactured the prior 6.41 / 2.27 chart failure. Do not refine it.
- The full independent handoff box fails that frozen-map retention test,
  limited by velocity 35.358 / 9.357. The 300 m*s integral ball is a dominant
  source. A physically reachable correlated entry set must be proved; a
  covariance ellipsoid is not automatically a true-error envelope. Its point
  requirements remain 4.974 / 27.452 initial sigma, not a certificate. The
  correlated relation derived above meets the chart threshold at one S=0
  cadence of dwell; the unconditional anchor is still open.
- `bias2-motion-gain.json`: BIAS2 was invoked on the corrected full21 graph;
  A21 point separation does not certify a uniform sector or strict master.
  Independent-port common gains and common-template/common-gain storages
  were quantitatively unusable. Neither grid extension nor interval depth
  repairs those lost dependencies.
- Scalar two-occurrence PE transport at the proposed 250 deg/s cap has a
  negative analytical margin. Tightening arithmetic cannot change its sign.
- Passive trace instrumentation must reproduce its matched baseline bitwise.
  The retained same-history attachment and all actual R_S events remain.

## Conditional P3 and physical-source scope

Conditional P3 passes H18 and A21 at delta=1e-18 under the declared
Normal-Live execution premises: acceleration <=4 m/s^2, rate <=30 deg/s,
Racc=.04 I, Rmag=.09 I, actual accepted-vector geometry/recurrence and the
complete covariance/event history. Retained H18 outward LDLT pivot:
4.987499868870966e-14; active A21 bias margin: 1.2499987189052501e-9.
The 17 conditional producer/promotion tests pass in the retained evidence.
Canonical outputs and hashes are in `reports/results/rao_stability/`.

The pinned reference is oceanography-waves-lib v1.2.1,
`sim-data-files-vessel-rao-28ft.zip`, SHA256
`6d6eb92db78e97f1d2456c6b92387993bf23be14a9a6cf26943f0a22fc6fee60`,
generator e442150682f560384be427df4cc7815956a091c5. Among 1,920,000
source samples, 10,205 exceed the acceleration cap and 31,937 the rate cap;
runtime Live leaves 10,002 / 31,589 violations. Actual replay Racc and Rmag
also differ from P3 premises. Complete deployment/source admission is false;
a sampled record or runtime Live flag is not an admission certificate.
BRMM-0 bounds same-history velocity and excludes fixed acceleration DC;
it does not alone bound position/S or imply vector PE.

Conditional mathematical work does not wait for device qualification. Its
execution/source premises and uniform nonlinear proof still must be explicit;
failed admission of these replay records is not a proof of instability for
the conditional source class. Physical qualification and P5 capture remain
separate from the post-entry mathematical implication.

## Performance evidence, unchanged by this proof work

The retained frequency, TFG, PII and selected OU-II default records pass.
OU-II's paired bias/drift study and 5 s stillness decay are documented in
`reports/results/sigma_horizon/study.md`; increasing the averaging horizon
alone did not fix low-wave accuracy. OU-III's two high-sea Z-bias gates and
NLO's low-wave PM-Stokes vertical gate remain unresolved. The nine-profile,
eight-seed study in `reports/results/rao_parameter_tuning/ou3_acc_z_bias/`
did not justify committing the default-seed-only passing profile. The RAO
propagation-to convention correction and its 11 tests are retained.
No performance coefficient or threshold is changed here.

## Which blocker is which failure class

`ou3_p4_blocker_falsification_classification.py` attributes every open gate
blocker to one of the declared classes A (an actual admissible trajectory or
source counterexample), B (rigorous source-uniform infeasibility), C (enclosure
/ conditioning / dependency loss), D (entry-set modeling) or E (missing source
qualification or materialization). Only A or a rigorous B would support saying
P4 is unprovable on the declared domain.

| Blocker | Class |
| --- | :---: |
| COMPLETE BRMM same-signal estimator/source cover | E |
| interdependent `(tau,sigma,T_S)->R_S` with the literal adaptation history | E |
| source-uniform same-cell Joseph correction/reset domain | C |
| source-uniform exact-graph endpoint augmented LDLT | E |
| source-uniform exact-graph every-prefix augmented LDLT | E |
| same exact graph every-prefix hard-domain retention | C |

No class A and no class B finding exists. The 4.5788 frozen-word Cayley value
is NOT promoted to a lower bound on the true nonlinear trajectory: outside its
chart the frozen expansion is not valid, so it is evidence that the independent
integral ball destroys the certificate, never a disproof of P4.

### The correction/reset blocker is a dependency loss, not infeasibility

The same-cell magnitude certificate is
`||d_theta||^2 <= trace(P^-_theta) * y^T R^{-1} y`. Two quantities decide it.

The integral entry ball is NOT its limiter. Feeding the correlated 3.2020 m*s
relation instead of the 300 m*s ball moves the S=0 ceiling from 9.6223e7 to
9.9481e5, and the ACCELEROMETER event then limits at 2.4091e6. Both remain far
outside the exact reset utility domain 3.0, so the entry set is not what blocks
this obligation.

The attitude covariance envelope is its limiter, and it is a conditioning
artifact. The retained endpoint-referenced envelope gives an attitude variance
upper of 3.99983e8 rad^2 per axis, 1.19995e9 as a trace. But the accelerometer
Joseph event admits a prior-INDEPENDENT posterior cap: for a scalar angle
measurement with gain `|f|` and noise variance `r`,
`P^+ = P^- r/(P^- |f|^2 + r) <= r/|f|^2` whatever `P^-` was, and the declared
Normal-Live invariant executes that update at EVERY valid IMU sample. With
`R_acc = .04` and `|f| >= 5.80665` this caps the two directions transverse to
the specific force at 1.18634e-3 rad^2 each. The envelope in use therefore
exceeds what the deployed event structure permits by **6.4485e11**, which is
the endpoint-referenced Lagrange/Vandermonde inversion losing conditioning over
its .775 s observation window, not a physical bound.

One lossy step then remains in the magnitude route itself. At the transverse
cap the `R^{-1}` relaxation gives an accelerometer ceiling of 3.3876, just above
the reset utility limit 3.0, while retaining the same-cell
`S^{-1} = (H P H^T + R)^{-1}` gives 2.3954, because at the cap `H P H^T` equals
`R` exactly and the innovation covariance is twice `R`. So the route closes on
the two transverse directions once `S^{-1}` is kept. The third direction,
rotation about the specific force, is not observed by the accelerometer at all
and needs the asynchronous magnetometer plus gyro-bias transport; that is not
established and stays open.

### Finite-precision status

The conditional mathematical binary32 additive ISS enclosure is closed and is
the only finite-precision object P4 depends on. Target-toolchain qualification
of the declared binary32 LDLT/libm forward postconditions is a separate
deployment blocker. Neither may masquerade as the other: a missing toolchain
qualification is not a mathematical P4 failure, and a host-only arithmetic
check is not deployment qualification.

## The one object the remaining blockers reduce to

Six gate blockers remain. Four consume the same missing object, the
**source-uniform COMPLETE BRMM cover**: the same-signal estimator/source cover,
the interdependent adaptation history and the endpoint and every-prefix
augmented LDLT each need a source-uniform cell family that has not been
materialized. The cover needs the estimator-owned transition operator over
every admitted BRMM continuation, every hard-entry radial segment and every
correlated Joseph cell, with the coefficient image proved inside the target
cell. A captured word cannot supply it and the contract forbids trying:
`point_trace_can_promote_source_uniform_cover` and
`trajectory_replay_or_pinned_generator_may_establish_uniform_cover` are false.

The other two, the correction/reset domain and the every-prefix hard-domain
retention, are now separately attributed. Neither is limited by the cover
alone: the correction domain is limited by an attitude covariance envelope
6.4485e11 above what the deployed accelerometer event permits, plus the
`R^{-1}` relaxation of the same-cell `S^{-1}`; the retention is limited by the
independent integral entry ball, which the correlated relation replaces on the
chart threshold. Both are class C, and both have a stated constructive repair
rather than an infeasibility.

**P4_MOTION_PASS=false, P4_PASS=false. P5-motion and P5 may not start.**
