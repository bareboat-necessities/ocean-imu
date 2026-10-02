## Complete regular-word joint zero-action kernel CLOSED (nullity <=1)

The authoritative qualitative reduction is now formalized without attitude-only Schur compression. Zero fresh process action plus four distinct applied S rows (extended-Chebyshev system for {1,t,t^2,psi_tau}) kills the independent homogeneous LIN/AW root. Zero homogeneous magnetic loss together with MAGNETIC SERVICE reduces the remaining root attitude space to dimension <=1. On active A21 the literal accelerometer Jacobian has J_ba=I, so a pure BA kernel is impossible and, once an attitude amplitude is chosen, one accepted accelerometer row uniquely determines the complete BA root; additional rows can only remove that line. Therefore every complete regular word has joint zero-action kernel dimension <=1. If nonzero it is the word-dependent physical tilt/BA compatibility line; exact existence of a line is not asserted on every word. Its persistence through qualified MOVING superwords is separately excluded by the closed A* physical source theorem. Continuous magnetic service rows are covered by their SPD aggregate Gramian and regular prediction/reset transports preserve nullity, so no callback enumeration is required. literal_event_strata_regressed is now TRUE in this qualitative sense.

## Correction: A* attitude block has nullity zero, not one

The proposed target nullity(C_theta,W)=1 is false for the aggregate attitude-only block under Corollary A*. A* strictness is precisely positive definiteness of the three attitude-column Gram, so the attitude row stack has rank 3 and nullity 0. The familiar one-dimensional field-axis tilt/BA compatibility line appears only in the JOINT accelerometer attitude+BA equations after allowing BA nuisance compensation; it is not an attitude-only kernel. Consequently the graph form [[C_theta,0],[A_W,I]] is not the physical compatibility construction when C_theta is the A* attitude block: with C_theta full rank that stack is full rank. Added a regression that fails the attitude-nullity-one claim under an A*-positive block. The literal kernel theorem must instead retain the joint acc/BA block through nuisance elimination, intersect it with magnetic/gyro/process/S constraints, and then prove that JOINT kernel has dimension <=1. Do not promote literal_event_strata_regressed from an attitude-only rank calculation.

## Canonical kernel correction: compatibility is a word-dependent graph

There is no single fixed rational compatibility vector valid for all complete words. The authoritative line is r_W=(a,-A_W a), where A_W depends on the actual nominal-force rows, resets, held/active BA transport and coupled chronology. A global fixed-vector rank regression would therefore be synthetic. The exact canonical lemma is graph-valued: after magnetic service eliminates its two service coordinates and process/S constraints eliminate independent LIN complements, let C_theta,W be the remaining attitude constraint and let A_W map attitude into forced BA. If nullity(C_theta,W)=1, then the exact stack [[C_theta,W,0],[A_W,I]] has nullity one for EVERY A_W and ker={ (a,-A_W a): a in ker C_theta,W }. Added exact rational graph_kernel_certificate and regression with arbitrary A. Thus BA graph coupling itself needs no interval cover; the remaining substantive statement is one-dimensionality of the reduced attitude constraint, followed separately by A* exclusion of persistent nonzero a.

## Continuous MAGNETIC-SERVICE/transport cover collapsed analytically

For the qualitative kernel theorem no finite magnetic row cover is needed. MAGNETIC SERVICE states that on every 1-s window the aggregate transported/whitened heading--axial-BG Gramian satisfies sum G_i'G_i >= mu_M I_2 with mu_M=1. Hence intersection_i ker G_i={0} in the two service coordinates for arbitrary admitted event count, timing and row orientation. This is a complete continuous magnetic-row kernel cover. Likewise regular prediction and reset transports cannot create kernel dimension: OU/attitude prediction is invertible on the regular domain (phi>0 and triangular nonzero diagonal), and G_reset=I+[d]/2 has det=1+|d|^2/4>0. Thus dt/tau/rotation/reset cells collapse by invertible coordinate transport. The only nontrivial exact rank regression left is the canonical transported intersection of S/LIN, aggregate accelerometer geometry, gyro/process and BA compatibility rows after magnetic service coordinates have been eliminated. This is one canonical algebraic stratum, not an event-schedule enumeration.

## Literal nullspace stratum assembly: covariance blocker removed, continuous cover still open

Added nullspace_strata.py. For qualitative zero action, correction gains and innovation covariance are not needed: exact Joseph zero loss implies H e=0. The literal stack therefore consumes only process/sync homogeneous rows, S rows, applied magnetic Jacobian rows, accelerometer Jacobian rows, gyro/process rows and BA compatibility rows, all transported to common root coordinates. Exact rank/nullity is then checked by complete_word_nullspace.py. This removes the open full-P seed blocker from the qualitative kernel lemma. However the current MAGNETIC SERVICE contract has no finite callback-pattern cover, and the repository still lacks a theorem-domain continuous cover of the admissible applied magnetic/service-row transports. Therefore literal_event_strata_regressed remains FALSE. The remaining task is specifically to construct that continuous service-row/transport cover (not enumerate schedules), then feed each interval/rational stratum to the exact kernel certificate.

## Literal nullspace strata: scheduler enumeration replaced by service-row parameterization

The existing magnetic seed audit proves there is no finite callback-pattern cover under MAGNETIC SERVICE: accepted information is constrained in aggregate and regular callbacks may occur at arbitrarily many 4--6 ms patterns. Therefore "enumerate every reachable event stratum" cannot legitimately mean enumerate scheduler bit patterns. For the qualitative zero-action lemma gains/full P are also unnecessary: exact Joseph zero loss implies H e=0 at an applied correction. The literal nullspace stack can therefore be parameterized by the actual applied measurement Jacobian/service rows and continuous prediction/reset transports, avoiding the open full-P seed blocker. The exact-rational kernel lemma now has a service-row stratum interface and explicitly refuses callback-pattern enumeration. What remains is a theorem-domain cover of the continuous service-row/factor parameters, not a finite schedule list.

## Exact block-nullspace lemma formalized

Added complete_word_nullspace.py. For the authoritative literal zero-action row blocks C_proc/S, C_mag, C_acc, C_gyro/process and C_BA, the lemma is exact: if a supplied nonzero compatibility vector r satisfies Cr=0 and exact rational elimination gives rank(C)=n-1, then ker C=span(r). No SVD tolerance, independent observability lemma or four-S=>AW=0 shortcut is used. Regression includes a deliberate hidden-second-kernel case and fails closed. The remaining step is literal event-stratum assembly: export the actual complete-word homogeneous zero-action matrices into this stack and check rank/nullity on every reachable regular stratum. A* persistence exclusion remains a separate physical theorem applied only after the literal kernel equality.

## Local affine defect exporter in progress

The retrospective word defect is being replaced by literal boundary composition. The temporary proof header now snapshots xext/qref at every internal prediction, accepted acc/S/mag correction and reset; sync/sync-completion are classified mean-neutral. Consecutive mean-operation pre-snapshots therefore provide exact post/pre boundaries, including quaternion-reference changes across resets, at one unchanged physical sample time. For each boundary d_k=e_{k+1}-A_k e_k and b_{k+1}=A_k b_k+d_k. The endpoint residual e_N-M e_0 is reserved only for verification. Algebra and boundary-pairing regressions are committed; native execution/parity must pass before the locally composed b_W is used in linked chi/gamma.

## Linked entry ratio: correct optimization and scope

For one literal word the relevant feasibility value is inf_{0<gamma<lambda_min(Delta)} chi_gamma/gamma, not the historical choice gamma=lambda_min(Delta)/2. The carried linked diagnostic now performs the one-dimensional optimization on its exact same-word M,b,J0,JN. This is still retrospective because b=e_N-M e_0 is reconstructed from the realized endpoint; it is not a source-uniform endpoint-defect map valid under varying root error. No carried pass/fail is promoted until b is generated from the same physical SLOW+FAST/process/nonlinear source history. The theorem target remains sup_history inf_gamma chi_gamma/gamma < .0225 with gamma chosen measurably/continuously inside the same word's homogeneous loss cone, plus every-prefix retention.

## Outer-to-inner entry: exact status after qualitative strictness

The compact outer release region now exists at regular A21 roots: C_out=sup_Rrelease V<infinity. Qualitative finite-superword homogeneous strictness gives, on every compact annulus r_in^2<=V<=C_out, a positive homogeneous loss margin gamma_ann>0 after choosing a sufficiently long common superword. This is NOT yet finite entry under physical forcing. For the same reachable word the exact completed-square identity requires chi_gamma/gamma<r_in^2, with r_in<.15, and every prefix needs its own linked bound. The candidate SLOW+FAST source theorem proves the geometry that makes gamma_ann positive but does not numerically bound the linked chi_gamma supremum. Compactness gives only finite Chi_ann=sup chi_gamma, not Chi_ann/gamma_ann<.0225. Therefore outer retention/finite entry into the local ball remains open at one precise linked finite-error inequality; do not infer it from LaSalle strictness alone.

## Outer A21 release region: existence closed at regular roots

Captured-domain H18 release is compact in mean/tuner/covariance/scheduler/reference coordinates, and held-H18 LIN BIBO is source-uniform. After the first complete 16-s regular A21 LIN window, the recurring root covariance certificate gives a uniform Loewner lower P_root>=L_root>0, hence J_root<=L_root^{-1}. Therefore V=e'J e is continuous and uniformly bounded on the compact reachable release image at those roots. Define C_out=sup_Rrelease V<infinity and R_outer={reachable regular A21 roots with V<=C_out plus the compact carried auxiliary state}. This is an existence-level outer region; no numerical C_out is asserted. BA is eliminated only through the correct marginal P_oo when an outer coordinate reduction is used.

## Exact complete-word nullspace intersection: corrected formulation

Do NOT use the retracted statement that four S atoms force base AW=0. The authoritative q=0 classification is: complete-word zero action iff one literal homogeneous trajectory has zero fresh process/sync action, zero S action, magnetic-axis AG compatibility and zero accelerometer compatibility action at every event. Estimator-level forward-compatible ZG trajectories can satisfy these equations locally; four-S alone does not exclude them. The new ingredient is external to that algebra but same-history: the candidate SLOW+FAST + bounded-velocity theorem closes the real-arithmetic source transfer and proves the Corollary-A* aggregate nominal force cannot remain in the field-compatible zero-action geometry on every required moving window. Therefore a shipping-closed physical ZG execution persisting through a complete qualified moving superword is excluded. Combining this with MAGNETIC SERVICE leaves only the declared gauge/tilt-BA compatibility line on finite prefixes; persistence of that line through every qualified moving superword is excluded by the same A* physical source margin. This is the proper LaSalle invariant-set exclusion. It does not assert four-S pointwise AW nullity or the stronger numerical G0 floor.

## Coupled zero-action kernel condition after A* closure

A* alone does not prove the nuisance-projected VF-A modulus: root AW/BA columns can mimic accelerometer attitude rows in a Schur complement. The valid zero-action argument must impose all complete-word zero conditions simultaneously. In a zero-action limit, four-S plus LIN/AW process action forces the independent LIN/AW homogeneous component into its established null (zero on the regular four-S block); BA process/projection leaves only the held compatibility component; MAGNETIC SERVICE restricts AG to the transported field-compatible attitude/BG line; accelerometer zero action then sees the A* aggregate attitude block with no independent AW nuisance left to mimic it. Hence any surviving zero-action vector lies in the established field-axis tilt/BA compatibility line, not an additional normal direction. This is a qualitative kernel statement; it does not assert a numerical Schur modulus. The remaining proof obligation is to write this intersection as an exact block-nullspace lemma against the literal complete-word matrices and regress its rank on all event strata.

## Complete linked matrix revisited after A* source closure

The historical direct-information attempt failed at VF-A because no theorem then separated the estimator nominal accelerometer coefficient from the magnetic-compatible nuisance span. That premise is no longer absent: the new source-uniform real-arithmetic AW transfer proves Corollary A*'s nominal signed-mean separation under the candidate MARINE/SLOW+FAST qualification. This does not localize to same-cell magnetic epochs, but it does apply to the complete aggregate accelerometer block, exactly where the linked reduced-matrix/complete-word formulation needs it. Therefore the same-cell c,b0 route is retired and the complete-word zero-action nullspace argument is reopened with A* as an established aggregate premise. The next target is qualitative: any normalized zero-action sequence must have magnetic-compatible AG direction, zero four-S/LIN-AW action, zero accelerometer aggregate attitude action, and zero gyro/process action; A* excludes a nonzero attitude component transverse to the field, MAGNETIC SERVICE excludes heading/axial-BG complements, and four-S/process kills independent LIN/AW complements. The only allowed common null is the already identified physical tilt/BA compatibility line. A quantitative rate is still separate.

## Same-cell c,b0 route fails at c localization

The proposed exact E=0 same-prediction-cell acc->mag construction cannot be made source-uniform from the current premises. MAGNETIC SERVICE guarantees informative applied magnetic rows on every 1-s interval, but it does not require those rows to occur at accelerometer cells whose nominal force is noncollinear with the reset-pulled magnetic direction. The newly closed Corollary-A* premise is an aggregate signed-mean statement; it proves a positive attitude-column Gram over the window, not a pointwise/same-cell angle at the magnetic epochs. Therefore no uniform c>0 follows for the same-cell groups. Per the failure protocol this route stops here; do not infer c from average force separation. The valid continuation is the aggregate linked six-column reduced matrix, which retains all accelerometer rows, actual magnetic service rows and chronological gyro transport jointly.

## Real-arithmetic source-only AW transfer CLOSED for Corollary A* positivity

Under the user-qualified MARINE profile and candidate SLOW+FAST temporal caps, the conservative source-only charge is 1.07943484027807165 m/s^2 < g/5=1.96133, leaving .88189515972192835 m/s^2 before float32. The finite-error reset terms d^2 v/6+c(r)v^2 are downstream nonlinear supply, not source. The normalized quaternion polynomial is already included in the qualified gyro prediction transport angle/singular-floor certificate and is not charged a second time. Therefore the real-arithmetic source-uniform Corollary-A* signed-mean premise is CLOSED. This establishes positive attitude-column Gram pointwise; compactness of the qualified regular word class yields existence of a uniform positive attitude-Gram minimum. It does NOT promote the historical explicit G0 floor, whose stronger m_perp<=.4 and u1<=1.2 premises remain unproved.

## Corollary A* versus explicit G0 constants

Do not conflate the two thresholds. The recovered physical/source budget targets Corollary A* positivity: m_perp<g/5=1.96133 m/s^2. This is enough for a positive attitude-column Gram and, with compactness plus magnetic service, qualitative strict geometry. The published explicit Theorem G0 floor 1.486786e-3 additionally assumes the much stronger numerical window premises m_perp<=.4 m/s^2 and nominal force L1<=1.2. The carried worst (.348,1.091) satisfies them diagnostically, but the new .895538 source allowance does not prove them source-uniformly. Therefore do not promote the old explicit G0 constant from the new Corollary-A* bridge. Use the qualitative positive-Gram route for source-uniform existence; a quantitative contraction constant still needs stronger window-statistic enclosure or a recomputed G0 at the proved bounds.

## Jump/twist source remainder resolved structurally

The exact reset remainder is d^2 v/6 + c(r)v^2 plus the literal quaternion polynomial defect. The d^2 v and v^2 pieces depend on homogeneous finite error and therefore are NOT part of the source-only physical-to-nominal AW transfer; they belong to the existing linked nonlinear/prefix-retention supply. The linear correction jumps remain in the signed ordered injection product. At zero homogeneous error the only reset source defect is the literal small-angle quaternion polynomial/arithmetic term; the polynomial defect is 4(d^6/46080+d^7/645120) for d<.01 and zero on the exact-real trig branch. Float32 remains a separate implementation obligation. Thus the previously named source-uniform correction-jump/twist obstruction is removed from the real-arithmetic source-only AW transfer rather than bounded by an invented per-update NIS ceiling.

## Candidate FAST gyro removes the historical continuous field-axis obstruction

The old continuous obstruction used n_g=.02 sin(t)b and required about .04 rad signed accumulation over a half-cycle. The candidate C_g=.002 rad excludes it by a factor 20. In the exact 17-s field-axis Stieltjes telescope, the inherited slow charge is 2 P_max B_g,s/T + P_max D_g,s = .0129961764705882 m/s^2 and the candidate FAST primitive contributes at most P_max C_g/T = .000647058823529412 m/s^2, total .0136432352941176 m/s^2 before coefficient-rotation and exact SO(3) correction-jump/twist terms. Thus the historical .055 m/s^2 witness is no longer admissible. The controlling correlated remainder is now the source-uniform sum of coefficient rotation plus exact correction jump/twist action; do not fall back to independent gain or total-NIS norms.

## Signed AW reader: uniformity status

The carried six-profile stress family now passes the recovered AWB2 allowance: worst signed 16-s error .371143 m/s^2 versus .89553839501604595, leaving .52439539501604595 m/s^2. Compactness/continuity of the qualified regular 16-s same-history word class proves existence of a finite source-uniform B_AW=sup|Phi|. Strict numerical closure is reduced to one dependency-preserving finite-cover inequality L_Phi delta_F + E_nl + E_f32 < .52439539501604595. L_Phi and delta_F are not yet certified; independent gain/tuner/covariance boxes are invalid because they destroy the paired cancellation. This is now the controlling quantitative AW-transfer obligation.

## Recovered non-recurrent-acceleration bridge

The earlier useful mechanism is recovered explicitly. It does not assume recurring acceleration. Bounded physical velocity gives the signed physical mean charge 2 V_max/L, jerk/sample fidelity adds J_max h_acc/4, and Corollary A* needs only a nominal signed AW mean transverse to the magnetic field. At L=16 s the physical charge is 0.8375 m/s^2. Under the candidate H_a=60 s,C_a=.05 m/s profile, the FAST signed mean over 16 s is <=.05/16=.003125 m/s^2. Even conservatively charging the entire slow accel amplitude .22516660498395405 leaves 1.96133-.8375-.003125-.22516660498395405 = .89553839501604595 m/s^2 of field-separation margin at the sensor-history level. This is source-uniform and uses no recurrent-acceleration premise. The remaining bridge is specifically the literal time-varying acc/S gain loop: sensor signed mean does not automatically bound nominal AW signed mean because gain rectification is possible. The exact AW-loop/two-Abel identity must transfer this .895538395 margin without replacing the gain-weighted history by independent boxes.

## Candidate assembled-device FAST profile

Use H_a=H_g=60 s, C_a=.05 m/s, C_g=.002 rad as the current candidate qualification. On any 30-s MARINE window the signed FAST accumulation caps remain .05 m/s and .002 rad. The known phi=.01 sin(t/2) witness needs at least .370655073155 m/s and .0197081716368 rad after sampling charge, so this candidate excludes that witness by either sensor separately, with reserves .320655073155 m/s and .0177081716368 rad. This is not yet the general same-history gauge theorem; hardware validation of these candidate caps also remains required.

## User-qualified MARINE excitation constants (2026-10-01)

For the controlling theorem, use T_E=T_P=30 s, theta_E=2 deg, P_E=.03 m. These are theorem qualification premises supplied by the user; they are not inferred from RAO data. IMU fast H_a,C_a,H_g,C_g remain OPEN. Consequences: P_E/T_P=.001 m/s; theta_E/T_E=.00116355283466 rad/s. The generic signed acceleration chord lower max(0,P_E-T_P V_max) remains zero, so no pointwise acceleration floor follows.

# OU-III proof: SLOW + FAST physical IMU qualification

Base main: `eef30627130f434eb14ea9f42f4df140a291d341` (2026-10-01).
Controlling derivation: [two-timescale IMU model](ou3-imu-two-timescale.md).
The previous ledger is preserved **verbatim** in
[the pre-migration history](ou3-proof-research-state-before-slow-fast.md).
Its suggested LF cutoff, 1-degree/60-second budget and stronger excitation
are research proposals, not established device qualifications or premises.

## Current hypothesis

Keep the single path construction -> capture -> magnetically informed H18 ->
refinement/release -> recurring A21 -> regional practical stability. Keep the
literal coupled (tau, sigma_aw, R_S, T_S), prediction/due-S/acc chronology,
full-covariance dissipativity and same-history LaSalle structure. No shipping
code, calibration behavior, MARINE MOTION, MAGNETIC SERVICE or gate is changed.

Physical calibrated errors are e_a=b_a,s+b_a,f and e_g=b_g,s+b_g,f.
Slow components have norm/rate bounds and linked increments min(2B_s,D_s*h).
Fast components have amplitude bounds plus signed accumulation on EVERY placed
window 0<T<=H: |integral b_f^hold|<=min(B_f*T,C), C<B_f*H. This is a simple
conditional qualification interface, not an empirical claim of cancellation.
Both H,C pairs are OPEN. The six inherited numerical amplitude/rate values are
explicit candidate deployment budgets; calibration-fit statistics, process
noise, typical sensor densities and estimator projections do not prove them.
One decomposition must satisfy the complete history, including all transitions.

## Evidence and derived inequalities

`imu_temporal.py` implements conditional same-history calculus. Distinct raw
sample epochs on complete hold cells have the global error-class supremum

    Delta_i(h)=min(2B_i,s,D_i,s*h)+c_i(dt_0)+c_i(dt_1),
    c_i(dt)=min(B_i,f,C_i/min(H_i,dt)).

The raw fast term can still be 2B_f. Averaging/cancellation never justifies a
smaller raw-point bound without the cell qualification. With K*(T) the tiled
window cap, shifted boxcar differences have fast bound
2 min(K*(L),K*(h))/L. Matrix Abel identities retain signed gain/frame variation,
the endpoint primitive and division by sample duration. They are support
bounds on one history, not a reason to independently maximize word loss and
forcing. Unknown temporal profiles fail closed in the physical validator and
in the sensor-supply helper. A finite prefix never certifies an all-time device.

The exact linked completed square remains valid:

    G=J0-M' JN M-gamma J0>0, z=M' JN b,
    chi=b' JN b+z' G^-1 z,
    V_N=(1-gamma)V_0+chi-(e_0-G^-1 z)'G(e_0-G^-1 z).

The required supremum is over the SAME reachable slow/fast histories that
produce M, b, covariance, tuner and scheduler. Physical BA/BG coordinates now
refer to the slow components; fast errors stay in the sensor forcing. The
H18 and active-A21 bias prediction identities keep their literal factors.
Every-prefix retention and chi/gamma<r_in^2 remain unproved.

## Counterexamples: classify by model, do not erase

The old V3 norm-only witness phi=.01 sin(t/2), with arbitrary n_a/n_g and a
constant BA offset, remains valid under its explicitly frozen OLD model.
Its all-slow realization requires .04903325 m/s^3 and .0025 rad/s^2: respectively
49.03325 and 250 times the inherited candidate slow-rate budgets. This excludes
only the all-slow realization. For ANY split, opposing half-windows of length
2*pi require continuous K_a>=.3725224 m/s and K_g>=.0198026 rad. Charging the
qualified 6-ms sampling defect gives conservative necessary requirements
K_a>=.3706550 m/s and K_g>=.0197081 rad. These are witness requirements, NOT
chosen sensor specifications. Missing H,C means current-model admission and
exclusion are both OPEN. The old sqrt(V)>=.4 bound cannot be transferred by
changing which part of the error is called physical slow BA.

The former phi=.001 sin(t/40)^3 construction has p=v=a=0 and therefore does NOT satisfy the amended MOVING displacement-span condition. It is retained only as a counterexample to the previous attitude-only MARINE contract. Charge the tiny gravity-representation
constant to slow accelerometer bias. Its slow amplitude/rate bounds obey the
inherited candidates, and the quiet execution retains MAGNETIC SERVICE.
It fits symbolic MARINE T_E=80*pi, theta_E=.002 on complete moving windows.
Thus that zero-translation construction no longer blocks the amended MARINE contract. Full moving compatibility must now include the nonzero displacement span together with attitude, SLOW+FAST IMU and magnetic service. This does
not refute every fixed deployment pair (T_E,theta_E), nor practical stability
relative to the compatible class. No EXCITED_MOVING strengthening is adopted.

## Retained facts and failed approaches

Retain the operation energy identities, linked completed square, signed
variation of constants, BA marginal elimination, source covariance comparisons
within their stated profile, and conditional radius-local field-axis theorem.
The previous .15 local result still requires its ACTUAL local-error premise.
Do not upgrade it to outer entry from physical tilt span alone.

**DEAD ENDS:** unrestricted residual boxes labelled fast; deriving deterministic
all-history caps from RMS/Allan/typical noise alone; fitting a cutoff to reject a
witness; deleting gyro or frame-weight variation; using only aligned-window
means; resetting the split at proof boundaries; inferring zero base innovations
from zero homogeneous action; S-only LIN contraction substituted for actual
interleaved H18; lower covariance bound substituted for full compactness; old
O1/O2 kernel-ceiling architecture used as the controlling route. Archived algebra
may remain correct as algebra or an explicitly amplitude-relaxed outer bound.

## Exact joint gauge qualification

The controlling quantifiers are now written explicitly in [joint SLOW+FAST physical gauge separation](ou3-joint-gauge-separation.md). For fixed qualified MARINE tuple m=(T_E,theta_E,T_P,P_E), Q_IMU(m) is the set of fast horizon/cap tuples for which the infimum of the SAME-history stacked acc+gyro causal compatibility residual, restricted by actual MAGNETIC SERVICE, is strictly positive. Process/accelerometer columns are paired and physical p/v/S/a boundaries telescope before norms. Because the four MARINE excitation constants are also numerically OPEN on main, the current empirical qualification is Q_phys over both MARINE and IMU parameters; no four-dimensional numerical Q_IMU is asserted before m is fixed.

## Current limiter and alternatives

Numerical assembled-device slow/fast qualification is OPEN for both sensors.
The complete physical joint acc/gyro/magnetic reachable ambiguity must then be
compared with the EXISTING fixed MARINE constants, not a bias-only threshold.
Held-H18 LIN BIBO and captured-domain release compactness are CLOSED qualitatively/source-uniformly by the post-#640 certificate; outer retention,
physical-to-nominal transfer, linked finite supply, nonlinear/prefix retention,
regime composition and float32 totality also remain OPEN. All end-to-end theorem
flags remain false. No finite carried replay supplies the missing uniformity.

An observable-consistency-class practical target remains a legitimate option,
but is not substituted for the current theorem without stating the difference.
No behavior change or automatic STILL detector is justified here.

## Next falsifiable calculation

Independently qualify both delivered-stream fast horizon/cap pairs and the
slow budgets with residual/reference evidence. For a declared (T_E,theta_E),
intersect BOTH exact sensor compatibility equations, their transported gyro
integral, the physical kinematic continuation and actual MAGNETIC SERVICE.
Evaluate the linked chi and every-prefix bounds on that same reachable class.
Do not spend more enclosure effort on the refuted old norm-only .15 claim or
try to close a source theorem by multiplying independent norm suprema.

## Validation scope

See `imu-two-timescale-certificate.json` for exact rational scalar consequences
and `test_ou3_imu_temporal.py` for finite-window/cross-boundary algebra tests.
Finite-prefix and synthetic-test success is not device qualification, all-time
float32 validation, a uniform contraction certificate or theorem completion.

## CI: information-ratio replay comparison

Failed quantity: `information_ratio_source_diagnostic --expect` on the
`theorem` job compared `diameter.relative_difference` (committed 1.03e-5 max,
fresh 2.3e-6 on the same case) with a 5e-3 relative tolerance. Classification:
implementation/CI. That leaf is the residual of an identity between two equal
quantities, so its value is float32 replay noise and is not reproducible
across toolchains. Invalidated hypothesis: every recorded leaf is well posed
for relative comparison. Retained facts: the identity contract
`relative_difference < 1e-3` is still enforced by `verify_diagnostic`; the
replay comparison now only requires both records below 1e-4; every other leaf
keeps the 5e-3 relative comparison. Current limiter and next falsifiable
experiment are unchanged from the sections above.

## CI: held-BA LIN word diagnostic did not execute

Failed quantity: `held_ba_lin_word_diagnostic` on the `theorem` job raised
`JSONDecodeError` (`"root_covariance":,`). Classification: implementation.
Its driver patch removed the root-covariance capture, compiled against the
untapped shipping header (zero recorded events), and recorded past first BA
activation, so the "last 3400 predictions" would not have been the held
window. Invalidated hypothesis: none; the diagnostic had not produced a
number. Fix: capture the root at the first Live sample, compile with the
`ag_readout_source_diagnostic.instrument` tapped header copy, and stop
recording permanently at the first `acc_bias_updates_enabled` transition.
Retained facts (finite carried evidence only): quiet live 18051 / active
24064, rho_LIN=0.03953, sigma_max=0.3258; wave live 6368 / active 24016,
rho_LIN=0.006467, sigma_max=0.03199. This does not establish source-uniform
held-BA LIN BIBO stability; the next falsifiable experiment above stands.

## CI: captured-domain H18 refinement margin was negative at 7 degrees

Failed quantity: `refinement_gate_margins` horizontal-fraction margin at the
captured tilt bound, `test_ou3_h18_release` on `ou-evidence / commit`. With
|B| in [20,75] uT, B_h>=15 uT, residual <=2 uT and tilt 7 degrees the rotation
loss is 2*75*sin(3.5 deg)=9.157 uT, so the horizontal fraction lower bound is
3.843/77=0.049905 < 0.05 (margin -9.46e-5). Classification: mathematical
(premise constant). Invalidated hypothesis: the 7-degree captured domain gives
positive margin at the literal 5% horizontal-fraction gate. The largest
admissible tilt is `refinement_tilt_limit_rad` = 6.9944 degrees. Retained
facts: on the <=6.9-degree domain the margin is +1.60e-3 and the norm-ratio
margin +0.1278; the conservative captured-domain release bound is unchanged
(469 s); compact A21 release remains conditional on capture, now into the
<=6.9-degree domain. Current limiter and next falsifiable experiment are
unchanged; general capture must now reach the smaller domain.

## CI: source readout correction events lacked the physical lift

Failed quantity: `ag_readout_source_diagnostic` on the `theorem` job raised
`event 2 missing physical lift fields` in `validate_physical_lift`.
Classification: implementation. The tapped correction observer captured
`estimator_state`/`estimator_quaternion` before `xext += K r`, but
`readout_correction` in `tools/stability/ag_readout_source.cpp` did not append
`physical_json()`, so correction boundaries carried neither the replay truth
nor the pre-correction estimator snapshot required by
`carried_local_defect_certificate`. Invalidated hypothesis: none; the
diagnostic had not produced a number. Fix: correction events now export the
same physical/estimator snapshot as prediction, sync and reset events.
Retained facts: terminal parity between observed and untapped shipping output
is still enforced, and the local forcing endpoint-parity gate is unchanged.
The committed signed-adjoint and signed-balance diagnostics fingerprint this
driver and had not been regenerated after the earlier observer snapshot
changes; they were regenerated. Their charges moved by about 1e-8 relative,
and the signed-balance margins stay negative (S-interval -19.03910780414,
rotation triangle -30.45431596760). Current limiter and next falsifiable
experiment are unchanged.
