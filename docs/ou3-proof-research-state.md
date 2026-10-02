## FAST signed-primitive tightening for guard: adjacent-sample route does not improve the bound

The requested signed-primitive calculation was carried out for the vibration-guard one-pole kernel. With piecewise-held FAST input u_k, U_n=h sum_{i<=n} u_i and every placed primitive bounded by C_a over H_a, discrete summation by parts gives a finite-window filtered-output bound of the form
|LP(u)| <= (1-alpha) B_f + (1-alpha) C_a/h + alpha^(H_a/h) B_f
(componentwise; later positive unit-gain LP poles do not increase the L-infinity bound). Therefore
|LP(u)-u| <= B_f + that quantity.
The tail is negligible for H_a=60 s.

Numerically at h=.005 s and fc=3 Hz, alpha=exp(-2*pi*3*.005)=0.91005724. The primitive-derived unweighted first-pole term is about 0.9264 m/s^2 before adding the raw B_f in LP-minus-raw form, whereas the simple pointwise difference bound 2 B_f is 0.6 m/s^2. Thus C_a=.05 m/s does NOT tighten the adjacent-sample LP/raw difference at this sampling/corner; its C_a/h conversion is too expensive. Replacing 2B_f by this result would make the certificate worse.

guard_fast_primitive.py records the derivation and selects min(pointwise,primitive), then applies the actual same-history guard weight w. This identifies the correct next use of the signed primitive: not as a surrogate adjacent jump bound Delta_k, but directly in the weighted convolution w_k(LP-I)u over the guard slew horizon, where w changes on 5 s and the Abel kernel can exploit C_a without paying C_a/h at every sample. No new assumption is needed.

This is a useful negative calculation: the previously proposed adjacent FAST primitive tightening is analytically unproductive and should not be pursued. Keep the valid 2B_f term in the raw-increment recurrence for now; attack the weighted FAST correction as one temporal operator if CI shows it is the limiter.

## Guard LP/raw correlation recurrence connected to physical history

The vibration-guard norm proof now uses the literal convex-blend structure instead of the decorrelated component perturbation. For each low-pass pole define d_p,k=stage_p,k-f_k. Shipping initialization gives d_p,0=0 exactly. With alpha=exp(-2*pi*f_c*dt) and Delta_k=||f_k-f_{k-1}||,

||d_0,k|| <= alpha (||d_0,k-1||+Delta_k),

||d_p,k|| <= (1-alpha)||d_p-1,k|| + alpha(||d_p,k-1||+Delta_k), p>0.

Since f_c,k=f_k+w_k d_LP,k and 0<=w<=1,

||f_c,k|| >= ||f_k||-w_hi ||d_LP,k||.

The proof guard carries this scalar coupled recurrence separately from the component interval LP state. Thus the huge XYZ diameter no longer directly destroys the Mahony normalization floor.

The same-history adjacent raw increment is connected from existing physical/source assumptions. The current conservative bound is
Delta_k <= (J_max + g Omega_max + A_max Omega_max + D_a,s) dt + 2 B_a,f.
The first four terms respectively cover translational jerk, rotating gravity, rotating translational acceleration and slow accelerometer-bias drift. The final FAST term is a conservative adjacent delivered-value jump; it preserves the existing FAST amplitude qualification but does not yet exploit the much stronger signed primitive cap C_a=.05 m/s. Therefore this is rigorous but potentially loose. If it limits the recurrence, the next tightening should derive adjacent/LP-convolved FAST contribution directly from the carried signed primitive history, not add a new physical assumption.

Tests pin d_p,0=0 and the weighted LP-difference floor. No estimator behavior changed.

## First executed source-uniform release failure: sample-0 Mahony norm dependency; coupled norm repair

CI run 37039885269 is the first run that reached the new interval-release diagnostic. It failed at history-leaf sample 0 with "normalization denominator not separated from zero". This is classification D, not a shipping counterexample. The cause is dependency loss: component intervals for the accelerometer contain the zero vector although the physical specific-force norm cannot vanish under the declared marine/bias envelopes.

The leaf also exposed a modeling omission: Mahony must receive delivered specific force including gravity in body coordinates, not translational physical_a+bias. The proof leaf now constructs a conservative body component hull g*gravity_dir + R^T a_phys + slow + fast, using |R^T a_phys| component <= A_max. Crucially it retains the coupled reverse-triangle norm floor
g-A_max-B_a_s-B_a_f = 9.80665-8.8-0.22516660498395405-0.3 = 0.4814833950160458 m/s^2.
MahonyBox accepts this coupled norm floor when normalizing, preventing the component hull from manufacturing a false zero denominator.

Because AccelVibrationGuard conditions the vector before Mahony, the raw norm floor is not reused blindly. causal_vibration_guard_interval now exports an upper bound on ||f_conditioned-f_raw||, and the leaf supplies max-information floor f_raw_lower-delta_guard to Mahony. If that quantity becomes nonpositive the enclosure again fails closed at the actual sample; no floor is invented.

This repair does not assume recurrent acceleration, absolute yaw, or an independent Mahony input. It is a rigorous relaxation of body translational direction using the existing A_max norm plus carried gravity-direction/bias histories. It may be wide; any later failure is a D-level enclosure-width issue unless an admitted execution is constructed.

## Interval closed-loop mean state implemented; full-horizon execution seam reached

closed_loop_release_stream.py now carries the release mean as (x_mid,x_rad). Literal OU prediction hulls phi_va/phi_pa/phi_Sa/alpha over the same generated tau interval. The same AW cell constructs the world-frame H_acc interval; the same physical-acceleration/conditioned-sensor/AW/BA cells construct r_acc; covariance-derived interval K is immediately multiplied by that residual cell; x is updated by the resulting interval K r; and the same dtheta interval constructs G=I+0.5[dtheta]x for covariance congruence. S correction follows the same pattern. Each accelerometer event records dominant radius attribution among AW, BA, residual, K, reset and guard.

A branch regression in literal_history_leaf_propagator was found before running this: executable text contained a literal backslash-n sequence and the code still used an independent A_max acceleration cube despite the prior intended repair. Both are corrected on the current branch; the leaf again uses physical_acceleration_history(cell,...) derived from velocity knots plus jerk and exports the same FAST/guard/Racc streams.

The interval driver fails closed when the tau-dependent S period is not point-resolved, innovation inversion fails, or a mean radius becomes nonfinite. It returns first_failure_sample/reason and radius attribution instead of independently widening generated outputs. Magnetic callbacks remain aggregate.

Execution limitation on this turn: repository tools available here can edit/read GitHub and inspect Actions but cannot execute the Python theorem modules or dispatch a new arbitrary workflow. Therefore no fabricated first-failure sample is recorded. The code path and tests are committed; the next actual CI/native execution must run generate_interval_stream on a source-uniform root subdivision and report the first returned failure. If the unsplit root fails immediately from S-period width, split the causal physical history coordinate identified by tuner sensitivity; do not split tau itself.

Structures preserved: literal OU small-x coefficients, same-history tau/sigma/R_S/guard state, velocity-derived physical acceleration, conditioned accelerometer, Racc, covariance-derived K, joint K r, reset from same dtheta, persistent S scheduler, aggregate magnetic service. Relaxations: coefficient hull over tau endpoints assumes endpoint extrema for the OU coefficients and requires a monotonicity audit before final promotion; if not proved, replace with derivative/interval evaluation. Failure classification remains D for enclosure blow-up, not shipping instability, unless an admitted same-history counterexample is constructed.

Next calculation after execution: feed each verified prediction/S/acc/reset event directly into structured gain export, vector AW recurrence and covariance-defect recurrence to obtain the actual (A_k,a_k,b_k,c_k,D_k) prefix sequence. The architectural seam is now implemented.

## Causal AccelVibrationGuard interval and joint interval K r implemented

causal_vibration_guard_interval.py now mirrors the literal guard state sample-by-sample: conditioning low-pass cascade, two detector high-pass stages, removed-ms EMA, RMS/excess calculation, engagement target/clamp, slew, rail parking and conditioned output. Generated guard state is persistent in the HistoryCell leaf. Clamp/rail crossings are enclosed by interval hulls; no midpoint branch is selected. The same conditioned accelerometer is fed to the causal adaptation path, and the same excess-RMS interval generates Racc through sigma_eff^2=(sigma_base*scale)^2+(gain*excess)^2. literal_history_leaf_propagator now exports conditioned_accel and Racc_interval from that same history.

interval_kr.py implements K*r directly for interval K and interval residual, including the bilinear radius term. closed_loop_release_stream exposes interval_correction_product(P,H,R,r_mid,r_rad), which derives K from the current covariance/measurement cell and immediately multiplies the same residual cell. It exports dtheta from that same product. This removes the previous point-leaf-only K*r blocker and does not replace it by ||K|| ||r||.

Important remaining issue before claiming the full release iteration: the closed-loop generator still needs its 21-state mean represented as an interval vector and its accelerometer H generated from the same AW interval. Point generate_point_stream remains diagnostic only. The new guard/Racc and K*r primitives are interval-capable, but the sample loop has not yet been converted to update (x_mid,x_rad,P) after each correction/reset. Also the current root's 30-s physical knots can make the guard interval broad; this is a cover-width issue, not a new physical assumption, and should be handled by source-history subdivision/sensitivity rather than boxing generated guard state independently.

Structures preserved: literal guard equations and persistent state; conditioned sample shared with adaptation; exact Racc algebra; covariance-derived interval K; same-cell bilinear K*r; dtheta from K*r; no independent generated-output boxes. Relaxations: ordinary component interval arithmetic in the guard loses cross-axis/state correlation but is a rigorous causal outer; branch hulls may be wide. No shipping counterexample.

Next calculation: intervalize the closed-loop mean state itself. Carry x=(mid,rad), build H_acc from the same AW cell, build r_acc=a-aw+sensor-ba from the same source/mean cells, apply interval_correction_product, update x, construct the interval reset G from dtheta, and update P. Then feed every resulting operation into the existing A/a/b/c/D recurrence across the release horizon. If innovation or reset cells blow up, report the first time and dominant radius instead of widening generated outputs independently.

## Closed-loop release operation-stream generator implemented; Racc state is the next exact seam

closed_loop_release_stream.py now starts from the explicit goLive covariance/zero mean and advances the literal [v,p,S,a_w] OU mean using the shipping phi_va, small-x/exact phi_pa, phi_Sa and alpha branches, plus BA OU decay. It computes covariance-derived K from the same P/H/R boundary, applies dx=K r to the same mean, then applies the literal first-order attitude covariance reset and zeroes the local injected theta coordinates. S due events use the progress-preserving scheduler and the same applied tau. The accelerometer world row uses the attitude-free Lemma-W representation with AW and active BA columns. Magnetics remain aggregate; no callback pattern is invented.

The stream is deliberately point-leaf/fail-closed at this stage: any non-point generated tau/sigma or physical/source interval must be split/enclosed before promotion; midpoint substitution is rejected. literal_history_leaf_propagator now exports the linked physical_acceleration and FAST-acceleration streams needed by the generator.

A shipping term caught during implementation is vibration-conditioned Racc. The actual wrapper commands sigma_eff=sqrt((sigma_base*scale)^2+(gain*excess)^2) from AccelVibrationGuard state. The closed-loop stream now refuses to silently use nominal 0.2^2 I and requires a same-history Racc stream. causal_racc.py exposes the exact algebra if the guard excess-RMS state is supplied and the already-qualified std<=0.3010398644698074 ceiling as a rigorous fallback. The fallback is source-uniform but loses residual/gain correlation, so it is not automatically promoted as the desired final linked calculation.

Thus the requested closed-loop architecture now exists: seed -> prediction -> due S correction/reset -> accelerometer correction/reset -> next sample, with K and r generated at the same boundary. The final source-uniform numerical run is currently stopped specifically by the absence of a causal interval replica/export of AccelVibrationGuard's detector/excess state (and, for non-point root cells, interval correction products rather than point-leaf splitting). This is narrower than the previous mean/covariance-operation-stream gap.

Structures preserved: literal OU mean coefficients including small-x branch; explicit goLive P0; covariance-derived K; K r mean correction; first-order covariance reset; progress-preserving S scheduler; active BA accelerometer column; attitude-free world row; aggregate magnetic service; linked physical acceleration and SLOW+FAST source. Relaxations: causal_racc qualified-ceiling fallback is marked unlinked and non-promoting for the final shaped result. No shipping counterexample.

Next calculation: intervalize the closed-loop K r operation on one HistoryCell and propagate the AccelVibrationGuard detector/excess state causally (or prove that using the qualified Racc ceiling is monotone/conservative for the particular A,b,c,D certificate despite lost correlation). Then execute the full release horizon and feed each operation directly into structured gain/AW/covariance-defect recurrence.

## Same-history geometry/residual propagator: algebra closed; physical acceleration linkage repaired

same_history_geometry_residual.py now constructs, on one operation boundary, the literal accelerometer residual, Jacobians, local affine defect d_acc=r-H e, S residual r_S=-S_hat=e_S-S_phys, and reset injection dtheta=(K r)_theta from the SAME gain/residual product. Tests pin zero-error residual, S identity, and injection identity.

Absolute attitude is deliberately not added as a new independent HistoryCell coordinate. Existing Lemma W permits the accelerometer attitude geometry to be represented in world coordinates by f_world=a_w_hat-g e_z and [f_world]x, eliminating absolute yaw from the source-uniform row. The new world_boundary path uses the signed source physical_a-a_w_hat+sensor_error-(same-boundary BA/lever transports). This is the admissible route for the current root, whose gravity_dir knots cannot determine yaw.

While wiring this path, a more basic source-link defect was found in literal_history_leaf_propagator: physical acceleration had been replaced independently at every sample by the full [-A_max,A_max]^3 cube. That destroys the signed p'=v, v'=a/jerk linkage and makes any r_acc enclosure too loose by construction. history_interpolation.py now derives a same-history acceleration outer from the physical-velocity knots and J_max: for each signed secant over h,
|a(t)-(v(t+h)-v(t))/h| <= J_max |h|/2,
then intersects all available secants with A_max. The causal leaf now feeds this linked acceleration into the delivered accelerometer/tuner path instead of the independent A_max cube.

This completes the algebraic geometry/residual operation map, but not yet a promoted full release replay. The remaining source-uniform numerical gap is estimator-mean propagation inside the HistoryCell leaf: a_w_hat, S_hat, BA mean and structured covariance/gain must be advanced together so the world residual source, K r correction, dtheta reset and next covariance coefficients are produced at the same boundary. The existing golive_to_a21_release driver propagates an affine mean interval from an operation_stream, but that operation_stream is exactly what must now be generated causally rather than supplied externally.

Structures preserved: world-frame attitude-free row factorization; same-boundary local defect; literal K r reset; S identity; physical velocity/acceleration/jerk linkage; SLOW+FAST causal tuner path. Relaxations: the velocity-knot acceleration formula is a rigorous outer secant/Jerk enclosure and may be loose between 30-s knots, but does not introduce independent acceleration controls. No shipping counterexample.

Next calculation: implement the closed-loop estimator mean/covariance operation-stream generator starting from golive_release_seed, using each sampled physical/bias history row, causal tuner tuple and literal S due branch. It should emit prediction -> due S -> acc -> reset in shipping order (plus aggregate magnetic action separately), updating a_w_hat/S_hat/BA and structured covariance before the next sample. That stream is the final missing input to the (A_k,a_k,b_k,c_k,D_k) iterator.

## Source-uniform release seam connected to explicit seed and literal S scheduler

The explicit goLive construction seed is now connected to the causal history leaf instead of being reported as missing: zero 21-state mean, diagonal covariance upper seed with AW variance 16.48 and BA variance 0.004^2, zero cross covariance at handoff, plus the existing BA-graph numerical interval. The source-uniform release chronology also carries the shipping progress-preserving S scheduler. Its period is generated causally from the same applied tau via T_S=clamp(c_T tau,.005,.15); if a tau interval straddles a due/not-due boundary, the leaf fails closed and must split a shared history coordinate rather than choosing a scheduler branch.

Important MAGNETIC SERVICE correction: the existing proof already establishes that there is no finite callback-pattern cover under the aggregate service premise. Therefore source_uniform_release_chronology.py deliberately does not enumerate magnetic callbacks. It retains AggregateMagneticService(T_M,mu_M) and requires the later homogeneous/source certificate to consume that applied-information operator. Inventing a periodic magnetic callback schedule would be a shipping-invalid strengthening.

The old literal_history_leaf_propagator failure reason is corrected. P0/mean/BA-graph seed is available. The current missing connection is same-history attitude/measurement geometry and residual/local-defect stream: R_wb/force geometry, accepted accelerometer row, d_acc, physical/estimator S, structured covariance coefficients, and reset dtheta=Ktheta r at the same correction boundary. Until those are propagated from each HistoryCell, the requested full numerical (A_k,a_k,b_k,c_k,D_k) iteration remains fail-closed rather than populated with independent extrema.

Structures preserved: explicit shipping goLive seed; persistent scheduler progress; causal coupled tau->T_S; aggregate applied magnetic service; same-history requirement for residual/reset inputs. Relaxations introduced: none in the chronology seam. No shipping counterexample. Next calculation: construct the missing geometry/residual stream from the physical history interpolation plus estimator mean/covariance propagation, then feed structured gain, vector AW and covariance-defect iterators sample by sample.

## Structured covariance defect recurrence added after reset audit

To continue the requested release iteration without dropping the literal MEKF reset, write P_k=P_struct,k+E_k and carry D_k>=||E_k||_2. The exact reset split supplies the new source-linked non-INJ remainder. Between resets the remainder is propagated in literal order. Prediction gives
D^- <= ||F||^2 D + D_axis + D_Q.
For a Riccati correction C(P)=P-PH'(HPH'+R)^-1HP, with structured innovation S0 and q=||S0^-1|| ||H||^2 D<1, Banach/resolvent expansion gives an explicit finite D^+ bound; if q>=1 the sufficient envelope fails closed and is classification D, not instability. Reset gives
D^+ <= ||G||^2 D^- + D_reset,
where D_reset is obtained from the same-history dtheta=K_theta r split, not an independent reset-angle box.

structured_covariance_defect.py implements this chronology. This makes the intended iteration conceptually (A_k, structured a/b/c coefficients, D_k), with A_k updated by ordered vector K r products and D_k carrying both changing-axis and reset-generated non-INJ covariance remainder. The next numerical promotion step is to feed this iterator from the certified release leaf: literal structured P coefficients, H/S inverses, same-boundary residual vectors/local defects, accepted S/mag scheduler events, and dtheta from each accepted correction. The current source-uniform history propagator still fails closed before those release P0/mean/BA-graph inputs are connected, so no source-uniform numeric radius is claimed yet.

Structures preserved: literal reset congruence, correction order, covariance-derived gains, same-history dtheta=Ktheta r, accepted S/mag events, and coupled tuner history. Relaxation introduced: D uses spectral-norm perturbation/resolvent bounds around the structured covariance. A miss of q<1 is D-level sufficient-bound failure. No shipping counterexample.

## Literal gain export corrected; reset creates a fourth structured covariance component

The requested gain export exposed two shipping-faithfulness corrections. First, after accelerometer coupling the S pseudo gain is not generally scalar*I: the exact structured expression is
K_aw,S = P_aw,S (P_SS + R_S)^-1,
with all three factors carried in the INJ algebra. Second, when the shipping lever-arm Jacobian is enabled, the accelerometer AW numerator also contains P_aw,bg J_bg'. The exporter now retains both terms and also exposes the accepted-magnetometer AW gain P_aw,theta J_mag' S_mag^-1. No entrywise K box is introduced.

The source link is now explicit at the same pre-correction boundary:
r_acc = H_acc e + d_acc,
where d_acc is the existing local affine measurement defect, and
r_S = -S_hat = e_S - S_phys.
Thus the existing local-defect bound constrains d_acc, not r_acc independently; the state/source correlation must be retained. The AW mean iterator now applies ordered vector products K r for accepted acc/S/mag corrections before taking ||a_w||. The previous scalar sum of gain-norm times residual-norm remains only an explicitly labelled triangle relaxation.

Attempting the requested scalar-coefficient iteration through the literal release chronology found a missing reset term. Shipping applyQuaternionCorrectionFromErrorState applies the covariance reset G=I+0.5[dtheta]x to the attitude block/cross-covariances. This is not a common orthogonal rotation. For an AW-theta INJ block B=aI+b nn'+c[n]x and dtheta=d_parallel n+d_perp,

B G' = B (I - 0.5 d_parallel [n]x) - 0.5 B [d_perp]x.

The first term is exactly INJ and updates (a,b,c) algebraically; the second is a genuine non-INJ matrix defect. structured_reset_split.py carries this residual exactly and verifies the bound
||E_reset||_2 <= 0.5 ||B||_2 ||d_perp||.
Therefore an exact shipping release iteration cannot consist only of (A_k,b_k,c_k,D_axis,k) unless D_k is enlarged to carry this reset-generated matrix defect and its subsequent effect on innovation/gain/covariance. Treating reset as a common rotation would be a dropped shipping structure.

Current limiter: propagate the reset matrix defect through the next prediction/acc/S/mag Riccati operations (or prove a closed finite-dimensional augmented algebra containing n and d_perp), while keeping the vector AW source chronology linked. Then iterate A_k together with structured coefficients and the accumulated covariance defect through the complete release schedule.

Structures preserved: literal acc/S/mag gain columns, lever-arm BG column, non-isotropic S covariance, vector correction order, local-defect/state source identity, first-order MEKF covariance reset. Relaxations introduced: only the optional scalar triangle AW recurrence and optional submultiplicative reset-defect bound; exact vector/matrix paths are retained. Failure classification: no theorem failure; the prior scalar-only reset-closure premise is an E-level proof-model omission and is corrected here. No shipping counterexample.

## Structured AW gains and literal source identities exported

For an accelerometer correction, the AW gain is now represented directly in INJ coefficients:
K_aw=(P_aw,theta J_att' + P_aw,aw J_aw' + [active BA] P_aw,ba) S_acc^-1.
All products and S^-1 stay in the same INJ algebra, so K_aw is an INJ block. Its exact spectral norm is max(|i+n|,sqrt(i^2+j^2)); no entrywise K interval is introduced. For a due S=0 correction while the LIN covariance blocks are scalar*I, K_aw,S=(p_aw,S/(p_SS+r_S)) I exactly.

The source side is linked by literal identities rather than new boxes: r_acc=f_meas-[R_wb(a_w-g)+lever+b_a(temp)], and r_S=-S_mean. Existing native local-defect instrumentation already exports literal r, K, estimator state and physical truth at each accepted correction, providing a carried parity target for the structured formulas.

The prospective scalar iteration is therefore now algebraically complete in form:
A_next <= phi A + ||K_aw||||r_acc|| + ||K_aw,S||||S_mean||,
delta=delta(A,phi),
D=|b|sin(delta)+2|c|sin(delta/2).
What remains OPEN is the source-uniform numeric iteration: the structured Riccati engine must actually propagate the block coefficients P_aw,theta, P_aw,aw, P_aw,ba, P_aw,S, P_SS and S_acc through prediction/corrections/resets, while the causal history leaf supplies dependency-preserving r_acc and S_mean bounds. Do not fill these with independent global maxima.

Structures preserved: exact gain algebra and literal residual definitions. Relaxation: scalar norm triangle in A recurrence. No shipping counterexample.


## Structured AW mean recurrence linked to changing-axis defect

The literal mean chronology is now encoded prospectively. Prediction gives a_w^- = phi a_w. Accepted accelerometer correction gives a_w^+ = a_w^- + K_aw r. A due S=0 correction subsequently contributes K_aw,S (-S). Therefore a source-uniform norm image satisfies A_next <= phi A + ||K_aw||||r|| + ||K_aw,S||||S||, with all gains/residual/means required from the SAME structured history leaf. No independent A_max is introduced.

Combining with the exact literal axis transport gives delta(A,phi) from w=a_w-g -> phi a_w-g and D=|b|sin(delta)+2|c|sin(delta/2), with pi/worst-basis fallback when w can approach zero. `structured_aw_mean_recurrence.py` implements this prospective recurrence.

This calculation also corrects the role of Mahony: after goLive the private Mahony observer does not drive qref. MEKF gyro propagation and accepted attitude injections are common rotations of the structured axis and cancel exactly; watchdog re-lock is a separate reset branch. Thus the missing numerical inputs are not a Mahony angle box but SAME-history structured bounds for K_aw, acc residual, K_aw,S, S mean, and INJ b,c.

Next calculation: extend the symbolic/block Riccati step to export K_aw and K_aw,S norms directly from its scalar coefficients, and connect existing physical/local-defect source bounds to ||r|| and ||S||. Then iterate A_k,D_k alongside covariance coefficients.

Structures preserved: literal mean correction order and common rotations. Relaxations: norm triangle inequality in the scalar A recurrence; source/gain correlation should be retained where possible in the final shaped certificate. Failure would be D. No shipping counterexample.


## Literal axis transport correction: gyro/reset rotations cancel; OU latent evolution is the continuous mismatch

Source audit of the shipping accelerometer row shows f_cog,b=R_wb(a_w-g). The private Mahony observer is measurement-only after goLive and does not propagate qref; MEKF qref is propagated by bias-corrected gyro and corrected by left attitude injections. If the structured covariance axis is carried through each injection/reset by the same left rotation, both the correction/reset and the next gyro prediction are exact common conjugations and preserve INJ coefficients. Therefore the continuous inter-row axis mismatch reduces exactly to the world latent vectors w_k=a_w,k-g and w_{k+1}=phi_k a_w,k-g during prediction (before measurement correction), not to physical jerk plus an independent Mahony/gyro term. Watchdog re-lock remains a separate shipping reset event.

For one OU prediction, ||w_+-w||=(1-phi)||a_w|| and min(||w||,||w_+||)>=g-||a_w|| when ||a_w||<g. Thus delta<=2 asin((1-phi)A/[2(g-A)]) when the chord is below diameter; otherwise use pi and weighted defect. Worst timing/tau corner h=.006,tau=.02 gives phi=.740818. Diagnostic deltas: A=.1 -> .153 deg; .5 -> .798 deg; 1 ->1.686 deg; 2 ->3.805 deg; 4 ->10.243 deg; 8 ->70.04 deg. Near A~g the direction can become arbitrary, so the proof must again weight by the actual structured anisotropic coefficients rather than require uniform small delta.

This supersedes the earlier physical-jerk axis diagnostic for the covariance H-axis: physical jerk belongs to source/physical-history bounds, but the literal H nominal axis is estimator latent a_w-g. No Mahony continuous term is admissible here because shipping does not feed Mahony attitude into the live MEKF.

Next quantitative calculation: derive/propagate a numerical same-history bound on ||a_w|| (and the relevant INJ b_k,c_k coefficients) from the literal OU mean + acc/S corrections over the release image, then evaluate the weighted rebase defect |b_k|sin(delta_k)+2|c_k|sin(delta_k/2). Include watchdog re-lock as its own reset branch.

Structures preserved: literal f_cog,b model, gyro qref propagation, correction injection, OU mean. Relaxations: norm-only A bound in the diagnostic. No shipping counterexample.


## Same-history per-sample force-axis mismatch quantified without low-force singularity

Using the literal marine jerk ceiling J_max=100 m/s^3 and qualified sample gap h<=.006 s gives ||Delta f||<=.6 m/s^2 from physical jerk over one sample. For nonzero endpoint magnitudes r=min(rho_k,rho_{k+1}) with .6<2r, the exact chord bound is delta_k<=2 asin(.6/(2r)). This gives 1.84767 deg at rho=18.60665, 3.43826 deg at rho=10, 6.87963 deg at rho=5, 34.9152 deg at rho=1. Near zero the angular bound correctly becomes vacuous.

The apparent singularity is harmless when axis mismatch is weighted by the anisotropic projector coefficient b=p_theta rho^2: the worst-basis projector cost collapses quadratically as rho->0 (0.222 at rho=.3, .02467 at rho=.1, exactly 0 at rho=0 in the goLive diagnostic). Therefore the structured covariance proof should propagate weighted two-axis defects directly, not seek a uniform small delta_k.

Important remaining shipping term: the current .6 bound is for physical world-force evolution from jerk. The requested delta_k compares n_{k+1} with U_k n_k, so the final same-history bound must add/derive the frame-transport discrepancy from the literal gyro/Mahony/reference update. Do not append an independent Omega*h box; derive the relative transport from the same history. This is now the next missing term.

Structures preserved: temporal jerk, sample-gap qualification, low-force projector scaling. Relaxation: diagnostic uses the goLive b=p_theta rho^2 to illustrate weighting; later covariance b_k must be carried from the structured Riccati coefficients. No shipping counterexample.


## Changing-axis structured covariance transition derived

For exact common rotation U, the INJ algebra transports with no coefficient inflation: if n_+=Un, then U(aI+b nn'+c[n]x)U'=aI+b n_+n_+'+c[n_+]x. Thus predictor rotation itself is not a source of block-algebra growth.

For the actual same-history mismatch m=n_{k+1} versus Un_k, the exact basis differences satisfy ||nn'-mm'||2=sin(delta) and ||[n]x-[m]x||2=2 sin(delta/2). Re-basing an INJ block therefore costs at most |b|sin(delta)+2|c|sin(delta/2). For any nonzero mismatch the minimal natural two-axis span used here has rank 6 (I, nn', mm', nm'+mn', [n]x, [m]x), versus rank 3 when axes coincide. CI validates the identities.

This identifies the next quantitative bridge: derive a SAME-HISTORY per-sample/window bound on delta_k=angle(n_{k+1},U_k n_k) from the literal physical acceleration/jerk plus Mahony/gyro propagation. Do not independently box delta. If delta_k=O(h J/|f| + gyro/reference mismatch) with useful constants, structured covariance can be propagated as INJ coefficients plus a controlled two-axis defect rather than arbitrary 3x3 blocks.

Structures preserved: exact common rotation and same-history axis relation. Relaxation: none in transition identities. No shipping counterexample.


## Structured covariance class identified: INJ block algebra

The literal predictor has a natural 3-axis block structure: LIN [v,p,S,aw] is a scalar 4x4 chain tensored with I3; BA is scalar*I3; AG is a 2x2 matrix of 3-D rotation/integral operators. For a fixed shared force direction n, define N=nn' and J=[n]x. The three-dimensional algebra A_n=span{I,N,J} is closed because N^2=N, J^2=N-I, NJ=JN=0 and J'=-J. It is also closed under inversion when nonsingular: for A=iI+nN+jJ, A^-1 has coefficients i/(i^2+j^2), 1/(i+n)-i/(i^2+j^2), -j/(i^2+j^2).

Therefore an accelerometer Riccati correction preserves blockwise A_n structure PROVIDED the common force direction symbol is carried through all relevant blocks: HPH', PH', S^-1 and P-PH'S^-1HP are sums/products in the same algebra. This reduces each 3x3 covariance block from 9 arbitrary entries to 3 scalars.

Important limitation: shipping prediction rotates attitude with gyro while the physical force direction can change. Thus a single fixed n algebra is not automatically invariant across multiple samples. Full repeated closure requires a transition map between A_{n_k} and A_{n_{k+1}} using the SAME Mahony/physical history and MARINE bounds, not independently reselecting n each sample. The theorem-status flag closes only the algebraic primitive, not full shipping chronology.

Structures preserved: common force symbol, exact skew/projector identities, literal block predictor topology. Relaxation: none in algebra closure; cross-sample changing-axis closure remains open. No shipping counterexample.


## Direction-free goLive accelerometer innovation inverse CLOSED

From the exact common-R Gram identity, write f=rho n and P_n=I-nn'. Then S=cI+aP_n with c=R_acc+p_aw+p_ba and a=p_theta rho^2. Since P_n is an orthogonal projector,
S^-1 = c^-1 nn' + (c+a)^-1(I-nn').
This is exact and eliminates both attitude and force-direction partitions at goLive. Over rho in [0,18.60665] and literal isotropic R_acc variance [0.04,0.090625], the verified eigen floor is c>=16.520016 and ||S^-1||2<=0.0605326291. At rho_max and minimum R, transverse inverse eigenvalue is .00114842656 while parallel is .0605326290. The theorem-status flag `golive_direction_free_projector_inverse` is closed.

This closes only the literal block-diagonal goLive correction. After one prediction/correction, P develops cross-block structure, so the same projector formula cannot simply be reused. The next calculation should derive whether the shipping prediction preserves a structured covariance class (isotropic 3x3 Kronecker/block form under common rotations) sufficiently to propagate a symbolic/block Riccati envelope without reverting to entrywise P/K intervals.

Structures preserved: exact projector/eigenstructure, full force range, isotropic R_acc. Relaxations: none in the inverse identity. No shipping counterexample.


## Common-R block Gram identity removes attitude partition at goLive

Keeping the shared R symbol through HPH' before intervalization yields an exact simplification for the literal block-diagonal goLive covariance:
HPH' = p_theta (||f||^2 I - f f') + (p_aw+p_ba) I.
The identities [f]x[f]x'=||f||^2I-ff' and R R'=I eliminate ALL R_wb dependence. Thus the earlier 1021 attitude cells are unnecessary for the goLive accelerometer innovation itself; their near-threshold residual was purely an enclosure artifact.

The blockwise probe verifies the exact-force case with residual 0.001529888 and inverse norm bound 0.06053263. Force uncertainty remains: 0.1-deg direction gives residual 1.08287 while .01 deg gives .109664 under the current Lipschitz force-Gram bound. The next tightening should exploit the eigenstructure of A(f)=||f||^2I-ff' directly (eigenvalues 0,rho^2,rho^2) rather than entrywise/spectral perturbation about a force direction. Since innovation eigenvalues depend on rho but orientation only rotates the rank-one null direction, an isotropic R_acc may allow a direction-free inverse bound and eliminate force angular covering too.

Structures preserved: literal block-diagonal goLive P, common R, exact skew identity. Relaxations: current force-cone Lipschitz perturbation is conservative. No shipping counterexample.


## Rotated accelerometer innovation frame tested: algebraic AW cancellation does not improve the factor enclosure

The exact transformation Q=R_wb' gives Q J_aw=I, Q J_ba=R_wb', and Q J_att=-[f_b]x R_wb'. This is algebraically lossless with isotropic R_acc and was implemented before intervalization. However, over the complete 1021-cell attitude cover the factor certificate is worse: even exact worst-magnitude local-body force has worst residual 7.955848 (Gram perturbation 131.606807), with all 1021 cells failing; the full force ball gives residual 340.418. The attitude uncertainty removed from AW reappears in BA and the right factor of J_att, and the current interval product pays for those dependencies separately.

Classification D: useful algebraic transformation, but no useful enclosure under the present factor norm. Do not pursue this rotated-frame route as controlling proof. Retained fact: Q J_aw=I exactly, but cancellation alone is insufficient. No shipping counterexample.

The next stronger route should preserve the shared rotation symbol across BOTH transformed attitude and BA blocks (symbolic/affine arithmetic or a block quadratic Gram identity), rather than replacing R' by independent entry radii in each occurrence.


## Joint attitude-cover x force-cone threshold is numerically impractical with current factor norm

Using the full verified 1021-cell attitude cover, the worst-cell joint threshold was solved before generating a force mesh. At exact force direction, the outer radial-shell width must be below roughly 0.0034--0.0040 m/s^2 (passing 0.00340698 gives residual .998298; failing .00397481 gives 1.003819). At exact force magnitude, force-direction half-angle must be below roughly .00549--.00610 deg (passing .00549316 gives residual .999944; failing .00610352 gives 1.003807). The worst attitude cells are near the factor threshold, so almost no force uncertainty remains.

Therefore a literal Cartesian product of 1021 attitude cells with a force spherical cover at this resolution would be enormous and is not the preferred constructive route. This is a D-level conservatism of the current factor perturbation norm, not a shipping issue. The next calculation should tighten the JOINT attitude+force geometry rather than independently budgeting their perturbations: exploit that accelerometer attitude block -[f]x and AW block R_wb arise from the same physical/attitude history, and preferably rotate the innovation into the local body/force frame before taking norms. Do not generate millions of force cones from these thresholds.

Structures preserved: complete attitude cover and exact force-cone identities. Relaxation causing failure: additive independent norm budget for attitude-induced R uncertainty and force-induced skew uncertainty inside Y. No shipping counterexample.


## Captured 6.9-deg attitude ball factor cover CLOSED

A deterministic cubic-lattice cover in rotation-vector coordinates now covers all of R^3 with Euclidean covering radius 1.10 deg; retaining cells whose 1.10-deg balls intersect the captured 6.9-deg ball gives 1021 local attitude cells. The analytic full-space lattice covering argument proves no holes; the local radius has 0.05147686 deg margin below the measured factor threshold 1.15147686 deg.

Each cell uses its own Rodrigues center R(theta_c) and multiplicative local rotation perturbation, not a zero-centered approximation. Running the PSD-factor innovation certificate over all 1021 cells with exact worst-magnitude force produced zero failures. Worst cell: residual 0.965173718, Gram perturbation 15.94380372, inverse norm bound 1.73547152. Thus the complete captured attitude domain admits a narrow factor-innovation inverse under this local cover.

Caveat: this run used exact force direction/magnitude to isolate and close the attitude cover. The final joint history cover must combine these attitude cells with force-direction/magnitude cells; the earlier whole-force cone failure means that joint refinement remains necessary. This does not change the 6.9-deg theorem domain.

Structures preserved: full captured attitude ball, center-dependent SO(3), PSD factor innovation. Relaxation: lattice cells overlap and include points outside the captured ball, a safe outer cover. No shipping counterexample.


## Force-angle threshold test: no alpha_* exists until attitude geometry is refined

The requested residual-vs-force-cone-angle calculation was completed with rho fixed exactly at 18.60665 m/s^2. There is no passing force half-angle under the current 6.9-deg attitude ball: even alpha=0 gives Gram perturbation 105.695718 and residual 6.389781>1. Therefore no force-direction spherical cover, however fine, can close the current factor inverse certificate.

A follow-up isolation sweep with exact force direction/magnitude shows the residual crosses 1 as the attitude/rotation cell radius shrinks. Verified bisection gives beta_* in [1.1514768600, 1.1514778137] deg, with residuals 0.999999662 and 1.000000500 respectively. Thus the current factor certificate requires attitude-history cells of radius <~1.1514769 deg. This does NOT strengthen the physical captured assumption: the full 6.9-deg captured ball must be covered by multiple attitude cells whose union contains it, with centers generated from the same history geometry. A fixed smaller attitude assumption is forbidden.

Structures preserved: exact force skew norm and PSD-factor innovation. Relaxation: factor perturbation norm remains conservative. Failure classification D. No shipping counterexample. Next constructive calculation: cover the 6.9-deg rotation-vector ball by <=1.15-deg local attitude cells (or derive center-dependent SO(3) cells) and retest factor inverse jointly with force cones.


## PSD-factor innovation refinement localizes remaining width to force geometry

The release innovation has been reformulated as P=L L', Y=H L, S=R_acc+Y Y'. `psd_factor_innovation.py` bounds the Gram perturbation at Y-level before inversion. On the entire captured attitude ball plus entire force ball, the first factor probe gives Gram perturbation inf-norm 4088.3910 and Neumann residual 247.10397, so a narrow inverse is not available. This is class D and sharply localizes the width to H/force geometry; covariance PSD and inverse existence are no longer the issue.

An adaptive force-ball sector cover was launched, but Cartesian force boxes are expensive and discard the norm-ball geometry. If that cover is large, the next preferred calculation is an analytic magnitude+direction/conic force parameterization for [f]x, preserving ||f||<=g+Amax directly rather than enumerating an axis cube.

Structures preserved: P factor, Y=HL, PSD Gram, attitude ball, force norm premise. Relaxation: interval Y factor radii remain conservative. No shipping counterexample.


## Repeated release propagation: invertibility solved, gain enclosure is the next limiter

The staged goLive->A21 PSD-aware propagation was executed after the one-step innovation gate closed. The center-zero inverse fallback causes interval Joseph arithmetic to become non-finite before the first useful staged checkpoint. This is class D: the inverse existence theorem is valid, but the entrywise gain enclosure is too coarse.

The fallback has therefore been tightened without changing the admissible shipping set: retain the actual midpoint innovation inverse as center and use the PSD/noise-floor theorem only for a rigorous perturbation radius. The primary Neumann/Krawczyk inverse remains preferred whenever its interval residual verifies. A new staged run is queued to determine whether this centered PSD enclosure survives 10/100/1000+ samples.

Structures preserved: P>=0, actual midpoint S, measurement floor, literal Joseph form. Relaxation: coarse norm perturbation radius only; failure is D. No shipping counterexample.


## Release innovation gate CLOSED source-uniformly on the captured ball

The captured attitude domain is now represented by the theorem's rotation-vector ball ||dtheta||<=6.9 deg, not the artificial Cartesian cube. Rodrigues bounds enclose the whole SO(3) image. A force-ball cover retains ||f||<=g+Amax.

Entrywise innovation inversion still failed after force subdivision because it discarded the structural PSD fact. The shipping covariance set satisfies P>=0 at every Riccati step, hence HPH'>=0 and S=HPH'+R_acc >= R_acc. `verified_joseph_update_psd` now uses this lossless covariance-set invariant: if the narrow approximate-inverse residual certificate fails, the strict measurement-noise floor supplies a rigorous center-zero inverse entry box with ||S^-1||<=1/lambda_min(R_acc). No midpoint gain is used.

CI result for the complete one-step captured attitude/force domain: VERIFIED, leaf_count=1, unresolved=0, max_depth=0. Thus no geometry subdivision is needed to prove innovation invertibility. The next release calculation is repeated PSD-aware prediction/correction propagation from goLive toward A21 activation; gain tightness, rather than invertibility, may become the limiter.

Structures preserved: attitude norm ball, force norm ball, covariance PSD invariant, literal R_acc floor. Relaxations: the noise-floor inverse entry box is coarse but lossless. Failure of its downstream Joseph enclosure would be D. No shipping counterexample.


## First release interval probe: broad H_acc box fails verified innovation inverse

After two CI-only class-E issues were fixed (workflow heredoc indentation and literal escaped newlines in theorem_status.py), the first actual one-sample goLive release probe ran. Result: `one_step_verified=false`, reason `innovation inverse not verified`. This is class D: the existing accelerometer H interval independently boxes skew-force and rotation entries and is too broad even at the numerical goLive covariance seed. It is not a shipping counterexample and does not invalidate the goLive-image construction.

The next refinement is shipping dependency, not smaller assumptions. `accel_geometry_cell.py` now constructs H_acc from one shared force-vector cell and one SO(3) rotation cell; J_att=-[f]x, J_aw=R_wb and J_ba=I remain coupled. Small-angle rotation cells include a rigorous quadratic remainder and fail closed when too wide. The release cover must subdivide the shared attitude/specific-force history coordinates until the innovation inverse verifies, rather than independently splitting H or K.

Structures preserved: literal acc geometry and common attitude/force dependency. Relaxation: local SO(3) enclosure is conservative with a rigorous remainder. Failure classification D. No shipping counterexample.


## goLive -> A21 finite release image propagation is now executable

`golive_to_a21_release.py` propagates the explicit zero-mean goLive seed through a supplied literal operation stream, carrying the full verified 21-state covariance and affine mean interval and emitting the A21 BA-graph enclosure. `release_interval_propagation.py` supplies the shipping-specific interval prediction and accelerometer-correction covariance map from the declared dt/tau/omega/noise ranges; it retains the AW covariance floor and deliberately does not invent S or magnetic callback cadences. Aggregate magnetic information remains a theorem action, not an event schedule.

The release construction no longer needs an existential standalone LIN radius: the numerical LIN release box is the finite image of the zero goLive LIN mean. The remaining source-uniform integration is coefficient-cell subdivision for the literal operation stream. A dedicated CI probe now tests whether the broad root coefficient cell can verify even one innovation inverse before attempting the 469-s image. If it fails, the correct response is shared-history subdivision (attitude/specific-force/tuner coordinates), not midpoint K or an independent gain box.

Structures preserved: literal goLive seed, full P cross-covariance propagation, acc correction chronology, AW sync increase, aggregate magnetic service. Relaxations: optional S/mag covariance reductions are omitted only for the release covariance upper comparison; their mean/action effects remain required in the literal operation stream used for the final leaf. No shipping counterexample.


## Explicit release numerics: AG/BG/BA closed; LIN radius eliminated as a prerequisite

The literal H18/refinement horizon certificate gives a conservative A21 activation time of 469 s from captured-service origin. Using the worst shipping attitude reseed variance 1.5708^2, Pb0=1e-6, AtomS3R shipping gyro noise density .00135 rad/s/sqrt(Hz-equivalent discrete input) at 200 Hz, and gyro-bias driving variance 1e-11, an observation-free AG covariance comparison gives P_AG <= 8.06571980866 I_6 over that finite horizon. This deliberately drops covariance-reducing acc/mag corrections. The shipping gyro-bias projection radius .5 plus physical slow-bias radius .02 gives ||e_bg||<=.52 rad/s. Since the active-A21 accelerometer Jacobian has J_ba=I, the zero-action BA graph obeys ||A_ba||<=g+A_max=18.60665, hence each graph entry lies in [-18.60665,18.60665].

A separate numerical LIN mean radius is no longer required to seed the constructive proof. Shipping bootstrap leaves the MEKF mean untouched; at goLive all 21 mean coordinates are zero, while covariance has an explicit diagonal outer seed: attitude <=1.5708^2, BG=1e-6, v=1, p=400, S=2500, AW<=16.48, BA=.004^2. `golive_release_seed.py` records this. The release set is therefore defined constructively as the finite literal image of this numerical goLive seed under the same 469-s causal history propagation. This is stronger than converting qualitative H18 BIBO compactness into an arbitrary LIN radius.

MAGNETIC SERVICE is now consumed directly as its actual transported/innovation-whitened information Gramian; no extra magnetic gain weight is needed and no callback schedule is enumerated.

Open integration: the sample history propagator must now propagate the goLive mean/covariance through literal MEKF prediction/corrections to the A21 activation event. Until that finite image is interval-computed, `release_state_enclosure_connected_to_history_cover` remains false. No missing standalone LIN constant remains as a theorem premise.

Structures preserved: literal goLive initialization, finite H18/refinement chronology, aggregate actual-innovation MAGNETIC SERVICE. Relaxations: AG covariance comparison drops corrections only in the safe covariance-increasing direction. No shipping counterexample.


## Aggregate magnetic service connected; explicit release-set audit exposes the true remaining quantitative gap

`aggregate_magnetic_service.py` now consumes MAGNETIC SERVICE exactly as stated: `sum G_i'G_i >= mu_M I_2` over each service window. It produces a root-coordinate PSD action lower bound once the SAME-history magnetic loss weight is certified; no callback schedule or service-row enumeration is introduced. `literal_history_leaf_propagator.py` now carries this aggregate operator instead of synthetic magnetic event boxes. Thus the prior event-schedule blocker is removed by using more shipping/theorem structure, not less.

`explicit_release_set.py` audits the H18/A21 release theorem against the needs of a constructive leaf. The existing repository proves captured-domain release compactness, BG/LIN mean compactness and a numerical LIN/BA covariance upper comparison, but it does NOT export numerical AG covariance upper bounds, BG/LIN mean radii, or a BA-graph interval. Consequently there is currently no rigorous explicit full-P0 interval seed. This is now the single release connector blocker. Qualitative compactness is deliberately not converted into an invented box.

Structures preserved: aggregate MAGNETIC SERVICE and authoritative release theorem. Relaxations introduced: none. Failure classification: missing quantitative release constants is an open proof obligation; treating compactness as numerical bounds would be E. No shipping counterexample.

## HistoryCell sample propagator implemented; two required release/service inputs exposed

`history_interpolation.py` rigorously expands sparse theorem-input knots: SLOW uses all knot/rate constraints; physical p/v use derivative envelopes; FAST uses the all-window signed primitive cap directly; gravity/magnetic vectors remain full envelopes between sparse knots unless causal information narrows them. `literal_history_leaf_propagator.py` expands a 60/100-s cell to delivered 5-ms samples and runs the closed Mahony->WavePeriod->band/variance->staged tuner chain sample-by-sample. It emits timed physical/gravity/magnetic boxes and the generated adaptation trace.

The integration attempt exposed that the current source-uniform root is insufficient to construct a complete `ConstructiveLeaf` without inventing generated initial/service state. Two authoritative inputs must be connected from already-existing proof machinery: (1) the certified compact A21 release-state enclosure, including initial P and the word-dependent BA compatibility graph; (2) accepted MAGNETIC SERVICE event strata/rows, not arbitrary magnetic-field knots. The propagator therefore returns `complete_constructive_leaf=false` until these are supplied. This is fail-closed shipping faithfulness, not a new mathematical blocker. Theorem status records `literal_history_cell_sample_propagator=true` but the two connections false.

Structures preserved: temporal interpolation outer cover and literal sample-by-sample adaptation. Relaxations: sparse gravity/magnetic knots are deliberately not used to certify excitation/service. Failure classification E/integration for missing authoritative release/service connectors, not D and not shipping instability. No shipping counterexample.

## Source-uniform causal history-cover generator implemented

`source_uniform_history_cover.py` now performs adaptive fail-closed branch-and-bound exclusively on shared theorem-input history coordinates. Generated Mahony/frequency/variance/tuner/covariance/gain/scheduler/reference/BA states cannot be split. Each cell is delegated to one causal leaf propagator and the completed nonlinear/kernel/QCQP leaf certificate; verified leaves contribute their rigorous ratio, unresolved leaves subdivide the most causally sensitive admissible input, and `max_ratio` emits a number only when the entire admissible cover verifies. `source_uniform_root.py` supplies finite-dimensional 60/100-s conditional root parameterizations in physical/SLOW/FAST/delivered-vector knot coordinates, never generated filter coordinates.

A final integration distinction is now explicit in theorem status: the cover GENERATOR exists, but `literal_history_cell_to_leaf_propagator` remains false. We have implemented its submaps (sensor->Mahony/WavePeriod/tuner; interval covariance/gain; temporal source domain; witness builder), but not one single function that takes a root knot cell, reconstructs/encloses every delivered sample over 60/100 s, applies temporal interpolation/jerk/MARINE constraints, and emits the complete `ConstructiveLeaf`. Without that function the cover cannot actually run over the theorem class, so no numerical max_i s_i/d_i is claimed. This is integration work, not a new proof architecture.

Structures preserved: shared input-only subdivision and all generated-state causality. Relaxation: finite-dimensional knot parameterization must be proved to outer-cover all admissible continuous/held histories using jerk/SLOW/FAST constraints; until that interpolation theorem is wired, it is non-promoting. Failure classification D for an overly coarse knot cover, E for propagation mismatch. No shipping counterexample.

## Executable history-cover pipeline connected; interval B/C midpoint escape removed

The MARINE span callback now implements the theorem's existential diameter logic per complete window: at least one certified endpoint pair in each window suffices; UNKNOWN candidates force subdivision. `constructive_history_cover.py` connects one causal leaf witness, interval joint form, nonlinear premises, angle-aware kernel restriction and direct source QCQP. During this connection a midpoint escape was found: the QCQP consumed midpoint B/C. `verified_interval_supply_bnb` now adds rigorous coefficient-radius contributions and the cover consumes interval B/C directly. No midpoint source coefficients remain in the executable path.

The architecture can now compute a ratio for a fully populated leaf. The remaining blocker to an actual numeric cover is data population: there is not yet a source-uniform history-cell generator that emits a finite cover of timed physical/gravity/magnetic boxes, interval A0/T/A1, compatibility graph line boxes, interval joint form, and e/u boxes over the complete 60/100-s theorem class. Existing carried replays are not such a cover. Therefore max_i s_i/d_i is still not numerically claimed. Building that finite causal cover is now the sole integration obligation; any leaf that fails inverse/kernel/nonlinear/QCQP certification subdivides only its shared history inputs.

Structures preserved: existential MARINE semantics, interval source coefficients, shared leaf dependencies, all prior shipping couplings. Relaxations introduced: branch-and-bound boxes only. No shipping counterexample.

## Causal history-cell witness builder implemented; interval B/C QCQP repaired

`causal_history_witness.py` now carries one dependency token through physical v/p/a/j and gravity-direction boxes, applied magnetic-event boxes, word-dependent attitude/BA compatibility-line cones, temporal source symbols, nonlinear MARINE/MAGNETIC premise certification, kernel-restricted action and direct QCQP supply. It derives 30-s witness pairs from the propagated physical boxes and never fabricates excitation when a whole-box witness is unavailable. `certify_leaf` emits a ratio only when both certificate halves verify on the same token.

While joining the halves, a remaining midpoint relaxation was found in the QCQP numerator: interval B/C coefficient radii were not included. `verified_interval_supply_bnb` now rigorously adds those coefficient radii to every branch-and-bound upper bound, and the causal witness uses it. Thus midpoint B/C no longer blocks promotion.

The remaining blocker to an actual numerical leaf ratio is no longer witness plumbing. We need a concrete source-uniform history-cell generator that feeds this builder with interval physical trajectories, source-symbol boxes, A0/T/A1 shaped actions, and compatibility-line graph/radius enclosures from the causal shipping propagator. The APIs now join correctly, but no synthetic or replay leaf is promoted as the source-uniform cover. `actual_history_leaf_ratio_certificate` therefore remains false.

Structures preserved: same causal dependency token across nonlinear premises, kernel action and interval B/C supply. Relaxations introduced: none in the join; branch boxes are outer representations of the same history cell. No shipping counterexample.

## Causal history-cell witness builder implemented

`history_witness_builder.py` now constructs, from one propagated dependency root, the nonlinear MARINE/MAGNETIC witness and the joint attitude/BA compatibility-line cone. It carries timed interval v,p,a,j and gravity-direction samples, all certainly 30-s-separated candidate endpoint pairs, applied magnetic event boxes, and the graph line r=(a,-G a). The line component radius includes `|G| da + dG |a| + dG da`; cells whose line cone reaches zero fail closed. `history_leaf_certificate.py` connects that SAME witness to nonlinear premise certification, later-word angle-aware kernel restriction, and the cellwise ratio composer.

This completes the witness/connector plumbing, but actual rigorous ratios still require the upstream causal history propagator to emit these timed physical/gravity/magnetic boxes and the source QCQP to finish a verified supply certificate on each leaf. The current builder's 30-s pair set is a candidate set; the callback presently requires every supplied pair to satisfy the span, which is sufficient but stronger than the theorem's existential diameter condition within each window. Before promotion, the builder must group candidate pairs per 30-s window and certify at least one pair per window, rather than all pairs. This is a D-level conservatism, not a shipping claim.

Structures preserved: shared dependency root, physical chronology, magnetic event chronology, joint BA graph line. Relaxation: current all-pairs span test is conservative and must be changed to per-window existential witness groups before promotion. No shipping counterexample.

## Nonlinear MARINE/MAGNETIC QCQP callbacks implemented

`marine_magnetic_qcqp.py` now provides whole-history-box TRUE/FALSE/UNKNOWN certification for Euclidean velocity/position/acceleration/jerk envelopes, 30-s gravity-direction span witnesses, 30-s displacement-diameter witnesses, and recurring 1-s magnetic service with field/residual/informative-event qualification. UNKNOWN forces subdivision; midpoint satisfaction cannot promote a leaf. `verified_linked_qcqp.shipping_side_callback` wires these premises into the direct temporal QCQP.

One concrete plumbing obligation remains before actual leaves can simultaneously emit d_i and s_i: the shared-history propagator must construct `MarineMagneticWitness` objects from each causal cell, including certified witness endpoint pairs for every complete 30-s MOVING window and the applied magnetic-event time/field/residual boxes. The callback itself is implemented; theorem status keeps `history_cell_nonlinear_witness_builder=false` until that causal builder exists. Likewise the angle-aware kernel primitive still needs actual cell line cones to set `kernel_restricted_interval_action_numeric=true`. No numerical max_i s_i/d_i is claimed before both are supplied by the same leaves.

Structures preserved: nonlinear MARINE/MAGNETIC premises and actual recurring service logic. Relaxation: none promoted; witness-pair selection must itself be certified by the history builder. No shipping counterexample.

## Angle-aware projector/restriction enclosure implemented

For a history cell whose word-dependent compatibility line lies within angle delta of its certified midpoint direction, `angle_projector_enclosure.py` uses the exact rank-one projector identity `||P(r)-P(r0)||_2=sin(delta)` and a verified interval action norm to lower-bound the homogeneous action uniformly over the entire line cone. No uncertain orthonormal basis is selected. The bound is `lambda0(1-eps^2)-||A||(2 eps+eps^2)`, eps=sin(delta), where lambda0 is the interval midpoint-complement floor. `restrict_interval_action` now consumes a nonzero line component radius through this cone certificate and fails closed if the resulting floor is not positive. Cell-ratio plumbing accepts the resulting line-cone certificate.

This closes the angle-aware restriction primitive. `kernel_restricted_interval_action_numeric` remains false until actual reachable history cells provide their interval A_super and compatibility-line radii and obtain positive floors; no carried midpoint result is promoted. The remaining independent constructive blocker is still the nonlinear MARINE/MAGNETIC QCQP callback implementation.

Structures preserved: entire word-dependent line cone, later-word action, common-root transport. Relaxation: the projector perturbation inequality is conservative but valid for the shipping line family; failure is D, not instability. No shipping counterexample.

## Word-dependent kernel-restricted interval action implemented

`kernel_restricted_action.py` now forms the superword homogeneous action in common root coordinates, `A_super=A0+T' A1 T`, before restricting by the first complete word's actual compatibility line. The complement is constructed from that word-specific line; no fixed global compatibility vector is used. `compatibility_line_interval.py` encodes the joint attitude/BA graph line r=(a,-A_ba a) and provides an angle-radius certificate for cell uncertainty. The restriction deliberately fails closed when the line has nonzero interval radius until an angle-aware projector enclosure is supplied. `constructive_cell_ratio.py` is wired to accept this later-word restricted action.

This removes the conceptual kernel blocker but not yet the numeric cell certificate: actual history cells produce an interval family of compatibility lines, not an exact line. The next calculation is therefore the angle-aware projector/restriction enclosure using the certified line cone, so `kernel_restricted_interval_action_numeric` can become true without midpoint-line substitution. After that, the nonlinear MARINE/MAGNETIC QCQP callbacks remain the other blocker to max_i s_i/d_i.

Structures preserved: word-dependent joint attitude/BA line, rotating weak directions, later-word action and common-root transport. Relaxation: none promoted; midpoint line is diagnostic only when line radius is zero. No shipping counterexample.

## Direct temporal QCQP selected; cellwise max s_i/d_i plumbing complete

A direct verified branch-and-bound QCQP now bounds the linked B/C supply on the SAME temporal source variables, avoiding an artificial ellipsoid. Linear SLOW/FAST window facets prune source boxes; nonlinear MARINE/MAGNETIC constraints are callback obligations and UNKNOWN forces subdivision/failure. `constructive_cell_ratio.py` composes a supply certificate with a kernel-restricted homogeneous action and computes s_i/d_i; `cover_max` returns max_i s_i/d_i only when every leaf verifies.

Two concrete blockers remain before a numerical radius can be honestly emitted: (1) construct the interval homogeneous action on the quotient/restriction that removes the word-dependent joint attitude/BA compatibility line using the already-proved later-word/A* action; full-space Gershgorin is explicitly rejected; (2) implement the nonlinear MARINE/MAGNETIC QCQP callbacks (Euclidean norms, 30-s spans, jerk, recurring service) on the same history boxes. Until both exist, the solver fails closed and no numerical max ratio is claimed.

Structures preserved: direct temporal polytope, shared source symbols, interval shaped B/C, kernel-before-division discipline. Relaxation: box branch-and-bound is an outer numerical representation of the same history domain; no generated coefficient boxes or independent source extrema. Failure from leaf limits/coarse boxes is D, callback/arithmetic defects E. No shipping counterexample.

## Temporal MARINE + SLOW+FAST source domain attached to shaped symbols

`temporal_source_domain.py` now builds one persistent source timeline. SLOW accel/gyro symbols obey amplitude plus all pairwise `min(2B,D dt)` difference constraints. FAST accel/gyro symbols obey instantaneous theorem envelopes plus every placed discrete signed-integral window up to H=60 s with caps C_a=.05 m/s and C_g=.002 rad; histories are not reset at proof boundaries. Physical v,p,a boundary symbols use the declared marine envelope, while Euclidean norm coupling, 30-s attitude/displacement spans, jerk and 1-s magnetic service are retained as nonlinear same-cell side constraints rather than replaced by component boxes. Measurement-model limits are attached from constants.json. Hardware qualification remains conditional.

`kernel_linked_quotient.py` adds the verified linked Schur primitive for the SAME interval joint form and SAME source metric; it never divides by an unrelated global lambda_min. However, the complete constructive quotient is not yet promoted because the temporal domain is currently represented by linear window constraints plus nonlinear side constraints, not yet by a verified quadratic source metric R on the exact joint source-symbol basis expected by the Schur solver. Converting this polytope/nonlinear domain to a dependency-preserving quadratic certificate (or solving the QCQP directly) is the next calculation. An independently diagonal R is forbidden.

Structures preserved: temporal SLOW/FAST, MARINE physical boundaries, magnetic service metadata, shared source symbols and full interval shaped form. Relaxation: component linear inequalities are necessary outer facets while Euclidean/span constraints remain explicit side constraints; they are not used alone for promotion. Failure of a coarse quadratic outer certificate is D. No shipping counterexample.

## Midpoint-Q relaxation removed; first rigorous interval shaped form exists

The covariance/gain path now computes a verified full 21x21 precision enclosure `Q=P^-1` at every storage boundary using the exact-binary64 residual/Neumann inverse certificate. `IJointQuadratic` accumulates interval `Q0,Q1,A,D` directly, so the shaped action coefficients no longer depend on midpoint precision. Failure of the full inverse certificate stops the proof. Tests require verified full precision and interval joint algebra.

This is the first genuinely rigorous shaped-form infrastructure, but the final constructive certificate is intentionally still fail-closed: the interval A block can be inspected for homogeneous action, while the B/C source blocks cannot be converted to a supply bound until the SAME history cell carries a rigorous source-domain metric encoding SLOW rates, FAST primitives, MARINE physical boundaries, measurement/model defects and magnetic-service source coordinates. `rigorous_cell_from_interval_joint` therefore refuses to promote a quotient without that metric. This avoids reintroducing independent source suprema.

Structures preserved: verified P/K/Q chronology, full 21-state precision, shared history source symbols, interval joint cross terms. Relaxations introduced: none in precision inversion; source-domain metric remains an open constructive obligation rather than a relaxation. No shipping counterexample. Next calculation: attach the theorem's temporal SLOW+FAST/MARINE source-domain quadratic/functional to the same history symbols and solve the kernel-aware generalized cell quotient.

## Literal covariance/gain -> joint shaped-supply subchain implemented

`causal_covariance_supply.py` now takes one shared reachable-history cell and an A21 covariance enclosure, propagates literal prediction P-=FPF'+Q, computes each 3-D innovation and a verified interval inverse, generates K=P H' S^-1 internally, applies literal Joseph covariance, and feeds A=I-KH plus the same history-source columns immediately into `JointQuadratic`. Reset covariance is transported congruently. AW covariance sync and BA release remain covariance-only chronology events and do not invent mean supply. P and K are never admissible history-cell coordinates.

This closes the covariance/gain plumbing but not yet the complete constructive certificate. The current joint accumulator uses the covariance midpoint precision for its numerical quadratic coefficients while P/K themselves are interval-certified. Promotion therefore still requires interval/affine enclosure of Q=Pbar^-1 and of the joint quadratic coefficients, plus literal history-derived H/source factors for acc/mag/S and the exact OU process-Q constructor driven by the generated tau/sigma tuple. Until those are enclosed, this is proof infrastructure, not a source-uniform numerical margin.

Structures preserved: literal covariance/gain chronology, Joseph updates, covariance-only sync/release, shared source dependency token. Relaxation: midpoint Q inside the provisional joint numerical accumulator; any margin from it is non-promoting. Failure class for that provisional restriction is D; inverse/parity implementation failures E. No shipping counterexample.

## Sensor history -> joint tuner chain causally closed

The reachable-history interval propagator now generates `f_tune` internally through the literal `WavePeriodEstimator` recurrence driven by the Mahony vertical signal. It carries both high-pass stages, leaky velocity/elevation states, period-scaled moment EMAs, moment-ratio inversion, canonical log-period EMA, usable-period one-way gate, and fixed-prior fallback until that gate qualifies. `ClosedCausalAdaptationBox.step()` therefore has no wave-frequency input: one delivered gyro/accelerometer history generates Mahony -> vertical -> WavePeriodEstimator -> f_tune -> adaptive band -> variance -> staged tau/sigma/R_S/T_S. Ambiguous moment-start, variance, omega, horizon-clamp, log-horizon and usable-period gates fail closed and require subdivision of the same sensor-history cell.

Structures preserved: literal measurement-only frontend timing and full joint adaptation chain. Relaxations introduced: interval outer arithmetic on delivered-history cells; no generated frequency/tuner coordinate is independently boxed. Remaining causal-history work is downstream covariance/gain propagation and joint shaped supply, plus stillness/vibration strata where the theorem cell reaches them. No shipping counterexample.

## Causal Mahony -> variance -> tuner history propagation implemented

`causal_tuner_interval.py` now propagates one shared delivered gyro/accelerometer/wave-frequency history cell through the private Mahony quaternion/integral feedback, vertical acceleration, adaptive wave-band state and white-noise covariance, debiased first/second variance moments, operating-point targets, common tau/sigma EMA, deployed SpectralMSE R_S target, R_S smoothing, one-sample staged commit and T_S(tau) cadence. Nonlinear normalization and clamp-stratum ambiguity fail closed and require subdivision of the shared input history, never generated outputs. The implementation preserves the shipping timing that candidate y_k is committed before y_{k+1}. Tests cover point-band arithmetic, nonnegative variance, staging, generated R_S linkage and fail-closed Mahony normalization.

This closes the first missing causal submap, not the complete history propagator. Still open: rigorous WavePeriodEstimator frequency-cell propagation from the Mahony vertical signal (rather than supplying a shared wave-frequency input coordinate), stillness/gate strata, vibration conditioning where active, and then covariance/gain propagation into the joint shaped form.

Structures preserved: measurement-only Mahony, adaptive band, variance moments, SpectralMSE law, joint tau/sigma/R_S/T_S, one-sample staging. Relaxation: wave-frequency is presently a shared upstream history coordinate pending literal WavePeriodEstimator enclosure; it is not independently boxed from the same history in the final theorem. Failure of a coarse interval is D; arithmetic/chronology mismatch E. No shipping counterexample.

## Constructive reachable-graph enclosure is now the controlling calculation

The new interface in `reachable_history_enclosure.py` permits subdivision only in shared theorem-input history coordinates and structurally rejects generated tuner/covariance/gain coordinates. `joint_shaped_supply.py` accumulates the exact same-symbol quadratic form for storage loss minus physical supply, retaining source/loss cross terms. The constructive certificate uses the cellwise linked quotient `max_i s_i/d_i`, never `(max s_i)/(min d_i)`. The remaining open implementation is the literal causal interval/affine propagator from one history cell through frontend/Mahony, frequency/variance, staged tuner, covariance/gains and mean defects, followed by every-prefix joint-form enclosure. Existing interval matrix arithmetic may be reused only underneath this history interface.

Structures preserved: complete shipping causal graph, source/loss dependencies, temporal SLOW+FAST primitives, OU/S chain and adaptive chronology. Relaxation: numerical outer enclosure of theorem-input history cells only. Failure from coarse enclosure is D; propagator/parity defects are E. No shipping counterexample found.

## Common-Q and separated-LIN relaxations removed

The shaped-storage diagnostic now uses the full dense 12x12 LIN block of the
literal covariance metric after shipping normalization and permits root/terminal
storage to vary with the actual generated tuner/covariance state.  The earlier
seven scalar block weights and one-common-Q-across-sea-states ansatz are retired,
not controlling.  The complete carried test is
rho=lambda_max(Q0^-1/2 Mbar' QN Mbar Q0^-1/2) with Qk=(Dk Pk Dk')^-1.
This retains v/p/S/a_w cross terms and causal adaptation.  Remaining relaxation:
finite carried replay rather than source-uniform enclosure.  Failure remains
class D unless a literal shipping execution violates a theorem conclusion.

## Adaptive shaped-storage derivation advanced

The exact lossless coordinate identities are now implemented in adaptive_shaped_storage.py and regressed: the shipping integrated OU prediction normalizes to a dimensionless h/tau family, and generated tuner scaling changes are exact coboundaries that telescope across a persistent history. The carried source driver now exports the actually applied tau, sigma_aw, R_S, T_S, targets, variance/frequency state and proxy quaternion at each recorded sample; adaptive_shaped_storage_diagnostic.py consumes those joint tuples without Cartesian boxing. This is the first quantitative route that explicitly keeps the shipping tuner generator inside W_out.

Structures preserved: literal shipping OU primitives; generated tuner tuple; one-sample causal adaptation; persistent Mahony-driven frontend history; persistent filter/physical execution.
Relaxations introduced: current shaped-storage diagnostic is finite carried evidence only; common quadratic Q has not yet been asserted source-uniformly.
Failed calculations: none in this shaped-storage derivation yet. Earlier escaped-newline/tooling faults are class E; short-word scalar entry remains class D.
Genuine shipping counterexample found: no.
Established shipping results retained: A* bridge, complete-word kernel <=1, persistent compatibility exclusion, qualitative superword strictness, compact outer release and finite practical absorbing radius.
Next calculation using more shipping structure: form the literal normalized S-correction factors from actual covariance-derived K_S together with generated R_S,T_S, then solve the complete 60/100-s same-history quadratic dissipation problem before any enclosure.

## Controlling quantitative route: adaptive shaped storage

The controlling next proof is now `docs/ou3-adaptive-shaped-storage.md`, not the short-word scalar chi/gamma route. In shipping-generated coordinates y=D(tau,sigma)(v,p,S,a_w), the frozen-tuple literal integrated OU transition is dimensionless and depends only on h/tau; tuner commits enter as exact causal diagonal coboundaries D_{k+1}D_k^-1, which telescope over a word and must not be independently normed. Literal S=0 corrections are transformed with their actual covariance-derived gains, while R_S and T_S remain generated from the same SpectralMSE/tau-cadence history. The next certificate seeks a same-history W_N-W_0 <= -Dissipation + Supply inequality over 60/100 s, carrying AG/BA/magnetic/accelerometer loss and rotating compatibility jointly. A common-Q failure is class D, not filter instability.

Structures preserved: literal OU primitives, S feedback, covariance/gains, causal tuner generation/staging, Mahony-driven tuner input, persistent history, source/loss direction, rotating weak directions.
Relaxations introduced: quadratic Q is only a candidate storage form; it does not enlarge the shipping execution family.

# Shipping-faithfulness gate — READ BEFORE NEW MATHEMATICS

The permanent governing protocol is [`docs/ou3-shipping-faithfulness-protocol.md`](ou3-shipping-faithfulness-protocol.md): **WE ARE PROVING STABILITY OF THE ACTUAL SHIPPING OU-III FILTER.** Before promoting any important result, record `Structures preserved` and `Relaxations introduced`. Classify every failed attempt A/B/C/D/E. PR #643's short-word scalar completed-square failure is **D. SUFFICIENT-BOUND FAILURE**, not filter instability. No proof-word boundary may reset a shipping/adaptive/physical history that shipping itself preserves.

## Adaptive OU/S normalization opened as the shaped-storage candidate

For one generated tuner tuple define dimensionless chain coordinates a_bar=a_w/sigma, v_bar=v/(sigma tau), p_bar=p/(sigma tau^2), S_bar=S/(sigma tau^3) and normalized time s=t/tau. The homogeneous continuous OU/integrator chain then has parameter-free drift a_bar'=-a_bar, v_bar'=a_bar, p_bar'=v_bar, S_bar'=p_bar; the exact discrete phi_pa and phi_Sa are precisely its sampled primitives. Thus tau and sigma should enter a long-word storage through the causal coordinate scaling, not as independent nuisance extrema. For the deployed SpectralMSE law, before explicit cadence clamps,
r_S=C_J q_eff^(1/14) sigma_aB^(6/7) tau^(24/7)/sqrt(T_S),
sigma=c_sigma sigma_aB and T_S=c_T tau, hence
r_S/(sigma tau^3)=(C_J q_eff^(1/14)/(c_sigma sqrt(c_T))) sigma_aB^(-1/7) tau^(-1/14).
The normalized S regularization therefore depends only on weak -1/7 amplitude and -1/14 time-scale powers along one generated adaptation history. This is the main analytical candidate for W_out. Clamp branches must be carried literally; no claim of a global constant normalized gain is made yet. Time-varying scaling adds causal adaptation-rate/coboundary terms, which must be bounded from the actual EMA and one-sample commit chronology.

Mahony scope correction: the private VerticalAccelComplementary Mahony observer is measurement-only. In Live it continues to generate the levelled vertical signal used by the wave-period/variance/tuner path; it does NOT directly overwrite the MEKF attitude every sample. Main nominal attitude evolves through literal gyro prediction, accepted accelerometer/magnetic corrections and the tilt watchdog/reset path. Therefore the long proof should use Mahony as a constraint on the generated tuner/reference coefficients and startup handoff, not invent a direct Mahony feedback term in the Live MEKF attitude dynamics.

## Next falsifiable horizon test: 16 s / 30 s / 100 s linked ratio

The short 0.32-s completed-square handoff is quantitatively dead, but that does not determine the superword result. The current diagnostic now exports one continuous carried physical/filter history over 16 s, 30 s and 100 s from the same regular A21 root, composes every local defect prospectively, preserves literal covariance/tuner/scheduler chronology, and evaluates the exact same-word M_W,b_W,J_0,J_N with per-history optimized gamma. The decisive finite diagnostic is inf_gamma chi_gamma/gamma at each horizon, together with operation-class terminal-metric attribution and raw state-block attribution. If the ratio collapses as signed forcing cancels and homogeneous information accumulates, proceed to source-uniform superword enclosure. If it remains orders of magnitude above .0225 at 100 s, stop trying to use V=e'P^-1 e as the sole outer storage and construct an outer storage W_outer with a proved finite handoff W_outer -> V inside the local region. This is a formulation decision, not permission to change MARINE MOTION, MAGNETIC SERVICE, SLOW+FAST constants, estimator behavior, or quality gates.

## Current limiter: 0.15 linked completed-square entry is numerically refuted in the present storage/word formulation

The existing carried moving/wave word already gives the decisive lower bound independently of the historical gamma=lambda_min/2 choice. In terminal-whitened coordinates its exported endpoint forcing has ||b||=9.589840717765055942946383 and the homogeneous loss cone has lambda_min(Delta)=0.0006929476516703682766284703. Since chi_gamma=b'b+z'G_gamma^{-1}z >= b'b and every admissible gamma<lambda_min(Delta),

    inf_gamma chi_gamma/gamma >= ||b||^2/lambda_min(Delta)
                              = 132715.717804109...,

whereas the requested 0.15 handoff requires <0.0225. The ratio to the budget is >5.898e6. Thus per-word gamma optimization cannot rescue this particular complete-word storage inequality. This is a NON-PROMOTING finite carried falsification of the proposed quantitative handoff, not a source-uniform instability theorem and not a refutation of qualitative strict dissipativity.

Failure classification: quantitative formulation failure. Invalidated hypothesis: that optimizing gamma in the current 0.32-s literal-word completed-square bound could plausibly close the 0.15 outer-to-inner entry. Retained facts: qualitative zero-action exclusion, compact outer A21 release, finite practical absorbing radius, exact completed-square identity, and the prospective local-defect construction remain valid. Current limiter: identify which same-history local source/coordinate class creates the large terminal-whitened b and whether the storage/word horizon is the wrong quantitative handoff object. Next falsifiable calculation: after native local-defect parity, transport every d_k to the endpoint, group contributions by prediction, acc, S, mag and reset while preserving within-group cancellation, and report both terminal-metric norms and raw state-block contributions. Do not interval-enclose the present chi/gamma formulation unless that diagnosis yields a mathematically different linked storage construction with carried margin.

## Local affine defect exporter: corrected boundary semantics

The earlier consecutive-pre-snapshot pairing was invalid because physical truth evolves across prediction. The proof-only observer now exports both sides of each literal mean operation. Prediction uses the previous physical sample as its pre truth and the current physical sample as post truth; accepted acc/S/mag corrections and reset use the same physical epoch on both sides. Covariance-only sync operations remain chronological but mean-neutral. The prospective composition remains d_k=e_{k+1}-A_k e_k and b_{k+1}=A_k b_k+d_k; endpoint residual is verification only. Native parity is still a required diagnostic before source attribution is trusted.

## Complete regular-word joint zero-action kernel CLOSED (nullity <=1)

The authoritative qualitative reduction is now formalized without attitude-only Schur compression. Zero fresh process action plus four distinct applied S rows (extended-Chebyshev system for {1,t,t^2,psi_tau}) kills the independent homogeneous LIN/AW root. Zero homogeneous magnetic loss together with MAGNETIC SERVICE reduces the remaining root attitude space to dimension <=1. On active A21 the literal accelerometer Jacobian has J_ba=I, so a pure BA kernel is impossible and, once an attitude amplitude is chosen, one accepted accelerometer row uniquely determines the complete BA root; additional rows can only remove that line. Therefore every complete regular word has joint zero-action kernel dimension <=1. If nonzero it is the word-dependent physical tilt/BA compatibility line; exact existence of a line is not asserted on every word. Its persistence through qualified MOVING superwords is separately excluded by the closed A* physical source theorem. Continuous magnetic service rows are covered by their SPD aggregate Gramian and regular prediction/reset transports preserve nullity, so no callback enumeration is required. literal_event_strata_regressed is now TRUE in this qualitative sense.\n\n## Correction: A* attitude block has nullity zero, not one

The proposed target nullity(C_theta,W)=1 is false for the aggregate attitude-only block under Corollary A*. A* strictness is precisely positive definiteness of the three attitude-column Gram, so the attitude row stack has rank 3 and nullity 0. The familiar one-dimensional field-axis tilt/BA compatibility line appears only in the JOINT accelerometer attitude+BA equations after allowing BA nuisance compensation; it is not an attitude-only kernel. Consequently the graph form [[C_theta,0],[A_W,I]] is not the physical compatibility construction when C_theta is the A* attitude block: with C_theta full rank that stack is full rank. Added a regression that fails the attitude-nullity-one claim under an A*-positive block. The literal kernel theorem must instead retain the joint acc/BA block through nuisance elimination, intersect it with magnetic/gyro/process/S constraints, and then prove that JOINT kernel has dimension <=1. Do not promote literal_event_strata_regressed from an attitude-only rank calculation.\n\n## Canonical kernel correction: compatibility is a word-dependent graph

There is no single fixed rational compatibility vector valid for all complete words. The authoritative line is r_W=(a,-A_W a), where A_W depends on the actual nominal-force rows, resets, held/active BA transport and coupled chronology. A global fixed-vector rank regression would therefore be synthetic. The exact canonical lemma is graph-valued: after magnetic service eliminates its two service coordinates and process/S constraints eliminate independent LIN complements, let C_theta,W be the remaining attitude constraint and let A_W map attitude into forced BA. If nullity(C_theta,W)=1, then the exact stack [[C_theta,W,0],[A_W,I]] has nullity one for EVERY A_W and ker={ (a,-A_W a): a in ker C_theta,W }. Added exact rational graph_kernel_certificate and regression with arbitrary A. Thus BA graph coupling itself needs no interval cover; the remaining substantive statement is one-dimensionality of the reduced attitude constraint, followed separately by A* exclusion of persistent nonzero a.\n\n## Continuous MAGNETIC-SERVICE/transport cover collapsed analytically

For the qualitative kernel theorem no finite magnetic row cover is needed. MAGNETIC SERVICE states that on every 1-s window the aggregate transported/whitened heading--axial-BG Gramian satisfies sum G_i'G_i >= mu_M I_2 with mu_M=1. Hence intersection_i ker G_i={0} in the two service coordinates for arbitrary admitted event count, timing and row orientation. This is a complete continuous magnetic-row kernel cover. Likewise regular prediction and reset transports cannot create kernel dimension: OU/attitude prediction is invertible on the regular domain (phi>0 and triangular nonzero diagonal), and G_reset=I+[d]/2 has det=1+|d|^2/4>0. Thus dt/tau/rotation/reset cells collapse by invertible coordinate transport. The only nontrivial exact rank regression left is the canonical transported intersection of S/LIN, aggregate accelerometer geometry, gyro/process and BA compatibility rows after magnetic service coordinates have been eliminated. This is one canonical algebraic stratum, not an event-schedule enumeration.\n\n## Literal nullspace stratum assembly: covariance blocker removed, continuous cover still open

Added nullspace_strata.py. For qualitative zero action, correction gains and innovation covariance are not needed: exact Joseph zero loss implies H e=0. The literal stack therefore consumes only process/sync homogeneous rows, S rows, applied magnetic Jacobian rows, accelerometer Jacobian rows, gyro/process rows and BA compatibility rows, all transported to common root coordinates. Exact rank/nullity is then checked by complete_word_nullspace.py. This removes the open full-P seed blocker from the qualitative kernel lemma. However the current MAGNETIC SERVICE contract has no finite callback-pattern cover, and the repository still lacks a theorem-domain continuous cover of the admissible applied magnetic/service-row transports. Therefore literal_event_strata_regressed remains FALSE. The remaining task is specifically to construct that continuous service-row/transport cover (not enumerate schedules), then feed each interval/rational stratum to the exact kernel certificate.\n\n## Literal nullspace strata: scheduler enumeration replaced by service-row parameterization

The existing magnetic seed audit proves there is no finite callback-pattern cover under MAGNETIC SERVICE: accepted information is constrained in aggregate and regular callbacks may occur at arbitrarily many 4--6 ms patterns. Therefore "enumerate every reachable event stratum" cannot legitimately mean enumerate scheduler bit patterns. For the qualitative zero-action lemma gains/full P are also unnecessary: exact Joseph zero loss implies H e=0 at an applied correction. The literal nullspace stack can therefore be parameterized by the actual applied measurement Jacobian/service rows and continuous prediction/reset transports, avoiding the open full-P seed blocker. The exact-rational kernel lemma now has a service-row stratum interface and explicitly refuses callback-pattern enumeration. What remains is a theorem-domain cover of the continuous service-row/factor parameters, not a finite schedule list.\n\n## Exact block-nullspace lemma formalized

Added complete_word_nullspace.py. For the authoritative literal zero-action row blocks C_proc/S, C_mag, C_acc, C_gyro/process and C_BA, the lemma is exact: if a supplied nonzero compatibility vector r satisfies Cr=0 and exact rational elimination gives rank(C)=n-1, then ker C=span(r). No SVD tolerance, independent observability lemma or four-S=>AW=0 shortcut is used. Regression includes a deliberate hidden-second-kernel case and fails closed. The remaining step is literal event-stratum assembly: export the actual complete-word homogeneous zero-action matrices into this stack and check rank/nullity on every reachable regular stratum. A* persistence exclusion remains a separate physical theorem applied only after the literal kernel equality.\n\n## Local affine defect exporter in progress

The retrospective word defect is being replaced by literal boundary composition. The temporary proof header now snapshots xext/qref at every internal prediction, accepted acc/S/mag correction and reset; sync/sync-completion are classified mean-neutral. Consecutive mean-operation pre-snapshots therefore provide exact post/pre boundaries, including quaternion-reference changes across resets, at one unchanged physical sample time. For each boundary d_k=e_{k+1}-A_k e_k and b_{k+1}=A_k b_k+d_k. The endpoint residual e_N-M e_0 is reserved only for verification. Algebra and boundary-pairing regressions are committed; native execution/parity must pass before the locally composed b_W is used in linked chi/gamma.\n\n## Linked entry ratio: correct optimization and scope

For one literal word the relevant feasibility value is inf_{0<gamma<lambda_min(Delta)} chi_gamma/gamma, not the historical choice gamma=lambda_min(Delta)/2. The carried linked diagnostic now performs the one-dimensional optimization on its exact same-word M,b,J0,JN. This is still retrospective because b=e_N-M e_0 is reconstructed from the realized endpoint; it is not a source-uniform endpoint-defect map valid under varying root error. No carried pass/fail is promoted until b is generated from the same physical SLOW+FAST/process/nonlinear source history. The theorem target remains sup_history inf_gamma chi_gamma/gamma < .0225 with gamma chosen measurably/continuously inside the same word's homogeneous loss cone, plus every-prefix retention.\n\n## Outer-to-inner entry: exact status after qualitative strictness

The compact outer release region now exists at regular A21 roots: C_out=sup_Rrelease V<infinity. Qualitative finite-superword homogeneous strictness gives, on every compact annulus r_in^2<=V<=C_out, a positive homogeneous loss margin gamma_ann>0 after choosing a sufficiently long common superword. This is NOT yet finite entry under physical forcing. For the same reachable word the exact completed-square identity requires chi_gamma/gamma<r_in^2, with r_in<.15, and every prefix needs its own linked bound. The candidate SLOW+FAST source theorem proves the geometry that makes gamma_ann positive but does not numerically bound the linked chi_gamma supremum. Compactness gives only finite Chi_ann=sup chi_gamma, not Chi_ann/gamma_ann<.0225. Therefore outer retention/finite entry into the local ball remains open at one precise linked finite-error inequality; do not infer it from LaSalle strictness alone.\n\n## Outer A21 release region: existence closed at regular roots

Captured-domain H18 release is compact in mean/tuner/covariance/scheduler/reference coordinates, and held-H18 LIN BIBO is source-uniform. After the first complete 16-s regular A21 LIN window, the recurring root covariance certificate gives a uniform Loewner lower P_root>=L_root>0, hence J_root<=L_root^{-1}. Therefore V=e'J e is continuous and uniformly bounded on the compact reachable release image at those roots. Define C_out=sup_Rrelease V<infinity and R_outer={reachable regular A21 roots with V<=C_out plus the compact carried auxiliary state}. This is an existence-level outer region; no numerical C_out is asserted. BA is eliminated only through the correct marginal P_oo when an outer coordinate reduction is used.\n\n## Exact complete-word nullspace intersection: corrected formulation

Do NOT use the retracted statement that four S atoms force base AW=0. The authoritative q=0 classification is: complete-word zero action iff one literal homogeneous trajectory has zero fresh process/sync action, zero S action, magnetic-axis AG compatibility and zero accelerometer compatibility action at every event. Estimator-level forward-compatible ZG trajectories can satisfy these equations locally; four-S alone does not exclude them. The new ingredient is external to that algebra but same-history: the candidate SLOW+FAST + bounded-velocity theorem closes the real-arithmetic source transfer and proves the Corollary-A* aggregate nominal force cannot remain in the field-compatible zero-action geometry on every required moving window. Therefore a shipping-closed physical ZG execution persisting through a complete qualified moving superword is excluded. Combining this with MAGNETIC SERVICE leaves only the declared gauge/tilt-BA compatibility line on finite prefixes; persistence of that line through every qualified moving superword is excluded by the same A* physical source margin. This is the proper LaSalle invariant-set exclusion. It does not assert four-S pointwise AW nullity or the stronger numerical G0 floor.\n\n## Coupled zero-action kernel condition after A* closure

A* alone does not prove the nuisance-projected VF-A modulus: root AW/BA columns can mimic accelerometer attitude rows in a Schur complement. The valid zero-action argument must impose all complete-word zero conditions simultaneously. In a zero-action limit, four-S plus LIN/AW process action forces the independent LIN/AW homogeneous component into its established null (zero on the regular four-S block); BA process/projection leaves only the held compatibility component; MAGNETIC SERVICE restricts AG to the transported field-compatible attitude/BG line; accelerometer zero action then sees the A* aggregate attitude block with no independent AW nuisance left to mimic it. Hence any surviving zero-action vector lies in the established field-axis tilt/BA compatibility line, not an additional normal direction. This is a qualitative kernel statement; it does not assert a numerical Schur modulus. The remaining proof obligation is to write this intersection as an exact block-nullspace lemma against the literal complete-word matrices and regress its rank on all event strata.\n\n## Complete linked matrix revisited after A* source closure

The historical direct-information attempt failed at VF-A because no theorem then separated the estimator nominal accelerometer coefficient from the magnetic-compatible nuisance span. That premise is no longer absent: the new source-uniform real-arithmetic AW transfer proves Corollary A*'s nominal signed-mean separation under the candidate MARINE/SLOW+FAST qualification. This does not localize to same-cell magnetic epochs, but it does apply to the complete aggregate accelerometer block, exactly where the linked reduced-matrix/complete-word formulation needs it. Therefore the same-cell c,b0 route is retired and the complete-word zero-action nullspace argument is reopened with A* as an established aggregate premise. The next target is qualitative: any normalized zero-action sequence must have magnetic-compatible AG direction, zero four-S/LIN-AW action, zero accelerometer aggregate attitude action, and zero gyro/process action; A* excludes a nonzero attitude component transverse to the field, MAGNETIC SERVICE excludes heading/axial-BG complements, and four-S/process kills independent LIN/AW complements. The only allowed common null is the already identified physical tilt/BA compatibility line. A quantitative rate is still separate.\n\n## Same-cell c,b0 route fails at c localization

The proposed exact E=0 same-prediction-cell acc->mag construction cannot be made source-uniform from the current premises. MAGNETIC SERVICE guarantees informative applied magnetic rows on every 1-s interval, but it does not require those rows to occur at accelerometer cells whose nominal force is noncollinear with the reset-pulled magnetic direction. The newly closed Corollary-A* premise is an aggregate signed-mean statement; it proves a positive attitude-column Gram over the window, not a pointwise/same-cell angle at the magnetic epochs. Therefore no uniform c>0 follows for the same-cell groups. Per the failure protocol this route stops here; do not infer c from average force separation. The valid continuation is the aggregate linked six-column reduced matrix, which retains all accelerometer rows, actual magnetic service rows and chronological gyro transport jointly.\n\n## Real-arithmetic source-only AW transfer CLOSED for Corollary A* positivity\n\nUnder the user-qualified MARINE profile and candidate SLOW+FAST temporal caps, the conservative source-only charge is 1.07943484027807165 m/s^2 < g/5=1.96133, leaving .88189515972192835 m/s^2 before float32. The finite-error reset terms d^2 v/6+c(r)v^2 are downstream nonlinear supply, not source. The normalized quaternion polynomial is already included in the qualified gyro prediction transport angle/singular-floor certificate and is not charged a second time. Therefore the real-arithmetic source-uniform Corollary-A* signed-mean premise is CLOSED. This establishes positive attitude-column Gram pointwise; compactness of the qualified regular word class yields existence of a uniform positive attitude-Gram minimum. It does NOT promote the historical explicit G0 floor, whose stronger m_perp<=.4 and u1<=1.2 premises remain unproved.\n\n## Corollary A* versus explicit G0 constants

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

## User-qualified MARINE excitation constants (2026-10-01)\n\nFor the controlling theorem, use T_E=T_P=30 s, theta_E=2 deg, P_E=.03 m. These are theorem qualification premises supplied by the user; they are not inferred from RAO data. IMU fast H_a,C_a,H_g,C_g remain OPEN. Consequences: P_E/T_P=.001 m/s; theta_E/T_E=.00116355283466 rad/s. The generic signed acceleration chord lower max(0,P_E-T_P V_max) remains zero, so no pointwise acceleration floor follows.\n\n# OU-III proof: SLOW + FAST physical IMU qualification

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

## CI: main after PR #643 (Ruff, readout lift, stale provenance)

Failed quantities: `quality-gates / python` (59 Ruff findings, including
`F821 _out` undefined in `linked_guard_aw_reader.py`, so that module raised
`NameError` on every call); `ou3-stability-proof / theorem` at
`ag_readout_source_diagnostic` ("event 4 missing physical lift fields"); and
in the full evidence suite, stale provenance/status, the same-cell grouper
("unretained hard event in regular AG group"), the world-frame observer anchor,
and stale test expectations. Classification: E (implementation/CI) for every
item; nothing here is a mathematical result or a shipping counterexample.
Invalidated hypotheses: every native trace event is an operation boundary
(the driver now also records a read-only `adaptive_state` tuner snapshot);
the world-frame quaternion tap can anchor on the pre-#643 reset line.
Fixes: `_out` rounds positive radii outward and keeps exact zeros; the lift
check and same-cell grouper pass `adaptive_state` through like `sync`; the
quaternion tap follows the relocated reset tap (qref is identical at both
points); diagnostics, `theorem-status.json` and provenance were regenerated
from the current sources. Retained facts: all regenerated margins keep their
sign (signed-balance S-interval margin -19.0391078, rotation-triangle margin
-30.4543160; same-cell group values move by about 1e-9 with shifted event
indices). Current limiter and next falsifiable experiment are unchanged.
