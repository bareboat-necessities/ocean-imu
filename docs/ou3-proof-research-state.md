# OU-III proof: current continuation state

Source audit: PR #637 at `1bbd3b3bc58879d7f95863cef5726ee10c80616f`.
The complete preceding ledger is preserved byte-for-byte in
[the pre-field-axis-audit ledger](ou3-proof-research-state-before-field-axis-audit.md).
The long-form derivations remain in `ou3-corrected-word-proof.md` and
`ou3-linked-soft-return.md`; do not revive their explicitly retracted claims.

## Single proof path and unchanged scope

construction -> capture -> magnetically informed H18 -> refinement/release ->
recurring A21 -> regional practical stability. The contraction word remains
100 s; 17 s remains nuisance/root warm-up. Runtime, physical assumptions,
coupled tau/sigma_aw/r_S/T_S chronology, noise, and quality gates are unchanged.
Stillness is not removed from the physical class by this continuation.

## New analytical result

[Source-audited regional field-axis exclusion](ou3-field-axis-regional-exclusion.md)
(FA1--FA12) is the current calculation. It repairs two source mismatches in QR:
accelerometer updates occur every regular IMU sample (maximum qualified applied
gap 6 ms, not 40 ms); default AW synchronization is an additive PSD floor, not
congruence. Isotropy makes its marginal eigenvalues exactly max(old,target),
which supplies the needed ceiling without a full 21-state covariance ceiling.
An exact-rational audit of the literal Q polynomials supports the conservative
regular real-arithmetic bound P_aw,aw<=16.48 I, hence ||e_aw||<4.06 sqrt(V_base).

The sharp same-history trapezoidal inequality gives

    sum_i w_i ||P_bi(a_hat_i-g_model)||^2 >= max(M_ref(r),0)^2,
    M_ref(r)=G0-delta_g-eta_ref(Amax+g_model_max)
             -4.06(1+eta_ref)r-2Vmax/L-Jmax h_acc/4.

Every a_hat_i and V_base,i is evaluated at the actual PRE-ACCELEROMETER prefix.
For a fixed, correctly bounded nominal reference with inclination <=80 deg,
the 100-s pre-radius margin is 1.4429069015... m/s^2. At r=.25 the margin exceeds
.4279 m/s^2 before reference-error charges (.4278 after worst-case word-boundary trimming). This is a conditional local
exclusion of the particular field-aligned base trajectory, not stability.

## Scope correction / unresolved implication

The physical magnetic field and the learned nominal reference are different
objects. MagAutoTuner's initial 5% horizontal-fraction gate gives a weaker
positive fixed-reference margin, but the default continuous hard-iron path
rewrites the reference without reapplying that fraction gate. This does NOT
prove it fails in a shipping execution. It prevents treating the startup cone
or physical 80-degree cone as an already established all-time nominal cone.

Next: prove eta_ref (or a sufficient direct nominal cone/variation estimate)
from literal acquisition, refinement and hard-iron updates, then discharge
prefix retention in the same physical history. Do not invent zero defects,
freeze the default reference, disable continuous hard iron, replace actual
accepted-update times, or assume stability while proving retention.

## Classification and validation

PROVED CONDITIONAL: FA1--FA12 under their stated regular covariance, prefix,
reference and physical premises. The new bound is a base geometric row-energy
floor, not the full innovation-weighted nuisance-eliminated word action.

OPEN: unconditional shipping exclusion, capture/release, source-uniform linked
O2 return bound, numerical full-word rho<1, nonlinear/physical supply retention,
and target arithmetic. No theorem or quality flag is promoted. Zero homogeneous
action must never be substituted for zero BASE innovation. Regional r is a
full-state covariance-metric radius, not the retained six-degree attitude angle.

The 13 new exact scalar/rational regressions pass with
`python -m unittest discover -s tests/validation -p 'test_ou3_field_alignment_exclusion.py' -v`.
These validate the analytical substitutions, not a full trajectory enclosure.
Full native validation and CI were not run here. Inherited failures, numerical
experiments, unsuccessful routes, and exact previous status are retained in the
archived ledger and previous PR history; they are not claimed repaired.


## Shipping-reference and prefix-premise audit — 1c8f4e94

The requested follow-up is now source-audited in
`docs/ou3-shipping-reference-prefix-retention.md`.

The continuous hard-iron reference update has an exact same-statistics
Lipschitz law: with `L(b)=wbar-Abar*b` and `||Abar||<=1`, each canonical
horizontal/z reference component moves by at most the applied body-bias
increment, and the full canonical reference by at most `sqrt(2)` times it.
The 45-s slew gives the corresponding per-update increment. This is literal
shipping chronology, not a compatibility relaxation.

That result does NOT yield the needed all-time nominal cone from the present
contract. Startup/refinement enforce a 5% horizontal fraction, but the default
continuous path does not reapply that fraction gate; it only requires positive
horizontal magnitude above .001 uT. Its loose accepted-fit envelope
`.35*(75+5+2)=28.7 uT` exceeds the declared physical 15-uT horizontal
minimum. More decisively, its statistics use the private Mahony tilt, and no
all-time deterministic true-to-proxy tilt tube has yet been proved. Therefore
the physical-field cone cannot be silently transferred to the nominal
reference.

The FA12 full-storage premise also cannot be inherited from six-degree capture.
The fixed-reference field-exclusion radius must satisfy `r<.355396`, whereas
the handoff tilt covariance sigma is .035 rad, so a six-degree tilt error alone
has minimum covariance-metric radius `(pi/30)/.035=2.99199...`. This does not
say handoff has six-degree error; it proves the existing capture target does
not imply the small FA12 storage ball.

To remove that artificial circularity, the field-axis lemma is reformulated in
the exact component it needs. If
`||a_hat_w-a_phys||<=eps_aw` at every relevant prefix, the fixed-reference
100-s exclusion needs only
`eps_aw < g*sin(10deg)-.26 = 1.4429069015... m/s^2`.
The varying-reference version is given as PR10 in the new note. This is much
weaker than requiring the entire 21-state Mahalanobis error to be <=.25.

Prefix invariance remains a simultaneous fixed-point problem, not an upstream
premise: the same block must close kernel covariance return, word error supply,
and every-prefix supply. Proving `V<=r^2` first from covariance bounds would
be circular because the strict word return needed for that storage recursion is
the open O2/BP obligation.

Next decisive calculation: derive (or falsify) an all-time private-Mahony tilt
tube and an all-prefix AW tracking tube from their literal error equations on
the SAME physical history, retaining the coupled tau/sigma_aw/R_S/T_S
chronology; insert those component tubes into BP-4/BP-10 and solve the retained
rectangle simultaneously. No theorem flag is promoted.


## Correction: the 1.4429069-m/s^2 target is sufficient, not yet implied — current head

The immediate blocker is still ONLY the persistent nominal force/field
collinearity trajectory. The desired contradiction is

    |weighted mean P_B(a_hat_w-a_phys)| < 1.4429069015... m/s^2

for the 100-s, 6-ms, 80-degree specialization, because physical bounded
velocity+jerk contributes at most .26 m/s^2 while the transverse gravity
requirement is g sin(10 deg)=1.7029069015....

However, the existing sections 59--61 of
`docs/ou3-aw-adjoint-cancellation.md` already prove that this inequality
CANNOT be inferred merely by saying that a_hat_w, tau, sigma_aw, R_S and T_S
are coupled.

For a prescribed smooth periodic physical history, the measurement-only
front end and tuner determine one periodic applied coefficient word U_*.
Conditional on U_*, covariance propagation is independent of innovation
VALUES. After the stabilizing periodic Riccati orbit is fixed, the mean is an
affine periodic linear system. The finite Kalman map from a measurement word
to its innovation word is block lower triangular with identity diagonal and is
therefore invertible. The exact periodic compatibility equation is PR9/PR14.
Generically nonsingularity gives a UNIQUE compatible forced periodic orbit; it
does not exclude one.

Therefore the coupled tuning law is essential for fixing the coefficients, but
it is not an amplitude theorem for deterministic estimation error. Covariance
P_aw,aw<=16.48 I likewise bounds uncertainty/action geometry, not
|a_hat_w-a_phys| on an arbitrary deterministic forced execution.

This invalidates the proposed shortcut

    coupled tuner law => |mean(a_hat_w-a_phys)|<1.4429.

The 1.4429 number remains a correct sufficient threshold. To prove the
pathological trajectory inadmissible one must instead show that the UNIQUE
self-consistent PR14 solution violates an EXISTING admissibility condition
(physical p/v/a/jerk, bias/projection, gate, service, retained local angle), or
prove a sharper signed identity that forces such a violation.

The constructive audit has already made this falsifiable. The first
commensurate 6-s candidate is excluded because its literal late relative
attitude reaches 7.488 degrees >6 degrees. The analytically refined 12-s
candidate has pre-compensation physical amplitudes strictly inside the
declared envelopes and is the decisive next target. It has not yet been
certified as a shipping counterexample because its late periodic
tuner/covariance/mean orbit, six-degree bound, gates/service and complete
weighted functional remain to be enclosed.

Accordingly, do NOT claim the accel||mag pathology is excluded yet, and do NOT
spend the next calculation deriving a generic AW tracking tube from covariance.
The decisive calculation is the literal 12-s PR14 periodic orbit with outward
enclosure. If it violates an existing condition, extract that violation as the
analytical exclusion lemma. If it satisfies all conditions, the pathological
trajectory is admissible under the current theorem contract and the stability
proof must be reformulated; no additional physical assumption may be silently
introduced.


## Controlling LaSalle update — radius-local field-axis candidate excluded

The controlling local proof is now `docs/ou3-radius-local-field-alignment.md`.
It replaces the attempted global nominal-AW/innovation control for the
zero-dissipation invariant-set question.

Inside the retained storage ball, the current literal source audit gives
`||a_hat_w-a_phys||<4.06 r` from `P_aw,aw<=16.48 I`; the inherited BA
marginal gives `||e_ba||<=r/40`. BA is NOT added to the field-alignment
tube because the implemented accelerometer attitude Jacobian uses the nominal
CoG vector `a_hat_w-g`; BA is a separate measurement column and lever arm is
attitude-independent there.

On the explicit committed-field branch `sigma_w>=1/5`,
`||P_B g||>=9.80665/5=1.96133 m/s2`. For one MARINE history,
`||v||<=5.5` gives a continuous T-window point with
`||P_B(a-g)||>=1.96133-11/T`. The 100-m/s3 jerk bound and an ACTUALLY
APPLIED accelerometer gap <=.006 s transfer this to a sampled epoch with

    m_phys(T)=1.96133-11/T-.6.

At T=17 s, `m_phys=0.714271176470588...`. Hence persistent nominal
field alignment is impossible whenever
`4.06 r<m_phys`, i.e. `r<0.175928...`. Adopt the deliberately
conservative local exclusion radius

    r_FA=.15, T_FA=17 s,

with strict margin `0.105271176470588... m/s2`.

Conditional on the already established zero-dissipation classification
(leaving only the field-axis candidate), this proves

    Inv_MARINE({D=0}) intersect {V<=.15^2} = {0}

on the regular real-arithmetic retained branch. The existing closed-stratum
compactness argument then gives existence of finite `m` and `eta_D>0`
for homogeneous finite-window strict dissipation by contradiction/diagonal
extraction. No numerical eta_D is claimed.

OPEN after this local invariant-set closure: entry/every-prefix retention in
V<=.15^2, finite capture/H18/release, nonlinear/source supply over the
finite-window block, recurring transition budget, and full float32 totality.
The end-to-end regional practical-stability theorem is NOT promoted.

Historical O1/O2/kernel-ceiling and signed-reader calculations remain useful
research but are non-controlling for exclusion of the zero-dissipation
field-axis trajectory.


## Entry/retention continuation — staged A21 entry required

`docs/ou3-staged-entry-retention.md` is now controlling downstream of the
r_FA=.15 field-axis exclusion.

Direct H18 release into full V<=.15^2 is NOT a valid universal target.
Shipping held BA is decoupled with sigma_bacc0=.004 m/s2; enabling A21 merely
floors its diagonal variance to sigma_bacc0^2. The admitted physical BA norm
is .22516660498395405. A permitted held b_hat_a=0 therefore has decoupled
release contribution V_ba=(.22516660498395405/.004)^2>3168, versus .0225 for
the final local ball. No assumption or runtime change is made.

Homogeneous every-prefix retention IS closed once inside the local ball:
prediction/correction are covariance-metric nonexpansive and congruent resets
preserve storage; PSD covariance inflations cannot increase fixed-error
storage. Finite nonlinear/source/arithmetic residuals require an inner root
radius r_in<.15 with

    sup_k G_k(.15,d) <= .15-r_in,
    E_W(.15,d) <= (1-q) r_in,

q=sqrt(1-eta_D).

The correct entry path eliminates the held BA coordinate. At release it is
literally decoupled; after release use the BA Schur complement / quotient
storage V_o|ba. Prove quotient-local entry/retention at .15, apply the
field-axis LaSalle exclusion there, then use active A21 BA dynamics and the
dissipative projection sector to enter the full ball. This is the same
H18-complement architecture, not a new proof path.

Current limiter: source-uniform outer/BA-quotient release set inclusion and
quantitative quotient finite-window dissipation/supply. Do not retry universal
direct full-V release.


## Literal release outer-storage audit — attitude is the entry obstruction

A dedicated unchanged-header release snapshot now evaluates the carried
construction history at first A21 activation (step 36008, 180.039996 s).
Correct BA elimination uses the outer covariance marginal; at the release
boundary BA cross covariance is zero so it also equals the conditional outer
block there.

Release tilt is 8.144927 deg. The attitude covariance eigenvalues are
3.5352043e-6, 3.5550410e-6 and 1.2944541e-5. Minimizing the BA-eliminated
outer storage over every other non-attitude outer coordinate still gives

    V_outer,elim >= theta' P_theta^-1 theta
                  = 5716.295063726523,
    sqrt(V_outer,elim) >= 75.60618403098071.

Thus this literal release is nowhere near r_FA=.15; ATTITUDE alone consumes
the radius. At this covariance, V<=.15^2 would necessarily require tilt
<=.000539678 rad=.0309213 deg.

This is finite carried evidence only: the stress history does not certify
all-time MAGNETIC SERVICE and therefore is not an eventual-capture
counterexample. It does prove that stage flags/release mechanics themselves
do not imply entry into the tiny local storage ball. The controlling proof now
needs an outer A21 retained/capture region and finite entrance from that region
to the inner r_FA=.15 LaSalle ball. The r_FA field-axis inequality cannot
simply be enlarged to the observed release storage: its current source-audited
AW component conversion loses positivity above r~=.176.


## Literal release outer-entry calculation — AW is the blocker

The requested BA-eliminated release calculation is recorded in
`docs/ou3-release-outer-entry-audit.md`.

A dedicated unchanged-header diagnostic snapshots first A21 BA activation at
step 24016 (120.079997316 s) on the carried diagonal-wave construction. Using
the exact physical construction and outer marginal identity
`V_o=min_ba V=e_o'P_oo^-1 e_o`, the finite carried release has
`V_o=79300.039948`, sqrt=281.603, versus target .15.

Principal-block lower bounds after minimizing every other coordinate identify
AW as dominant: attitude 9.7976, bg 3.59e-7, v 331.982, p 1996.290,
S 84.977, AW 60867.918. The release AW error norm is about 8.385 m/s2.
This finite history is not an all-time service/capture certificate.

The second requested calculation has a negative but decisive result:
the present assumptions do NOT imply source-uniform pointwise AW entry.
An existing admitted A21 history already refutes pointwise physical-AW
tracking at 7.647 m/s2 on a 16-s window, whereas the r_FA=.15 local lemma
requires <.609 m/s2. P_aw,aw<=16.48 I only converts an ALREADY SMALL storage
to an AW component bound; it does not bound deterministic AW error before
entry. No AW mean projection supplies such a cap.

Therefore direct H18/release -> V_o<=.15^2 is not the correct bootstrap.
The local FA theorem remains valid after entry, but entry must use a shaped
signed/windowed AW functional or return to the global same-history FA
reachability/action calculation. Do not retry a source-uniform pointwise AW
tracking lemma; it is already falsified on the admitted class.


## Outer-entry no-go and controlling reformulation — current

The requested universal shaped-storage theorem `every certified H18/A21 release -> eventual retained V<=.15^2 entry` is analytically false under the present contract. A constant-rest member of the exact stationary attitude/BA gauge suffices: rotate the true attitude by alpha=1/1000 rad about B, set b_g=0 and b_a=g(Q_alpha e_z-e_z), and keep p=v=a=0. It has zero bias rates, ||b_a||<.009807 m/s2, exactly nominal accelerometer/gyro/magnetometer packets forever, and retains the quiet actual MAGNETIC SERVICE floor. Arbitrary complete stillness is explicitly admitted, so moving excitation is not owed. On the identical nominal filter history the proved quiet P_ba,ba<=I/1600 gives sqrt(V)>=||e_ba||/.025>.392>.15.

FAILED INEQUALITY: no source-uniform `W_out large => Delta W_out<=-epsilon` can hold toward the point physical-error set on every certified continuation, and no universal eventual retained point-entry time exists. Failure class: identifiability/theorem-target failure, not conditioning or numerical sharpness. Invalidated hypothesis: the full coupled shipping structure plus MAGNETIC SERVICE is sufficient to collapse stationary attitude/BA ambiguity. Retained facts: the r_FA=.15 local LaSalle exclusion after point entry, the 16.48 AW covariance ceiling, operation-wise homogeneous nonexpansion, correct BA storage elimination, and release audits remain valid.

CONTROLLING PATH: one nested architecture, but the outer set is distance to the stationary measurement-compatible attitude/BA class. Prove stationary practical retention to that class; then use an actually complete MARINE moving window to collapse the gauge and enter V<=r_in^2<.15^2 with retention; only then invoke the local LaSalle theorem. Indefinite physical rest has a consistency-class conclusion rather than an impossible point-error conclusion.

NEXT FALSIFIABLE CALCULATION: form the exact stationary gauge tangent K_stat at the carried A21 root, quotient the complete corrected-word action by K_stat, and evaluate the remaining gauge action over one complete same-history T_E moving window using the literal coupled tau/sigma_aw/R_S/T_S chronology. Do not refine a global shaped point-storage or retry pointwise AW tracking.


## Mahony/gauge quotient audit — current

Mahony is part of the literal coefficient chronology, not an independent observation. The shipping order conditions accelerometer input, advances the private measurement-only Mahony vertical observer, updates period/sigma tuning, stages the coupled online tuple, and commits it at the next IMU sample. Therefore any corrected-word quotient calculation must carry Mahony/front-end state and the resulting lagged coupled tau/sigma_aw/R_S/T_S sequence. Identical conditioned IMU histories imply identical Mahony and tuner histories.

At a stationary root, with body magnetic vector b=Q^T B, the acc/mag physical observation differential has the one-dimensional attitude/BA kernel

    K_stat = span{ (delta_theta=b, delta_ba=g[Q^T e_z]_x b) },

with all other error coordinates zero (up to the global attitude-error sign convention). For complete-word quadratic action J_W, quotient by choosing a root-metric complement Z and using Jbar_W=Z^T J_W Z, equivalently minimize the action over additions lambda*k_stat.

The hoped-for strictly positive action on the remaining gauge coordinate over one generic moving T_E window is NOT implied by the current contract. The numerical T_E and theta_E fields remain OPEN/null in constants.json. More strongly, the existing exact sin^3 rest/motion witness has positive gravity-direction span on each complete moving window while its accelerometer, gyro and magnetometer packets remain exactly nominal through compensating admissible physical biases. Hence Mahony and the complete coupled tuner word are also nominal. Symbolic positive attitude span alone therefore cannot give a positive quotient floor.

Current limiter: derive the largest gravity-direction span Theta_gauge(T_E) achievable by this exact packet-indistinguishable family under the existing B_a,D_a,B_g,D_g and Omega_max bounds. Only an independently certified MARINE pair satisfying theta_E>Theta_gauge(T_E) could exclude this gauge and justify a point-entry moving theorem. Otherwise MOVING also requires a consistency-class theorem. No theorem flag is promoted.


## Exact hidden-gauge span envelope — current

For the packet-indistinguishable family rotate the physical attitude by phi(t) about the (fixed world) magnetic axis and compensate the physical residual biases so that the measured gyro, accelerometer and magnetometer packets equal the nominal packets. The exact relations are

    ||b_a|| = 2 g |sin(phi/2)|,
    ||dot b_a|| = g |dot phi|,
    ||b_g|| = |dot phi|,
    ||dot b_g|| = |ddot phi|,
    ||omega|| = |dot phi|.

Therefore every such history satisfying the existing bounds obeys

    |phi| <= A_g := 2 asin(B_a/(2g)),
    |dot phi| <= L_g := min(D_a/g, B_g, Omega_max).

D_g constrains curvature but cannot improve the source-uniform range bound on an arbitrary interior T-window, because constant dot-phi is admissible and has dot-b_g=0. Hence the exact sharp envelope implied by these five scalar bounds is

    Theta_gauge(T) = min(2 A_g, L_g T)
                   = min(4 asin(B_a/(2g)),
                         T min(D_a/g,B_g,Omega_max)).

With current constants g=9.80665, B_a=0.22516660498395405, D_a=.001, B_g=.02, D_g=1e-5 and Omega_max=.6108652381980153,

    A_g = .0229611081599661 rad = 1.31557459051 deg,
    2 A_g = .0459222163199322 rad = 2.63114918102 deg,
    L_g = D_a/g = .000101971621297793 rad/s
        = .00584254353047 deg/s,
    T_sat = 2 A_g/L_g = 450.343102674 s.

Thus

    Theta_gauge(T) = min(.0459222163199322,
                         .000101971621297793 T) rad.

This bound is sharp for arbitrary interior windows under the listed scalar constraints: a constant-rate segment realizes the Lipschitz branch (with D_g charge zero), and sufficiently slow ramps plus a plateau approach the amplitude branch while respecting D_g. Join smoothness may reduce a particular boundary-crossing construction, but the MARINE excitation quantifier applies to every complete window contained in a moving episode and cannot assume a rest join at each window endpoint.

DECISIVE CONDITION: a numerically certified MARINE pair can exclude the exact attitude/BA packet gauge only if

    theta_E > Theta_gauge(T_E).

Equality is not enough because the excitation premise is >= theta_E. The current constants.json still has T_E and theta_E null/OPEN, so the comparison cannot yet be discharged. Mahony does not alter the envelope: packet equality makes its measurement-only trajectory and the complete staged coupled tuner chronology identical to nominal.

Next: obtain/derive the existing theorem-grade MARINE (T_E,theta_E) from admissible physical evidence without strengthening the assumption. If none is currently certified, point-entry on MOVING remains conditional on the displayed strict inequality; proceed with quotient-action positivity only after it is satisfied.


## MARINE excitation qualification audit — current

A complete repository audit found no existing theorem-grade numerical pair (T_E,theta_E) to populate constants.json. The controlling proof documents intentionally keep both symbolic, and constants.json marks numerical qualification OPEN. The pinned v1.2.1 28-ft vessel-RAO bundles are finite statistical replay evidence with provenance; they do not certify an all-time rolling minimum of gravity-direction span for every complete window of every admitted moving continuation.

Nor can a positive theta_E be derived from the other present MARINE/IMU bounds. The exact packet-indistinguishable family admits alpha>0 arbitrarily small and nu>0 sufficiently small while satisfying B_a,D_a,B_g,D_g,Omega_max, zero translation/jerk/primitive and magnetic service. Its complete-window span is 2 alpha>0 but tends to zero with alpha. Thus the infimum of admissible moving-window gravity span under the remaining assumptions is zero. Any positive numerical theta_E would be an additional quantitative excitation qualification, not a consequence of the currently numeric envelopes.

Consequently constants.json must remain null/OPEN: filling it from the finite RAO traces would promote statistical evidence into an unsupported source-uniform physical assumption. The quotient point-entry theorem is conditional on an independently justified pair satisfying

    theta_E > Theta_gauge(T_E)
            = min(4 asin(B_a/(2g)), T_E min(D_a/g,B_g,Omega_max)).

Without such a qualification, the source-uniform theorem conclusion must remain stability/retention relative to the measurement-compatible attitude/BA class even during MOVING; point convergence is not identifiable. This is now the controlling assumption gap, not a missing numerical calculation.


## Physical-excitation derivation attempt — no positive source-uniform pair

The requested derivation of a physically justified numerical (T_E,theta_E) from existing evidence was completed and is negative. Repository proof sources intentionally leave the pair symbolic. The pinned 28-ft vessel-RAO data are finite statistical response replays, not an all-time lower-envelope qualification. Linear RAO response scales with incident wave amplitude, so the dataset cannot imply a nonzero response floor for the broader MARINE class without a lower environmental wave-energy/amplitude premise that the theorem does not contain. Encounter-frequency degeneracy likewise prevents manufacturing a universal finite excitation period from vessel speed/heading alone.

Analytically, the existing physical envelopes provide only upper bounds. The exact indistinguishable family can scale alpha -> 0 and nu -> 0 while remaining strictly inside all current B_a,D_a,B_g,D_g,Omega_max, translation, jerk, primitive and magnetic-service bounds. Hence for every proposed T>0 and eps>0 there is an admitted moving family whose complete-window gravity span is positive but below eps (choosing a sufficiently small amplitude and sufficiently slow smooth periodic motion). Therefore no positive theta_E(T) is derivable from the current numerical assumptions: the source-uniform lower envelope is zero.

This means a numerical pair satisfying theta_E>Theta_gauge(T_E) cannot honestly be populated from current repository or generic RAO evidence. Such a pair requires an independently justified quantitative excitation premise (for example a certified operational sea/motion lower envelope), which would strengthen MARINE MOTION and is prohibited merely for proof convenience. Until such evidence is adopted by the theorem contract, retain T_E/theta_E as OPEN and formulate source-uniform stability relative to the measurement-compatible class.


## EXCITED_MOVING minimum physical premise — current

EXCITED_MOVING is now defined separately in ou3-regime-design.md as a proof-side physical subregime, not a runtime mode. It adds no wave-height, spectral, RAO, roll-RMS or estimator-derived condition. Let u_g be the true body gravity direction and define Gamma_ba(h)=2 asin(min(1,min(2 B_a,D_a h)/(2 g_min))). A window is gauge-breaking with margin delta_X>0 if it contains t1<t2 with angle(u_g(t1),u_g(t2)) >= Gamma_ba(t2-t1)+delta_X. This is the direct physical separation needed to defeat the exact attitude/BA packet ambiguity.

For finite entry from a compact outer annulus, one isolated window is not enough. The minimal recurrence premise is: while the same carried execution remains outside the inner target, every T_X interval contained in an EXCITED_MOVING episode contains such a gauge-breaking pair with one fixed positive delta_X. Excitation is not required after inner entry and is not imposed on STILL, TRANSITION or weak MOVING.

This premise is intentionally physical/reference-side; Mahony or OU-III output cannot certify it. Remaining requirements for entry are analytical obligations, not additional physical assumptions: derive compact outer release/retention from existing H18/A21 contracts; prove quotient corrected-word continuity/coercivity and extend the zero-action classification over the outer annulus with literal accepted-update and coupled tuner chronology. If zero quotient action reduces to the measurement-compatible gauge, delta_X excludes it; compactness then yields positive annular dissipation and recurrent windows give finite (possibly history-dependent) inner entry.


## EXCITED_MOVING controlling definition and entry continuation — current

Controlling engineering definition: EXCITED_MOVING is sustained vessel motion in which the vessel undergoes a non-negligible change in roll or pitch within a bounded time. Formally, fixed physical qualification constants T_X<infinity and theta_X>0 require every complete T_X interval contained in an EXCITED_MOVING episode to have true body-gravity-direction span at least theta_X. This is physical/reference-side and is not a runtime detector, Mahony output, OU-III state, minimum wave height, spectrum, heave or RAO-amplitude assumption. Generic weak MOVING remains admitted.

Sufficiency lemma, not definition: the exact packet-indistinguishable attitude/BA family has span at most Theta_gauge(T)=min(4 asin(B_a/(2g_min)), T min(D_a/g_min,B_g,Omega_max)). Therefore theta_X>Theta_gauge(T_X) excludes that exact gauge on every EXCITED_MOVING window. Sensor-bias constants belong here, not in the physical regime definition.

ENTRY CONTINUATION: any fixed strict margin delta_X=theta_X-Theta_gauge(T_X)>0 is topologically sufficient for an existence-level annular dissipation floor IF the following already-open analytical properties are proved on one compact retained outer class: (1) the literal complete-word action is continuous in the carried physical/filter history and root error after quotienting the stationary compatible class; (2) zero action on an EXCITED_MOVING word implies membership in the exact packet-compatible gauge; (3) the compact outer annulus and every-prefix retention are source-uniform. Reason: the closed subset satisfying span>=theta_X is separated from the closed zero-action gauge set by the strict margin; a continuous nonnegative complete-word action therefore attains a strictly positive minimum on each compact annulus. No numerical lower margin beyond positivity is needed for qualitative finite entry. A numerical margin will be needed later for explicit eta_out/finite-error robustness.

This does NOT yet prove entry: properties (1)-(3), especially the outer zero-action classification and compact retained release set, remain open. The next calculation is to prove property (2) for the literal corrected word while carrying Mahony/tuner chronology, not to tune theta_X.
