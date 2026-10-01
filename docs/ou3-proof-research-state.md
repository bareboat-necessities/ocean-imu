# OU-III proof: controlling state after linked finite-supply audit

Audit: PR #639 at `6c96893640abbbc49caeb68a9fe27a27680850b3`, with shipping
sources unchanged on main `17df0e8d64ca07dd6d35ee9a8fea72a85e3eb0d7`.
Current derivation: [linked finite-error supply](ou3-linked-finite-supply.md).
Historical detail remains in Git history and the existing archived ledgers;
retracted calculations are not premises. Runtime, tuning, gates and the
MARINE/IMU/MAGNETIC and residual numerical limits are unchanged.

## Current hypothesis

One path: construction -> capture -> magnetically informed H18 ->
refinement/release -> A21 -> regional practical stability. Retain the conditional
radius-local field-axis result in `ou3-radius-local-field-alignment.md`, with
r_FA=.15, T_FA=17 s and strict local margin .10527117647 m/s^2. Do not assert
source-uniform outer entry merely from this local result.

EXCITED_MOVING retains the engineering definition: sustained vessel motion
with a non-negligible roll/pitch change in a bounded time. Its physical
qualification is true body-gravity-direction span >=theta_X on every complete
T_X window contained in an excited episode, independently of the estimator.
Weak motion is distinct. A bias-only quiet-packet envelope is not a universal
finite-error sufficiency test; all residual supplies must be retained.

## Evidence

For one actual event history write e_N=M e_0+b, with b the signed transported
sum of the exact operation defects. Let J0=P0^-1, JN=PN^-1 and

    G=J0-M'JN M-gamma J0>0, z=M'JN b,
    chi_gamma=b'JN b+z'G^-1 z.

The new exact completed-square identity is

    V_N=(1-gamma)V_0+chi_gamma
          -(e_0-G^-1 z)'G(e_0-G^-1 z).

It couples forcing direction to the same word loss. For a source-uniform
chi upper bound, root retention at C requires chi<=gamma C; finite entry
below r_in requires STRICT reserve chi/gamma<r_in^2. Every-prefix retention
still needs the linked prefix inequalities, not just an endpoint test.

The 80-digit carried 225.00--225.32-s diagnostic preserves observer/control
terminal state and covariance exactly. Wave: rho=.9993070523483296,
gamma=.0003464738258352, chi=306.1085110793. The inner budget ratio is
chi/(gamma*.15^2)=39266523.7608. The fixed-forcing sufficient radius is
939.9451 instead of the separated bound 27673.6030. This is a failed finite
feasibility budget, NOT a source-uniform counterexample: b is retrospective,
and the wave replay does not certify all-time magnetic service.

A separate analytical norm-only finite-residual witness DOES exclude the
proposed full physical-error inner retention for an admitted excited example:
roll=.01 sin(t/2), zero translation, B=75 e_x, constant BA=.01 e_x,
compensating residuals |n_a|<.108067<.3 and |n_g|<=.005<.02. Packets are exactly
quiet while every 60-s window has tilt span .02 rad>1 degree. The actual quiet
MAGNETIC SERVICE certificate survives the true-axis projection. Nominal BA
stays zero and P_ba,ba<=I/1600 implies sqrt(V)>=.4>.15 forever on the regular
real-arithmetic tail. This uses the explicit packet/mean induction, not D=0.
No temporal cancellation of the bounded residuals is assumed in constants.json.
See the full note and `finite_residual_obstruction.py` for all-time bounds.

## Failed approaches / DEAD_ENDS

**Homogeneous/base conflation.** Failed implication: D_W=0 => r_base=0 =>
nominal AW=0 after four S atoms. The action concerns homogeneous measurement
variation, not the base innovation. ZF-5--ZF-8 already explicitly forbid this.
The literal small-x integrated-OU coefficients also are not exact exponential
integrals. Failure: mathematical premise substitution. The purported outer
zero-action exclusion and subsequent homogeneous entry closure are withdrawn.
Retained: correct homogeneous operation loss and the local conditional result.

**Wrong release boundary.** Failed implication: zero LIN at constructor/pre-Live
handoff => zero LIN at A21 release. H18 runs predictions and corrections before
BA activation. Carried example: Live 31.84 s, A21 120.08 s, release LIN norm
.711963. Failure: chronology, not a failed small constant. The source-uniform
H18/A21 release mean box remains open; finite per-history release times do not
supply one common compact release set.

**Unproved compactness.** Failed implication: P>=P_min and V<=C => compact full
history/error set. Scalar P=n^2,e=n gives V=1 with unbounded state/covariance.
Failure: missing coercivity/upper bounds. Retained: regular root covariance
lower comparison and conditional finite-horizon continuity. Every-prefix
outer retention and closed strata must still be proved.

**Linked finite budget.** Failed inequality on the carried wave word:
306.1085110793 <= .0003464738258352*.15^2. This is a feasibility failure of
this fixed-forcing bound. No interval refinement or independent TV/gain/BA
maxima is justified. Retained: exact sharp fixed-word completed square.

**Noise-blind excited entry.** Failed implication: true tilt span above the
bias-only envelope => robust full V entry below .15. The exact finite-residual
witness has sqrt(V)>=.4. Failure: physical indistinguishability under the
norm-only residual model. This is not observer divergence or a refutation of
local homogeneous contraction. A label such as "fast" adds no spectral or
zero-mean premise by itself.

## Retained facts

P_aw,aw<=16.48 I on the audited regular source branch; |e_aw|<4.06 sqrt(V)
only when the ACTUAL error storage is local. BA elimination is
min_b V=e_o'P_oo^-1 e_o, not conditioning. The source S scheduler has bounded
regular gaps; SPD innovation covariance gives accepted S corrections in exact
real arithmetic, not automatically in float32. Exact OU/S identities, signed
variation of constants, covariance energy identities and the coupled
(tau,sigma_aw,R_S,T_S) chronology remain useful. No finite carried replay or
arbitrary covariance box is promoted to a source theorem.

## Current limiter

The requested finite-error .15 target is incompatible with the explicit
norm-only residual witness for the tested EXCITED_MOVING qualification. In
addition the claimed outer homogeneous closure relied on invalid premises.
Source-uniform outer release/retention, correct zero-action/base transfer,
linked supply and every-prefix bounds, regime composition and float32 totality
remain open. All end-to-end theorem flags remain false.

## Alternatives

Use the linked matrix identity on a physically appropriate practical-error
set, retaining the observable consistency class and unavoidable residual tube.
The existing .15 theorem stays conditional where its base-error premise holds.
A stronger temporal/stochastic sensor qualification would be a separate change
requiring justification; none is introduced here. Do not enlarge a covariance
ball and claim compactness, or declare actual innovations zero from D=0.

## Validation and CI boundary

The new exact-rational linked-supply tests and finite-residual witness tests
pass. The native source diagnostic has exact observer/control terminal parity;
its completed-square residual is below 4.4e-78. All-time float32 service and
arithmetic are not certified. Main's existing CI/provenance repairs are
retained; full native `make all` and repository-wide CI are not claimed run by
this mathematical/documentation continuation.

## Next falsifiable calculation

Before another outer entry enclosure, specify a practical physical-error target
consistent with LS11 and prove its retained covariance/history domain. Evaluate
LS1--LS7 on that same-history target with exact OU/BA/S mismatch; a positive
homogeneous loss alone is insufficient. Do not spend enclosure effort trying
to force the already refuted norm-only full-error .15 retention statement.


## Fast-residual temporal qualification — proposed controlling form

The pointwise residual bounds ||n_a||<=.3 m/s2 and ||n_g||<=.02 rad/s are insufficient for physical point-entry. They admit persistent low-frequency residuals that exactly counterfeit genuine roll/pitch while remaining inside the amplitude boxes. A running-mean/DC condition alone is also insufficient: a compensating sinusoid can have arbitrarily small long-window mean while cancelling vessel motion sample by sample.

The weakest natural engineering qualification identified here is therefore a LOW-FREQUENCY RESIDUAL CONTENT envelope, not a smaller instantaneous amplitude. Keep the existing pointwise boxes for fast spikes/vibration, but decompose the already calibrated residual through one declared stable low-pass qualification operator L_X whose passband covers the EXCITED_MOVING attitude band. Require, on every qualified continuation,

    ||L_X n_a|| <= eps_a,LF,
    ||L_X n_g|| <= eps_g,LF,

with the complementary high-frequency residual retaining the existing .3/.02 pointwise/RMS qualification. L_X is a certification/analysis operator, not a shipping filter or estimator change. Its exact transfer function/cutoff and eps bounds must come from stationary/dynamic IMU characterization (PSD/Allan/time-record evidence), not be selected merely to make the proof close.

Why this is minimal: the finite-residual witness uses phi=.01 sin(.5t), i.e. f=.07958 Hz, and requires a compensating gyro residual of amplitude .005 rad/s and accelerometer residual of roughly g*.01=.0981 m/s2 plus the .01 m/s2 DC compensation. Any qualification that still permits those low-frequency components cannot exclude the witness. A pure window-mean bound can permit them. Conversely, bounding the residual after a low-pass that passes the vessel-attitude band directly limits exactly the part capable of masquerading as physical attitude; high-frequency vibration need not be tightened.

Equivalent certification forms are acceptable if proved to imply the same deterministic low-frequency envelope: (a) a PSD/integrated spectral-energy ceiling below a declared f_X, plus a deterministic conversion appropriate to the theorem class; (b) a bank of finite-window sinusoidal/correlation bounds covering [0,f_X]; or (c) a stable low-pass state-space filter with a source-uniform output bound. Allan deviation is useful engineering evidence for selecting/validating timescales and bias/noise decomposition, but by itself is statistical and does not imply the deterministic all-history bound required by the theorem.

NEXT CALCULATION: choose L_X from an independently meaningful vessel/IMU separation timescale and derive the exact modified gauge envelope Theta_gauge,res(T_X) including B_a,D_a,B_g,D_g and eps_a,LF/eps_g,LF. Then determine the maximum allowable eps_a,LF and eps_g,LF for simple candidate EXCITED_MOVING cutoffs such as 1 degree/60 s. These are qualification requirements to compare against real BMI270 data; they are not yet assumptions.


## Low-frequency residual cutoff for the 1-degree / 60-s candidate

This is a NECESSARY gauge-exclusion qualification, not yet the complete finite-error entry budget. For T_X=60 s and theta_X=pi/180, the previously derived bias-only ambiguity envelope is

    Theta_gauge(60)=60 D_a/g = .00611829727786757 rad
                   = .350552611828192 deg.

The remaining physical angular separation is

    delta_X = theta_X-Theta_gauge(60)
            = .0113349952420757 rad
            = .649447388171808 deg.

Let the certification low-pass L_X cover the EXCITED_MOVING band and satisfy deterministic component-independent norm envelopes ||L_X n_a||<=eps_a_LF and ||L_X n_g||<=eps_g_LF. A conservative same-history ambiguity charge is

    Theta_res = 2 asin(eps_a_LF/g_min) + T_X eps_g_LF.

The accelerometer term follows from the maximum gravity-direction chord generated by two low-frequency residual vectors of norm eps_a_LF; the gyro term is the maximum integrated low-frequency rate error over T_X. Therefore a sufficient strict exclusion condition is

    2 asin(eps_a_LF/g_min)+60 eps_g_LF < .0113349952420757.   (LF1)

Axis intercepts are

    eps_a_LF < g sin(delta_X/2) = .0555788680071 m/s2

if the low-frequency gyro residual is negligible, and

    eps_g_LF < delta_X/60 = .000188916587368 rad/s
                            = .0108241231362 deg/s

if the low-frequency accelerometer residual is negligible. A balanced 50/50 angular allocation gives the simple candidate pair

    eps_a_LF < .02778954558 m/s2,
    eps_g_LF < 9.44582937e-5 rad/s = .00541206 deg/s.

These are NOT proposed sensor specifications yet. They are the maximum deterministic low-frequency residual tradeoff implied by the simple 1-degree/60-s EXCITED_MOVING candidate after already charging the existing physical bias-rate envelope. Actual qualification should leave engineering/proof margin below the LF1 boundary and must be demonstrated on calibrated assembled BMI270 devices over temperature and mounting conditions. The Bosch typical broadband noise numbers are not deterministic low-frequency guarantees and cannot by themselves certify LF1.

LF1 excludes the explicit .01-rad/.5-rad-s finite-residual witness: that witness needs about .098 m/s2 low-frequency accelerometer compensation and .005 rad/s low-frequency gyro compensation, far outside the LF1 intercepts. It does not yet prove the linked chi_gamma/gamma<r_in^2 finite-error condition; remaining model/reference/arithmetic supplies still have to be carried in that same-word budget.


## Concrete low-frequency certification operator and BMI270 protocol

Recommended proof/certification operator L_X: an OFFLINE zero-phase low-pass with passband edge 0.20 Hz and stopband beginning 0.30 Hz (nominal separating frequency 0.25 Hz), applied after calibration/temperature compensation to residuals resampled on the qualified 166.7--250 Hz timestamp grid. Zero phase is deliberate: this is a certification operator, not shipping runtime, so phase/group delay should not consume the deterministic residual budget. The passband must have a declared minimum gain (recommend >=.99 on [0,.15] Hz); the stopband attenuation and transition band must be fixed in the qualification artifact. A forward-backward Butterworth/FIR implementation is acceptable only after its actual magnitude response is exported and used in the proof. The theorem should bind the operator coefficients/hash, not the name '0.25-Hz filter'.

Why 0.25 Hz: the explicit residual obstruction is .0796 Hz and therefore remains in-band; a 60-s 1-degree qualification has relevant timescales down to O(.02--.1 Hz). 0.25 Hz leaves substantial margin above these motions while remaining orders below the BMI270 hardware anti-noise bandwidth at the shipping ~200-Hz ODR. This value is an engineering starting point and must be checked against actual vessel roll/pitch spectra; it is not inferred from Bosch typical noise density.

Qualification protocol must estimate TRUE residual, not raw sensor output. Use the exact shipping sensor configuration, calibration, mounting, power and timestamp path. Required assembled-device tests:

1. Stationary six-orientation thermal soak: at least +/-X,+/-Y,+/-Z, covering the deployment temperature range or declared temperature bins. Reference angular rate is zero and reference specific force is gravity in the surveyed orientation. After applying the same calibration/temperature compensation as shipping, form n_a(t),n_g(t), apply L_X offline, and record sup norm and long-window maxima. This certifies low-frequency offset/thermal/creep behavior.
2. Slow single-axis rotation: precision rate table/encoder reference, roll and pitch separately, with sinusoidal/triangular motions spanning .02,.05,.08,.10,.15,.20 Hz and amplitudes including 1--10 degrees. Compute residual after subtracting reference kinematics/gravity, then L_X. This catches scale, cross-axis, phase, mounting and slow dynamic errors that stationary testing cannot.
3. Combined two-axis slow motion: at least representative roll+pitch trajectories in the same band, to prevent a per-axis qualification from missing norm/cross-axis coupling.
4. Temperature repetition during/after slow motion and multiple power cycles. Include assembled PCB/mounting stress, because the theorem is for the device, not a bare BMI270 typical part.
5. Duration: stationary records must be long enough to expose low-frequency drift well below .02 Hz (hours, not minutes); dynamic runs need many cycles per frequency and repeated runs. Exact duration/sample-count is a deployment qualification choice, not proved here.

Pass criterion for the 1-degree/60-s candidate is the JOINT deterministic tradeoff

    2 asin(eps_a_LF/g_min)+60 eps_g_LF < .0113349952420757.

Do not independently require the two axis intercepts. A qualification report should publish the measured worst assembled-device pair (eps_a_LF,eps_g_LF), its margin to this boundary, filter coefficients/response, temperature/mounting envelope, reference-instrument uncertainty, number of devices/runs, and raw-data hashes. Reference uncertainty must be charged into eps_a_LF/eps_g_LF rather than ignored.

For an initial engineering target with 50/50 angular allocation, demonstrate strictly below approximately

    eps_a_LF = .02779 m/s2,
    eps_g_LF = 9.45e-5 rad/s (.00541 deg/s),

and retain additional qualification margin. These are test targets, not yet theorem assumptions. If measured data fail them, first optimize the physically chosen T_X/theta_X pair or calibration characterization; do not silently alter runtime or residual bounds for proof convenience.


## Corrected EXCITED_MOVING homogeneous invariant-set audit

Re-proof starts from the homogeneous variation, not base innovations. On one literal carried base execution, D_W=0 means every HOMOGENEOUS fresh-source factor and corrected measurement variation is zero. It does NOT imply r_acc^base=r_mag^base=r_S^base=0.

The existing corrected-word machinery already gives the exact qualitative nullspace reduction. Zero fresh source action makes one deterministic homogeneous root nuisance trajectory. Four distinct applied S rows form an extended-Chebyshev system for {1,t,t^2,psi_tau}; hence zero homogeneous S loss kills the full homogeneous (v,p,S,a_w) root. Zero homogeneous magnetic loss plus service reduces the root attitude space K_M to dimension <=1. With homogeneous AW zero, one accelerometer row determines the complete BA root from that attitude amplitude; additional rows can only remove the line. Therefore for every complete regular word

    Null(Action_W) subset span(nu_W),

where nu_W=(theta_hat_W,0,...,0,-q_W) is the WORD-DEPENDENT common magnetic/accelerometer compatibility line, possibly trivial. This statement is valid with arbitrary nonzero BASE innovations.

What EXCITED_MOVING must exclude is persistence of nu_W, not the stationary packet gauge directly. Physical normalization shows that if a nontrivial exact compatibility line maps from word W to W+1, its attitude amplitude has |s_W|=1. Active BA OU decay instead forces

    Phi_b,W q_W = s_W q_(W+1),

so indefinite persistence drives |q_W| ->0 while leaving the physical homogeneous tilt amplitude undiminished. At accelerometer epochs this implies

    J_att,k F_k theta_hat_W ->0,

where J_att,k is built from NOMINAL (a_hat_w-g), not true physical force. Thus an indefinitely persistent homogeneous kernel approaches a nominal zero-BA field-axis collinearity condition.

EXCITED_MOVING + the new LF residual qualification constrain TRUE attitude and sensor residual content. They do not by themselves bound the nominal AW/J_att trajectory because base accelerometer/S innovations remain nonzero and can replenish nominal AW. Therefore the implication

    D_W=0 + EXCITED_MOVING + LF residual qualification => e=0

is NOT yet proved. The earlier argument that four S atoms set base AW to zero is withdrawn and must not be revived.

The remaining bridge is now minimal and precise: prove that on one same-history base execution satisfying EXCITED_MOVING and the LF residual envelope, the shipping closed-loop nominal force cannot remain asymptotically field-axis-collinear on every applied accelerometer epoch while a unit physical homogeneous field-axis tilt direction persists. This must use the actual base innovation/AW/S correction chronology and coupled tau,sigma_aw,R_S,T_S law. A sufficient theorem would be a windowed lower bound

    sum_{k in W} w_k ||P_{B_ref,k}(a_hat_w,k-g)||^2 >= c_X >0

on every complete EXCITED_MOVING word (or an equivalent linked action), derived from physical tilt excitation + LF residual bounds + bounded physical primitives and literal innovation dynamics. Pointwise AW tracking is neither needed nor allowed.

NEXT FALSIFIABLE CALCULATION: derive the exact base AW error recurrence at accepted accelerometer and S updates, project it transverse to the committed magnetic field, and combine it over a 60-s EXCITED_MOVING window with the LF residual charge. Test whether the coupled recurrence yields a positive integrated nominal transverse-force floor without replacing innovations by independent controls. If the high-precision same-history feasibility floor is zero/negative, identify the surviving base trajectory before attempting interval enclosure.


## Exact transverse base-AW recurrence and 60-s summation — current

Let b_k=B_ref,k/||B_ref,k|| and P_k=I-b_k b_k'. At the pre-accelerometer epoch define x_k^S=a_hat_w,k^S-g and u_k^S=P_k x_k^S. The literal mean chronology is:

    prediction:  a_hat_w,k^- = phi_k a_hat_w,k-1^+,
    due S row:   a_hat_w,k^S = a_hat_w,k^- - K_aw,S,k S_k^-,
    acc row:     a_hat_w,k^+ = a_hat_w,k^S + K_aw,a,k r_a,k,

with the actual base innovation

    r_a,k = f_meas,k - [Rhat_k x_k^S + lever_k + b_hat_a,temp,k].

Therefore exactly

    u_k^+ = P_k[phi_k a_hat_w,k-1^+ - g
                -K_aw,S,k S_k^- + K_aw,a,k r_a,k].          (TAW-1)

Writing x_(k-1)^+=u_(k-1)^+ + b_(k-1) alpha_(k-1),
alpha_(k-1)=b_(k-1)'x_(k-1)^+, gives the explicit projector-transport form

    u_k^+ = phi_k P_k u_(k-1)^+
            +phi_k alpha_(k-1) P_k b_(k-1)
            -(1-phi_k)P_k g
            -P_k K_aw,S,k S_k^-
            +P_k K_aw,a,k r_a,k.                           (TAW-2)

The second term is the committed-reference rotation charge. TAW-1/2 carry the actual S and accelerometer innovations; homogeneous D=0 does not remove them.

Substitute the physical calibrated measurement model into r_a,k:

    f_meas,k = Rtrue_k(a_phys,k-g_phys)+lever_true,k
               +b_a,phys,temp,k+n_a,k.

Then TAW-1 is a closed same-history recurrence driven by physical a/attitude/bias/residual and the estimator state/covariance/reference chronology. The LF qualification constrains only L_X n_a and L_X n_g; it does not directly constrain r_a because r_a also contains Rtrue a_phys-Rhat a_hat_w, BA mismatch, lever/reference mismatch and base state error.

For a 60-s window W, unrolling TAW-1 gives exactly

    u_N = Phi_(N,0) u_0 + sum_(j=1)^N Phi_(N,j) c_j,          (TAW-3)

where Phi_(N,j) is the ordered product of the literal transverse prediction/update maps and c_j contains the linked gravity-forgetting, reference-transport, S-innovation and accelerometer-innovation terms. The desired action is

    A_X(W)=sum_(k in W) w_k ||u_k^S||^2.                    (TAW-4)

No lower bound on A_X follows from the LF residual envelope alone: LF bounds only one component of c_j. In particular, EXCITED_MOVING constrains TRUE Q(t), whereas u_k is a nominal-force state and actual base innovations provide closed-loop feedback capable in principle of replenishing the gravity-scale AW component. Bounding n_a,n_g in the vessel band removes the explicit quiet-packet residual witness but does not algebraically prevent cancellation through the physical acceleration/BA/base-innovation terms.

Thus the proposed implication

    EXCITED_MOVING + LF residual envelope => A_X(W)>=c_X>0

is NOT established by TAW-1--4 without an additional already-existing physical primitive relation being used. The next valid calculation is to substitute r_a into TAW-3 and eliminate the physical acceleration contribution by the bounded-v/p/S primitives over the SAME 60-s word, while retaining K_aw,a, K_aw,S and the coupled tuner chronology. This is a closed-loop forced-response calculation, not pointwise AW tracking and not an independent-innovation bound. A non-promoting same-history diagnostic should evaluate the resulting signed functional before interval enclosure.


## 60-s physical-acceleration elimination in the transverse AW recurrence — route falsified

Substituting the exact base accelerometer innovation into TAW-3 isolates the physical acceleration contribution as a SAME-HISTORY signed gain-weighted sum

    A_phys(W)=sum_j beta_j a_phys(t_j),

where beta_j is the ordered transported product containing the literal P_B/reference projector, AW accelerometer gain K_aw,a,j, true/nominal frame convention, intervening S/acc/mag mean maps, and coupled tuner/covariance chronology. It is not a scalar averaging weight and is not independent of the physical history.

With h_j=t_j-t_(j-1), w_j=beta_j/h_j and v'=a, exact first Abel summation gives

    sum_j beta_j a(t_j)
      = w_n v_n-w_1 v_0
        -sum_(j=1)^(n-1)(w_(j+1)-w_j)v_j
        -sum_j w_j q_j,

    ||q_j|| <= (J_max/2) h_j^2.

Thus

    ||A_phys|| <= V_max[||w_1||+||w_n||+sum||Delta w_j||]
                  +(J_max/2)sum h_j||beta_j||.              (TAW-A1)

The tempting unweighted 60-s mean bound 2 V_max/60=.18333 m/s2 is therefore NOT applicable through the literal time-varying gains/projectors. It would be valid only for essentially constant scalar weights, which shipping does not supply.

Using p'=v performs a second Abel step and differentiates the gain-weight sequence again. Grouping coefficients inside actual S intervals before taking norms is the strongest already-permitted refinement of this mechanism. The existing same-history calculation gives velocity charge 23.98805965312 and jerk charge 3.82238280532, total 27.81044245844, against a recorded projected-gravity budget 8.77133455729 before root/sensor/BA/spline defects. The actual signed acceleration action on that carried word is about .005, proving that the failure is relaxation/correlation loss, not large physical acceleration.

Therefore bounded v,p,S primitives + EXCITED_MOVING + LF residual qualification do NOT close a positive nominal transverse-force floor through any coefficient-variation/Abel norm bound. Per AGENTS failure protocol this mechanism has had its motivated refinement and must stop. Do not try narrower interval subdivision, second/third Abel variation, or independent gain/projector maxima.

RETAINED EXACT STRUCTURE: TAW-1--3, the signed gain-weighted physical acceleration term, the LF residual qualification, bounded physical primitives, and the tiny carried signed action remain valid. The next route must preserve cancellation between A_phys and the OTHER terms driven by the same accelerometer innovation. In particular K_aw,a multiplies the full innovation

    r_a = Rtrue(a_phys-g_phys)-Rhat(a_hat_w-g_model)
          +BA/lever/residual terms.

Splitting K_aw,a Rtrue a_phys away from -K_aw,a Rhat a_hat_w destroys the closed-loop feedback cancellation. NEXT FALSIFIABLE CALCULATION: combine those two terms before summation and derive the exact affine closed-loop transverse map

    u_k^+ = A_cl,k u_k^S + Kbar_k Rtrue,k(a_phys,k-g_phys)
            + linked BA/lever/LF terms,

with A_cl,k containing I-P_k K_aw,a,k Rhat_k on the transverse subspace. Test passivity/contraction of the COMPLETE pair using the actual Riccati identity K S K'=P^- -P^+ and the S restoring identity, rather than variation of K. This is a new cancellation-preserving mechanism, not another Abel refinement.


## Closed-loop transverse AW passivity calculation — exact result

At an accepted accelerometer correction, with pre-row transverse projector P_B and nominal body rotation Rhat, the AW mean row is

    a_w^+ = a_w^S + K_aw r_a,
    r_a = Rhat(a_phys-a_w^S)+eta,

where eta keeps attitude/true-vs-nominal rotation, gravity mismatch, physical/estimated BA, lever and sensor residual in the chosen consistent convention. Therefore

    u^+ := P_B(a_w^+-g)
         = P_B(I-K_aw Rhat)(a_w^S-g)
           +P_B K_aw Rhat(a_phys-g)
           +P_B K_aw eta_g,                                (CL-1)

with eta_g adjusted so the identity is exact. This is the desired closed-loop map; prediction, reference transport and the due-S correction precede it exactly as in TAW-1/2.

However the Riccati identity does NOT make the transverse AW block passive by itself. Partition the full pre-correction covariance into the selected LIN/AW coordinates X and nuisance coordinates N with cross block C, and H=[H_L,H_n]. Schur elimination gives

    D_c=H_n C' X^-1,
    H_tilde=H_L+D_c,
    R_eff=R+H_n(N-C'X^-1 C)H_n' >0,

and the exact marginal gain K_L=X H_tilde'(H_tilde X H_tilde'+R_eff)^-1. The literal nominal map is A_L=I-K_L H_L=A_tilde+K_L D_c. For any adjoint lambda the exact Joseph balance is

    ||A_L'lambda||_X^2 + ||K_L'lambda||_(R_eff)^2
      = ||lambda||_(X+)^2
        +2 t' D_c X lambda_tilde + ||D_c' t||_X^2,          (CL-2)

where t=K_L'lambda and lambda_tilde=A_tilde'lambda. The RHS extra terms are the exact correlation supply from attitude/BA/other nuisance covariance. They have no fixed sign. Thus K S K'=P^- -P^+ proves passivity of the FULL correction/error storage, but not of transverse AW after deleting nuisance coordinates.

The due S correction is better: its measurement row has no nuisance part, so D_c=0 and

    A_S'(X_S^+)^-1 A_S
      =X_S^-1-H_S'(H_S X_S H_S'+R_S)^-1 H_S,               (CL-3)

an exact restoring loss. For physical error the pseudo target contributes supply nu=-S_phys, so even CL-3 becomes supply-minus-innovation loss rather than pure decay.

Therefore the hoped-for AW-only inequality obtained by combining CL-1 with KSK'=P^--P^+ is FALSE in general. Cross covariance can transfer correction storage between AW and attitude/BA. This is not a numerical looseness and should not be repaired with gain signs or marginal covariance boxes.

RETAINED ROUTE: keep the COMPLETE 21-state innovation storage at accelerometer rows, where the exact information identity is

    Delta V_acc = nu_acc' R_acc^-1 nu_acc
                  -(H e+nu_acc)' S_acc^-1(H e+nu_acc),      (CL-4)

and combine it with the exact S identity

    Delta V_S = nu_S' R_S^-1 nu_S
                -(H_S e+nu_S)' S_S^-1(H_S e+nu_S).          (CL-5)

The same innovations that replenish nominal AW are therefore charged in the full-state loss/supply ledger instead of treated as arbitrary AW forcing. EXCITED_MOVING + LF residual qualification must enter by proving that a persistent word-dependent attitude/BA compatibility direction requires a nonzero sequence of BASE innovation supplies whose linked physical supply is smaller than the corresponding innovation loss. This is the correct passivity target.

NEXT FALSIFIABLE CALCULATION: on the exact persistent compatibility-line ansatz from PT/PER, derive the minimum base accelerometer+S innovation action needed to keep J_att F theta_hat near zero while BA compatibility q_W decays. Compare that required innovation action with the maximum physical supply allowed by bounded v/p/S, LF residuals and BA rates using CL-4/5. This is scalar/line-constrained and preserves full-state Joseph passivity; it avoids AW marginal passivity and gain variation.


## Persistent-line minimum base-innovation action — proposed route fails at zero lower bound

Parameterize a nontrivial homogeneous compatibility line on word W by x=lambda nu_W with nu_W=(theta_hat_W,-q_W), ||theta_hat_W||=1. Exact persistence requires

    J_att,k F_k theta_hat_W = R_ba,k phi_b,k q_W              (PL-1)

at every applied accelerometer row, magnetic compatibility at the applied magnetic rows, and Phi_b,W q_W=s_W q_(W+1), |s_W|=1. These equations are derivatives of the measurement maps with respect to the HOMOGENEOUS perturbation around the carried base execution.

Crucially, PL-1 contains J_att,k (hence the base nominal pre-update a_hat_w) but contains NO base accelerometer innovation r_a,k. The base innovation enters the NOMINAL state recursion that generated J_att, but once that base trajectory is fixed, the homogeneous null condition is independent of the residual magnitude. The same is true for homogeneous S loss: H_S delta x=0 does not constrain r_S^base=-S_hat^base.

Therefore the optimization

    inf { sum r_a,k' S_a,k^-1 r_a,k + sum r_S,k' S_S,k^-1 r_S,k :
          persistent homogeneous compatibility line PL-1 }

has no positive lower bound from the compatibility equations themselves. Algebraically the lower bound is 0: a base trajectory that already lies on the required nominal-collinearity manifold may have arbitrarily small/zero base innovations while the homogeneous measurement derivative retains a nontrivial null line. Full-state Joseph passivity correctly charges innovations that occur, but it cannot create innovation action merely because the linearized observation map is rank deficient.

Hence the hoped-for inequality

    minimum required base innovation action > maximum physical supply

cannot be proved without first proving the missing CLOSED-LOOP REACHABILITY statement that EXCITED_MOVING + LF residual qualification prevents the base nominal trajectory from remaining on/near that collinearity manifold with small innovations. Using Joseph identities before that bridge is circular: they price innovation supply but do not force innovations to exist.

FAILURE CLASSIFICATION: structural/observability, not numerical looseness. Invalidated hypothesis: persistent homogeneous compatibility necessarily requires nonzero base innovation replenishment. Retained facts: full-state CL-4/5 passivity, persistent-line BA decay q_W->0, EXCITED_MOVING, LF residual qualification, and exact base AW recurrence remain valid.

NEXT FALSIFIABLE CALCULATION: attack the collinearity manifold directly with the ZERO/SMALL-BASE-INNOVATION specialization, which is the least favorable case for the proposed supply argument. Set r_a=r_S=0 (or tend them to zero) in the exact base recurrence and ask whether an EXCITED_MOVING physical history satisfying bounded v,p,S, BA rates and LF residual bounds can keep P_B(a_hat_w-g)->0 over successive words. If impossible, compactness/continuity can yield a positive distance/action from the manifold and only then Joseph passivity converts that distance into a finite-error supply bound. If possible, it provides the surviving pathological base execution and point-entry remains false under the strengthened assumptions.


## Magnetic-reference cone audit — canonical geometry proved; quantitative cone still conditional

Literal source sharply narrows the reference problem. Every continuous hard-iron reference write goes through ContinuousHardIronTracker::slewTowardEstimate. It does NOT adopt an arbitrary fitted 3-D field direction. The tracker preserves the canonical world gauge and writes

    B_ref=(h,0,z),
    h=h_anchor+(h_new-h_anchor),
    z=z_anchor+(z_new-z_anchor),

with h>MAG_INIT_MIN_MAG_NORM. The continuous estimator uses the independent yaw-stripped startup-proxy tilt, and its fit is gated by information, residual RMS and max hard-iron fraction. Therefore continuous adaptation cannot freely rotate the nominal reference azimuth; approach to gravity can occur only through collapse of h relative to |z|.

However the literal source floor is only MAG_INIT_MIN_MAG_NORM=1e-3 uT. By itself this yields no useful source-uniform lower cone. The physical MAGNETIC SERVICE constants B_h,true>=15 uT and 20<=|B_true|<=75 uT apply to the physical field/accepted information; they cannot simply be substituted for nominal h_ref.

A useful cone follows IF an already-qualified all-time reference-vector error bound ||B_ref-B_true||<=E_B is available. Then

    h_ref >= 15-E_B,
    |B_ref| <= 75+E_B,
    ||P_Bref g|| = g h_ref/|B_ref|
      >= g (15-E_B)/(75+E_B),                               (MR-1)

for E_B<15. If one were entitled to combine the existing 5-uT hard-iron residual and 2-uT measurement residual into E_B=7 uT, MR-1 would give

    ||P_Bref g|| >= 9.80665*(8/82) = .956746341463... m/s2.

But that implication is NOT presently proved. The 5-uT and 2-uT constants bound physical residual components used by MAGNETIC SERVICE/sampling certificates; they are not a theorem that the continuously learned world reference remains within 7 uT of the true field. The continuous tracker itself permits a fitted bias up to .35 times measured field norm and gates fit residual RMS at 3 uT; neither alone gives ||B_ref-B_true||<=7 because proxy-tilt error and identifiability enter the level-frame fit.

The new LF accelerometer/gyro residual qualification helps bound the startup-proxy tilt error in the vessel-motion band, but a quantitative all-time proxy-tilt/reference error transfer has not yet been derived. Thus the desired source-uniform positive nominal cone is CONDITIONAL, not closed.

NEXT FALSIFIABLE CALCULATION: derive E_B directly from the canonical reference update formula ref(b)=mean(R_i m_i)-mean(R_i)b. Write R_i=R_true,i Delta_i and m_i=R_true,i' B_true+b_true+n_m. Bound B_ref-B_true in terms of (i) proxy tilt error Delta_i, (ii) fitted/applied hard-iron error, and (iii) magnetic measurement residual. Use the estimator's 3-uT fit-residual and .35-field bias gates only where they mathematically bound those terms. Then combine the new LF IMU qualification with the proxy observer dynamics to see whether E_B<15 uT is actually certifiable. If not, the reference cone remains an explicit independent qualification requirement.


## Magnetic-reference error decomposition — exact; proxy-tilt bound is the limiter

For the continuous hard-iron accumulation let Q_i be the TRUE body->level tilt rotation and let the proxy leveling rotation be R_i=Delta_i Q_i, where Delta_i is the tilt-frame error. Use the physical magnetometer model

    m_i=Q_i' B_true + b_true + n_m,i

in the yaw-stripped level convention (any fixed yaw gauge is absorbed into B_true). For an applied hard-iron correction b_app, the estimator's exact level reference is

    B_ref = mean(R_i m_i)-mean(R_i)b_app.

Substitution gives the exact decomposition

    B_ref-B_true
      = mean[(Delta_i-I)B_true]
        +mean[Delta_i Q_i(b_true-b_app)]
        +mean[Delta_i Q_i n_m,i].                           (MR-2)

Because rotations preserve norm, if angle(Delta_i)<=delta_proxy on the accumulation window, ||b_true-b_app||<=E_HI and ||n_m||<=E_m, then

    ||B_ref-B_true||
      <= 2 B_max sin(delta_proxy/2)+E_HI+E_m
      =: E_B.                                               (MR-3)

No independence/zero-mean assumption is used. The same bound applies to a weighted mean with nonnegative normalized weights. Continuous partial application is covered by E_HI for the actually applied b_app. The canonical re-gauging in slewTowardEstimate preserves the horizontal/vertical level-field changes rather than arbitrary azimuth rotation; MR-3 is therefore the natural physical error budget to qualify that update.

Combining with the physical field envelope gives, whenever E_B<15 uT,

    ||P_Bref g|| >= g (15-E_B)/(75+E_B) >0.                 (MR-4)

Numerically, if HI and magnetic residual charges were temporarily zero, E_B<15 requires only

    delta_proxy < 2 asin(15/(2*75)) = .2003348423 rad
                = 11.47834 deg.

With a hypothetical combined E_HI+E_m=7 uT, the proxy requirement for merely positive cone is

    delta_proxy < 2 asin((15-7)/150)
                = .106717... rad = 6.115... deg,

and the resulting E_B=7 special case (if independently justified) would give ||P_Bref g||>=.956746 m/s2.

CURRENT LIMITER: the new LF IMU residual qualification bounds only sensor residual content. It does NOT by itself give angle(Delta_i)<=delta_proxy for the private Mahony proxy, because the proxy accelerometer sees genuine marine translational/wave acceleration as well as gravity. The proxy gains are deliberately below the wave band, but a deterministic source-uniform tilt-error transfer from bounded physical v/p/a/jerk + LF residuals to the Mahony observer has not been proved. Therefore MR-3/4 are conditional and E_B<15 is not yet certified.

Also, the existing constants hard_iron_residual_norm_max=5 uT and measurement_residual_norm_max=2 uT cannot simply be inserted as E_HI and E_m without checking their exact physical definitions against MR-2. In particular the continuous estimator's 3-uT fit RMS and .35 field-scale bias gate are estimator diagnostics, not direct pointwise bounds on b_true-b_app.

NEXT FALSIFIABLE CALCULATION: derive a deterministic low-frequency tilt-error bound for VerticalAccelComplementary/Mahony in the yaw-free two-axis subsystem. Treat physical horizontal acceleration as the disturbance, use bounded velocity/displacement to bound its low-frequency component through the proxy correction transfer, and add the certified LF accelerometer/gyro residual envelopes. The target need only be modest (roughly <6 deg if a 7-uT non-tilt reference budget is later justified), not sub-degree. If the proxy transfer has nonzero DC gain from arbitrary bounded marine acceleration that cannot be controlled by the existing primitives, the reference cone must remain an independent qualification rather than be derived from EXCITED_MOVING.


## Private Mahony proxy tilt transfer — linearized two-axis result

Shipping VerticalAccelComplementary uses the private IMU-only Mahony observer with twoKp=.2 and twoKi=.02. In Mahony_AHRS the cross-product error is a half-vector, so the small-angle yaw-free tilt subsystem has effective k_P=twoKp/2=.1 s^-1 and k_I=twoKi/2=.01 s^-2 (sign chosen for stable negative feedback). For one horizontal axis, with tilt error e, integral-feedback/bias state beta, horizontal specific-force disturbance u=a_h/g, gyro residual d_g and accelerometer-direction residual d_a,

    e_dot = -k_P e + beta + d_g + k_P(u+d_a),
    beta_dot = -k_I e.                                      (MP-1)

Thus the translational-acceleration-to-tilt transfer is

    H_a(s)=k_P s/(s^2+k_P s+k_I),                           (MP-2)

while gyro-rate residual enters through

    H_g(s)=s/(s^2+k_P s+k_I).                               (MP-3)

The PI loop rejects constant gyro bias but H_a has the expected low-frequency behavior: H_a(0)=0 with positive k_I; near sqrt(k_I)=.1 rad/s it can pass substantial slow acceleration. This is why the integral term can wind up against sustained horizontal acceleration, as the source comment states.

For sinusoidal physical horizontal acceleration at angular frequency omega, bounded physical velocity gives A_a<=omega V_max. Hence the corresponding linearized proxy tilt amplitude obeys

    |e_a| <= (V_max/g) omega |H_a(j omega)|.                 (MP-4)

Using V_max=5.5 m/s and the literal gains, the values for .02,.05,.08,.10,.15,.20 Hz are approximately 3.67,3.37,3.28,3.25,3.23,3.22 degrees. This is below the ~6.1-degree proxy target associated with a hypothetical 7-uT non-tilt magnetic-reference budget. It is a feasibility result only: a sinusoidal spectral component is not a deterministic arbitrary-history theorem.

The LF accelerometer residual enters MP-1 as direction error approximately eps_a_LF/g and the LF gyro residual through H_g. Their deterministic contribution can be added only after defining the same certification band/operator for the proxy theorem. The 50/50 test targets are small enough to evaluate, but no theorem is promoted here.

CURRENT LIMITER: existing MARINE bounds |v|<=V_max and |p|<=P_max do not directly imply a source-uniform L_infinity bound on the output of the stable convolution H_a applied to arbitrary a=v'. A crude L1 impulse-response times |a|<=A_max is too loose and ignores the derivative structure. The correct deterministic route is to integrate H_a(s)*s V(s) by parts: define G_v(s)=s H_a(s)=k_P s^2/(s^2+k_P s+k_I), separate its direct k_P term from the stable strictly-proper remainder, and exploit BOTH |v| and the finite-window/position primitive to control the low-frequency part. Alternatively formulate the proxy error directly as a stable state driven by bounded v through an integration-by-parts storage. Do not replace the arbitrary history by sinusoidal decomposition without a spectral norm theorem.

NEXT FALSIFIABLE CALCULATION: derive an induced bound from bounded v and p for MP-1 over the 30-s magnetic refinement / 600-s continuous-HI windows. Compute the exact impulse kernels for the literal k_P,k_I and evaluate the sharp endpoint + L1-kernel constants after integration by parts. If the resulting deterministic proxy bound plus LF sensor charges is < the MR-3 cone budget, close E_B; if it exceeds it badly, proxy tilt/reference cone must be independently qualified from device/vessel data rather than inferred from the broad MARINE envelopes.
