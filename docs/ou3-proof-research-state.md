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
