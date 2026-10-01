# OU-III shipping reference and prefix-premise audit

Source snapshot: PR #637 at `3476d2a7e8269c2811b01056baa9f6fb4fb67b66`.
This note attacks the two premises left by FA1--FA12. It does not change
runtime, physical assumptions, tuner laws, quality gates, or proof clocks.

## 1. Literal continuous-reference transport

After refinement, the default continuous hard-iron path is exogenous but not
constant. Let the estimator's exponentially weighted statistics at an accepted
solve be

    Abar = sum(lambda_i R_i)/sum(lambda_i),
    wbar = sum(lambda_i R_i m_i)/sum(lambda_i).

For an applied body offset b the estimator returns exactly

    L(b)=wbar-Abar b.                                      (MR1)

The wrapper anchors the currently stored canonical reference B_* and the
currently applied offset b_* when the first valid continuous estimate arrives.
For a later applied offset b, it writes

    B_ref(b) =
      ( B_*,h + ||L(b)_h||-||L(b_*)_h||, 0,
        B_*,z + L(b)_z-L(b_*)_z ).                         (MR2)

This is the literal source formula. Since every R_i is orthogonal and the
weights are nonnegative,

    ||Abar||_2 <= 1.                                       (MR3)

Therefore reverse triangle inequality gives the exact source-level Lipschitz
bound

    |h_ref(b)-h_ref(b_*)| <= ||b-b_*||,
    |z_ref(b)-z_ref(b_*)| <= ||b-b_*||,
    ||B_ref(b)-B_ref(b_*)|| <= sqrt(2)||b-b_*||.           (MR4)

The slew supplies a time-local version. If b_tgt is the current accepted
target and alpha=1-exp(-dt/45),

    ||b^+-b|| = alpha ||b_tgt-b||.                         (MR5)

MR1--MR5 are shipping bounds. They use the actual reference chronology and do
not introduce a compatibility control.

## 2. Why the present gates do not imply the needed all-time cone

The startup/refinement MagAutoTuner canonicalizes an accepted reference and
requires horizontal_fraction >= .05. Hence at those writes,

    h_ref >= .05 ||B_ref||.                                (MR6)

The continuous path does not reapply MR6. Its final wrapper check is only

    h_ref > MAG_INIT_MIN_MAG_NORM = .001 uT,               (MR7)

plus finiteness. The continuous estimator accepts a fitted body offset up to

    ||b_fit|| <= .35 field_scale,                           (MR8)

and the deployed application fraction is one. Its residual-RMS, information
and ridge gates can reject a bad fit, but none is an algebraic lower bound on
h_ref/||B_ref||. In particular MR4+MR8 alone cannot preserve even the declared
physical 15-uT horizontal lower bound: the theorem permits measured field
scale up to 75+5+2=82 uT, and .35*82=28.7 uT > 15 uT.

This is NOT a constructed physical counterexample. It proves that the literal
continuous-reference safety gates, by themselves, do not imply the nominal
cone needed by FA9.

More importantly, the continuous statistics use startupProxyTiltQuat(), the
private measurement-only Mahony tilt, for every R_i. The current theorem has
no all-time deterministic bound

    angle(R_proxy(t), R_true(t)) <= eps_proxy              (MR9)

after startup/refinement. The gravity gate is used for handoff; it is not an
all-time invariant gate on the continuous hard-iron accumulation. Therefore
the physical assumptions ||B||<=75, B_h>=15 and the 5+2-uT residual envelopes
cannot yet be transported through MR1 to an all-time bound on B_ref. Treating
R_i as true attitude would assume precisely the missing result.

A sufficient analytical bridge is explicit. If MR9 is proved with
2 sin(eps_proxy/2)<=delta_R, and calibrated magnetic data obey

    m_i = R_true,i^T B_0 + d_i,  ||d_i||<=eps_m,

then

    ||wbar-B_0|| <= 75 delta_R + eps_m,                    (MR10)

before the fitted-offset term. MR1 then gives, for any applied b,

    ||L(b)-B_0|| <= 75 delta_R + eps_m + ||b||.            (MR11)

The sharper version subtracts the actual physical hard-iron component before
bounding the fitted residual; that requires a declared decomposition/accuracy
bound, not merely the current aggregate 5-uT envelope. MR10--MR11 identify the
missing quantity exactly. They do not assign it a value.

## 3. Prefix retention: the FA12 radius cannot be inherited from six-degree capture

FA7 used a full-state covariance storage premise

    V_i=e_i'P_i^-1 e_i <= r^2

at EVERY pre-accelerometer prefix. In the fixed-reference 80-degree
specialization, positivity requires

    r < (g cos80 - .26)/4.06
      = 0.3553957885... .                                  (PR1)

The illustrative r=.25 gives the published .4279-m/s^2 margin.

At proxy handoff the shipping attitude covariance is seeded with tilt standard
deviation .035 rad and zeroed stale attitude cross-covariances. Even ignoring
all other state error, a tilt error theta has minimum storage

    sqrt(V) >= |theta|/.035.                                (PR2)

Thus the boundary of the existing six-degree capture region gives

    sqrt(V) >= (pi/30)/.035 = 2.99199... ,                 (PR3)

over eight times the largest radius allowed by PR1. This is not a statement
that handoff actually has six degrees of error; it proves that the existing
six-degree capture theorem target does NOT imply the FA12 storage premise.

The startup gravity gate also cannot repair PR3 by reinterpretation. Its
threshold .075 is a residual of a low-passed world-frame accelerometer in the
proxy frame, not a deterministic true-attitude error bound. The timeout path
requires only the aligned branch, not the trusted gate. Converting .075 to a
true tilt angle would discard physical acceleration/bias terms and change the
meaning of the implemented gate.

## 4. Why retention is logically downstream of the open word return

For an actual finite-error word the existing storage calculus has the form

    sqrt(V_{j+1}) <= q_j sqrt(V_j) + E_j,                  (PR4)

with corresponding prefix bound

    sup_{k in word j} sqrt(V_k) <= G_j(sqrt(V_j)).          (PR5)

The actual-gain operation inequalities needed to construct E_j and G_j exist,
but a source-uniform strict q_j<1 is exactly what O1/O2 and the linked block
return are still trying to prove. A bounded covariance marginal such as
P_aw,aw<=16.48 I does NOT bound deterministic estimation error and cannot be
used to replace PR4.

Therefore unconditional prefix invariance V<=r^2 cannot be proved first and
then used to prove the very word contraction needed for PR4. That would be
circular. The correct closure is simultaneous on a candidate rectangle, as
already required by the controlling obligations:

    D(c,r)<=c,
    q(c,r) r + E(c,r) <= r,
    G(c,r)<=r,                                             (PR6)

with the field-axis geometric lower bound inserted into the SAME word/block
loss that determines q(c,r).

## 5. Noncircular reformulation of the field-axis exclusion

FA7 actually needs only a deterministic AW-error tube, not full-state storage.
If at every relevant prefix

    ||a_hat_w-a_phys|| <= eps_aw,                           (PR7)

then the fixed-reference argument is immediately

    sum w_i c_i^2 >= [G-C_phys-eps_aw]_+^2.                (PR8)

For the 100-s, 6-ms, 80-degree specialization,

    eps_aw < 1.4429069015... m/s^2                         (PR9)

is sufficient. This is much weaker than V<=.25^2 through the generic
P_aw,aw<=16.48I conversion. For a varying reference, the same substitution
gives

    M_ref =
      G0-delta_g-eta_ref(Amax+g_model_max+eps_aw)
      -eps_aw-C_phys.                                      (PR10)

PR7 must still be proved from the shipping error dynamics; covariance alone
does not prove it. But PR8--PR10 remove the artificial requirement that the
ENTIRE 21-state Mahalanobis error already be tiny before field-axis geometry
can help create contraction.

## 6. Result and next decisive calculation

PROVED from shipping source:
- MR1--MR5, the exact continuous-reference Lipschitz/slew bounds;
- MR6--MR8, the exact distinction between startup and continuous gates;
- PR1--PR3, the incompatibility between six-degree capture as presently stated
  and the small FA12 full-storage radius;
- PR8--PR10, the component-tube reformulation needed for noncircular closure.

NOT PROVED by the current theorem contract:
- MR9, an all-time true-to-private-proxy tilt tube, or an equivalent direct
  all-time nominal-reference cone;
- PR7, an all-prefix deterministic AW tracking tube;
- the simultaneous retained rectangle PR6.

Accordingly the requested two unconditional premises cannot honestly be
promoted from the present assumptions. This is not because arbitrary
compatibility controls were considered: every equation above follows the
literal shipping reference update or the actual covariance/error storage.

The next calculation should NOT try another independent reference or covariance
box. It should compose the private Mahony error equation and the AW
accelerometer-correction error equation on the same physical history, deriving
MR9 and PR7 together with the already coupled tau/sigma_aw/R_S/T_S schedule.
Those component tubes then enter BP-4/BP-10 directly; solve PR6 simultaneously.
If the Mahony dynamics cannot supply MR9 under the existing marine/bias
contract, that is a genuine missing theorem premise or runtime certificate,
not a covariance algebra problem.
