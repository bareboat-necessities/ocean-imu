# Covariance-weighted AW adjoint cancellation

This is a subordinate lemma for `ou3-joint-aw-passivity.md`, on the same
construction/capture/H18/release/A21 proof path. It introduces no estimator
change or separate stability architecture. It is an algebraic result, not a
numerical diagnostic or a proof of the 1.96133 m/s^2 mean bound.

Scope: the qualified regular real-arithmetic profile, optimal accepted
corrections using their actual effective measurement covariance, and the
literal operation chronology. Finite-precision/source defects, projections
and reconfiguration retain their existing separate qualifications.

## 1. S corrections have exactly zero weighted-adjoint jump

Use the twelve-component mean z=(a_w,v,p,S), with each block a world vector,
and X=P_LL. At an S correction let H=E_S and

    Omega = H X H^T + R_S,
    K = X H^T Omega^-1,
    A = I-K H.

Here K is the LIN part of the actual full gain; because the S row observes
only LIN, it equals this expression even with arbitrary inherited cross
covariances. The marginal covariance update gives

    X^+ = X-K Omega K^T = A X = X A^T.

Between readout insertions the exact nominal-mean adjoint obeys
lambda^-=A^T lambda^+. Define chi=X lambda at the corresponding boundary.
Then

    chi^- = X A^T lambda^+ = X^+ lambda^+ = chi^+.

In particular,

    sup_(reachable applied S events) ||chi^+-chi^-|| = 0.    (SC)

No sign restriction on K_aw,S or P_aw,S is needed. This holds for every
positive definite X, so it holds on its shipping-reachable subset. It is
NOT a claim that beta itself has zero jump or small raw variation.
A readout atom c is a separate operation, adding X c to chi. An ordinary
attitude error reset leaves the LIN marginal and reduced mean unchanged;
its effect on other covariance blocks and later gains is still retained.

## 2. The exact nonzero jump at a nuisance-coupled observation

Partition the full prior covariance and its actually used Jacobian as

    P^- = [[X,C],[C^T,N]],           H=[H_L,H_n].

Let K_L be the actual LIN gain. Then

    K_L Omega = X H_L^T + C H_n^T,
    X^+ = X-K_L Omega K_L^T.

The nominal LIN mean, with the realized nuisance history kept in the forcing,
uses A_L=I-K_L H_L, not the Schur-adjusted observation multiplier. Consequently

    chi^- - chi^+
       = (K_L Omega-X H_L^T) K_L^T lambda^+
       = C H_n^T K_L^T lambda^+.                         (NC)

Proof: substitute lambda^-=A_L^T lambda^+ in X lambda^- and subtract
X^+ lambda^+. No bound or sign assumption is used.

For S, H_n=0 and (NC) recovers (SC). For active lever-arm-disabled
accelerometer updates,

    C H_n^T = P_L,theta J_theta^T + P_L,ba.

Held BA must use its literal covariance-update row. For a magnetic update,
H_L=0 and (NC) equals K_L,m Omega_m K_L,m^T lambda^+. The actual LIN mean
also receives K_L,m r_m; its transported source cannot be omitted.

For comparison, on the FULL homogeneous error adjoint p^-=(I-KH)^T p^+,

    P^- p^- = P^+ p^+,
    K^T p^+ = R^-1 H P^+ p^+.

These identities hold across every optimal measurement, but the full error
adjoint is not the reduced nominal lift with nuisance coordinates deleted.
Equation (NC) is precisely the term that deletion would hide.

## 3. Prediction, synchronization and scale commits

If the actual LIN covariance prediction is

    X_new = F X_old F^T + Q + Delta,
    lambda_old = F^T lambda_new,

with invertible F and the actual pending PSD AW increment Delta, then

    chi_old = F^-1 [chi_new-(Q+Delta)lambda_new].           (PC)

An isolated covariance increment leaves the nominal mean and lambda unchanged
and changes chi by Delta lambda. It is not a direct AW mean input.
Readout atoms, (NC), (PC) and literal reset/scale maps are the complete
weighted-adjoint chronology. No S-gain jump is charged independently of the
covariance change that exactly cancels it in (SC).

In the isotropic settled scalar normalization, put

    D=diag(sigma I3, sigma*tau I3, sigma*tau^2 I3, sigma*tau^3 I3),
    x=h/tau,      zeta=2 sigma^2 tau/r_a,      c_T=T_S/tau.

The ideal OU process and transition in these coordinates depend only on x.
For the SpectralMSE target law the normalized S covariance is

    Rtilde_S = r_S^2/(sigma^2 tau^6)
       = C_J^2 4^(1/7)/(c_T c_sigma^(12/7)) * zeta^(-1/7).

The realized applied values require the additional factor ell^2, where

    ell = r_S sqrt(T_S) /
          [C_J (2 r_a)^(1/14) (sigma/c_sigma)^(6/7) tau^(24/7)].

This is an exact definition of the target-to-applied discrepancy, constrained
by the carried smoothers, commits, clamps and floors; ell is not set to one.
For general realized scalar R_acc,

    Rtilde_acc = [h R_acc/r_a] * 2/(zeta x).

At a scale commit the normalized mean transforms by D_new^-1 D_old. This
known coordinate map must accompany the normalized covariance and adjoint.

## Consequence

The direct S-event contribution to chi variation is uniformly ZERO. The
physical beta variation still depends on (NC), (PC), actual measurement rows,
readout atoms, and applied schedule/scale changes. A small complete-word bound
on these terms has not been proved. In particular, a small observed signed
response does not bound raw first/second variation, and failure of a sufficient
variation certificate is not a construction-reachable instability example.


## 4. Complete covariance-weighted reader theorem

The remaining nominal-mean problem can be stated without raw beta variation.
Freeze the literal coefficients on one same-history word W, but do not restart
the root, tuner, scheduler or physical source.  Factor every covariance/action
channel exactly as in the existing finite-error composition:

    P_i = A_i P_(i-1) A_i^T + B_i B_i^T.

At prediction B_i is a factor of the literal Q_i (including its within-step
correlations); at an accepted correction its fresh factor is K_i R_i^(1/2);
at a congruent reset it is zero.  PSD synchronization increments are separate
prediction-like factors.  Let s be the single stacked vector of all these
whitened source coordinates, in chronological order.

For the scalar transverse nominal-mean readout q_W x (q_W includes the
convex AW-window weights and the fixed world direction u), variation of
constants gives exactly

    q_W x_W = q_W M_W e_0 + Z_W s + d_impl,                 (CR1)

where

    Z_W = [ q_W M_(W<-i) B_i ]_i.                           (CR2)

No independence or Gaussian premise is used: the factorization is matrix
algebra.  The same physical/source block is represented once, so all of its
multiple chronological appearances are summed before its norm is taken.

Let

    C_root^2 = q_W M_W P_0 M_W^T q_W^T,                    (CR3)
    C_src^2  = Z_W Z_W^T.                                  (CR4)

Then covariance Cauchy--Schwarz gives the exact reader inequality

    |q_W x_W|
      <= C_root sqrt(V_0) + C_src sqrt(A_W) + |d_impl|,     (CR5)

where A_W=s^T s is the complete stacked source action.  This is stronger than
splitting accelerometer and magnetic terms by the triangle inequality.

If bookkeeping requires source families, write the columns of Z_W as the
disjoint concatenation

    Z_W=[Z_phys,Z_acc,Z_mag,Z_proc,Z_sync,...].

For A_f=||s_f||^2,

    |Z_W s|
      <= sqrt( sum_f C_f^2 ) sqrt( sum_f A_f ),             (CR6)
    C_f^2=Z_f Z_f^T.

Alternatively the weaker requested display follows from familywise
Cauchy--Schwarz,

    |q_W x_W|
      <= C_root sqrt(V_0)
         + C_phys sqrt(A_phys)
         + C_acc sqrt(A_acc)
         + C_mag sqrt(A_mag)
         + C_proc sqrt(A_proc)
         + C_sync sqrt(A_sync)
         + C_impl.                                          (CR7)

The physical bounded-primitive contribution can instead be reduced by the
signed Abel/cycle identity before it is inserted into Z_phys; doing so yields
the requested notation C_phys(Vmax,Pmax,Jmax) without treating physical
acceleration as a stochastic OU source.

### Terminal-covariance normalization

The existing covariance recursion supplies

    P_W = M_W P_0 M_W^T + Z_state Z_state^T.

Whitening by P_W^(-1/2) gives a horizontal matrix with operator norm <=1.
For any terminal scalar row q,

    C_root^2 + C_src^2
      <= q P_W q^T                                         (CR8)

when q is applied to the corresponding terminal/readout augmented state.
For a multi-time AW mean, augment the state with the deterministic readout
accumulator; the same identity holds because readout insertion has no fresh
source.  Thus the root and source coefficients are not independent constants:
they share one covariance budget.

This is the key improvement over

    C_root sqrt(V0)+C_acc sqrt(A_acc)+C_mag sqrt(A_mag).

For any nonnegative V0,A_W,

    C_root sqrt(V0)+C_src sqrt(A_W)
      <= sqrt(C_root^2+C_src^2) sqrt(V0+A_W)
      <= sqrt(q P_W q^T) sqrt(V0+A_W).                      (CR9)

No square root of the number of corrections appears.

### Relation to the minimum-action historical reader

For the augmented chronological design y=O_h h0+A s and readout/terminal
quantity r=T_h h0+T s, every linear reader L satisfying L O_h=T_h leaves

    r-Ly = (T-LA)s.

For a scalar readout q, the minimum possible source coefficient is therefore

    C_min^2
      = q [ Pi + Ttilde I_eff^(-1) Ttilde^T ] q^T,          (CR10)

with

    Pi=T(I-A^T Sigma^(-1)A)T^T,
    Ttilde=T_h-T A^T Sigma^(-1)O_h,
    I_eff=O_h^T Sigma^(-1)O_h,
    Sigma=A A^T.

This is exactly the existing joint minimum-action reader B*.  Hence the
nominal-mean source theorem does not require a new observability architecture:
it is a scalar dual use of the already proved reader algebra.  The useful
bound is the signed factor-space evaluation of (T-L*A), not a product of
per-event gain norms.

## 5. How the requested source families enter

For an accepted correction with deterministic physical/model residual r_i,
Joseph form gives P_i >= K_i R_i K_i^T.  In the whitened correction channel
the source coordinate is

    s_i = R_i^(-1/2) r_i,

so its action is exactly

    A_i = r_i^T R_i^(-1) r_i.                              (CR11)

Therefore

    A_acc = sum_(i in acc) r_acc,i^T R_acc,i^(-1) r_acc,i,
    A_mag = sum_(i in mag) r_mag,i^T R_mag,i^(-1) r_mag,i. (CR12)

They are square-summed before multiplication by the complete reader.
The correlation term C H_n^T K_L^T lambda derived above is already contained
in the same full-state correction factor and must not be added again.

Prediction/model mismatch is handled identically with the least-norm
coordinate in a factor of Q_i.  PSD AW covariance synchronization has no mean
residual; it changes the reader/covariance chronology and contributes a factor
to the auxiliary covariance budget, but it is not a physical mean source.
Tuner lag likewise changes the literal coefficients.  Its convex-recursion
bounds qualify the coefficient set; it is not added as an independent AW
forcing unless the implementation differs from the committed coefficient.

Direct magnetic mean forcing K_m r_m is exactly the magnetic member of
(CR11--CR12), not an extra Euclidean gain term.

## 6. Physical acceleration is better treated deterministically

The physical acceleration entering the accelerometer is one bounded-primitive
history, not an independent white source.  Keep its complete signed coefficient
beta_j from the chronological reader and apply the exact identities already
proved:

    sum beta_j a_j
       = endpoint velocity terms
         + signed first-difference velocity term
         + jerk quadrature,

or the cycle/jerk decomposition for rapid scheduler modulation.

Define C_phys(W) to be the resulting deterministic bound using the SAME beta
rows generated by the complete reader.  Then

    |u^T mu_W|
      <= C_root sqrt(V_0)
         + C_phys(W;Vmax,Pmax,Jmax)
         + C_meas sqrt(A_acc+A_mag+A_other)
         + C_impl,                                          (CR13)

where C_meas may be taken as the Euclidean norm of the concatenated
measurement-source reader rows.  Familywise C_acc,C_mag is valid but weaker.

The carried root is not discarded.  On a retained region V_0<=r_0^2 its
contribution is at most C_root r_0.  If the word has homogeneous contraction
margin delta_W, the terminal-state part additionally satisfies the existing

    ||M_W e_0||_(P_W^-1) <= sqrt(1-delta_W) sqrt(V_0).

For the multi-time nominal mean use the augmented readout accumulator rather
than substituting this terminal inequality blindly.

## 7. Symbolic field threshold

Let b be the applicable unit reference field and

    sigma_w = ||e_z x b||.

Corollary A* requires exactly

    ||mu_W x b|| < g sigma_w.                              (CR14)

A sufficient scalar bound is ||mu_W||<g sigma_w.  Keep this symbolic in the
reader theorem.  If the theorem domain separately assumes sigma_w>=1/5, the
right side is at least

    g/5 = 1.96133 m/s^2.

If the only field-domain information is |I|<=80 degrees under
sigma_w=cos I, the uniform right side is instead

    g cos(80 deg) ~= 1.7029069 m/s^2.

Reference-field and magnetic-model defects belong once in the magnetic
residual/action or in the geometric field premise, according to the chosen
model; do not subtract the same defect twice.

## 8. Quantitative closure condition

The exact sufficient condition on every qualified word is

    C_root(W) r_0
      + C_phys(W)
      + C_src(W) sqrt(A_W)
      + C_impl(W)
        < g sigma_w,min,                                   (CR15)

or the sharper joint two-vector form

    sqrt(C_root(W)^2+C_src(W)^2)
      sqrt(r_0^2+A_W)
      + C_phys(W)+C_impl(W)
        < g sigma_w,min.                                   (CR16)

Here C_src is the minimum-action/signed-factor coefficient, not a maximum
single-event gain.  This is the requested complete covariance-weighted reader
theorem.

What is still OPEN is numerical/analytic evaluation of the source-uniform
suprema of C_root,C_src,C_phys,A_W,C_impl on the retained coupled
Riccati/tuner/physical trace class.  The reader formula itself, the
square-summed measurement action, the absorption of the nuisance-correlation
term, and the symbolic threshold are closed.  Existing O1/O2 reader machinery
can be reused to bound C_src because (CR10) is its scalar dual; no new
observability lemma is required.


## 9. Noncircular kernel decomposition and the actual remaining obstruction

The complete reader removes the artificial G0 circularity but does not, by
itself, make every source coefficient small.  Decompose the relative attitude
into the magnetically observed plane and the rotation about the applicable
unit field b.  The observed component is charged in the magnetic
measurement-action family.  For the field-axis component M_b b=b, the exact
world accelerometer identity contains

    M_b a_phys + g(e_z-M_b e_z) + beta_a,

where beta_a is the carried world accelerometer bias/sensor residual after
using one consistent true/nominal rotation convention.

The gravity part has the SHARP kernel bound

    ||g(e_z-M_b e_z)|| <= 2 g sigma_w sin(theta/2).          (K1)

This is smaller than the generic 2g sin(theta/2) by sigma_w.  At the retained
theta=6 deg and sigma_w=1/5 it is 0.2052961622 m/s^2.

The declared post-projection accelerometer-bias error and fast residual give

    B_ba + N_a = 0.625166604983954 + 0.3
               = 0.925166604983954 m/s^2.                   (K2)

Thus the non-physical kernel box excluding translated acceleration is

    K0 = 1.1304627671 m/s^2                                 (K3)

at sigma_w=1/5.  This is below g/5, so the static magnetic kernel alone does
not create the forbidden nominal mean.

The physical translated term is not bounded by A_max.  Since a=dv/dt,

    (1/T) int M_b a dt
      = [M_b v]_0^T/T - (1/T) int dot(M_b) v dt,             (K4)

and therefore

    ||mean(M_b a)||
      <= 2 Vmax/T + (Vmax/T) TV(M_b).                       (K5)

The trapezoidal sampled version adds the already proved J h_max/4 sampling
term plus the corresponding discrete variation of M_b.

This identifies the exact noncircular obstruction: the present retained
domain bounds sup angle(M_b)<=6 deg but does NOT bound TV(M_b) by a constant
independent of T.  Lemma I* bounds the NET ordered injection rotation, not
its total variation.  The estimator axial gyro-bias projection gives an
amplitude sector, not a signed temporal-variation bound.  Deterministic fast
sensor residuals likewise have an amplitude bound but no variation premise.

Consequently (K5) cannot be made small by MARINE MOTION/IMU BIAS/MAGNETIC
SERVICE boxes alone.  Replacing TV(M_b) by T times the gyro-sector amplitude
is finite but quantitatively useless.  This is not evidence that such a
field-axis oscillation is shipping-reachable: M_b is an estimator error, not
a freely selectable physical input.

### Why a new observability lemma would be circular

At the forbidden boundary,

    f_hat=a_hat-g e_z parallel b

is exactly the condition

    ||a_hat x b|| = g sigma_w

at its minimum-norm representative.  The accelerometer attitude row then
annihilates the field-axis direction.  Therefore any argument that first
assumes an accelerometer/magnetic six-column floor in order to control M_b
and then uses that control to prove ||mu_hat x b||<g sigma_w is circular.
The needed result must instead come from shipping reachability of the
field-axis error/AW/bias loop.

### Coupled reachability lemma that would close the gap

A sufficient noncircular lemma is a signed temporal bound

    TV_W(M_b)
      <= C_bg,0 + C_bg,V sqrt(V0)
          + C_bg,a sqrt(A_acc) + C_bg,S sqrt(A_S)
          + C_bg,impl,                                      (K6)

where every coefficient is obtained from the literal axial gyro-bias,
accelerometer-bias, AW/S and reset recursion WITHOUT using G0 or the nominal
AW mean premise.  Substitution in (K5), followed by (K1--K3), feeds directly
into (CR15/CR16).

The BA-rate law can help only jointly.  In a zero-residual compatibility
calculation, changing field-axis attitude changes the gravity compensation at
rate of order g sigma_w |dot theta_b|; D_a=0.001 would by itself force a very
slow theta_b.  But AW correction and the deterministic fast accelerometer
residual can share that compensation, so D_a alone is not a proof of (K6).
The literal coupled AW/S and axial-bg recursion must be retained.

### Minimum source action formulation

Equivalently, define the constrained complete-word action

    A_min(m) = inf A_W

over all shipping-reachable chronological traces with

    ||mu_hat x b|| >= m,

the carried root and tuner state fixed only by the retained set, and all
physical histories satisfying the declared contracts.  The desired theorem
is

    A_min(g sigma_w) > A_available.                         (K7)

The full-reader identities show that S corrections and optimal measurement
correlations are already represented without gain-sign assumptions.  OU
leakage alone cannot establish (K7): at tau=12 s the per-step correction
needed to replenish a g/5 DC AW component is only

    (1-exp(-h/12)) g/5
      in [0.00065367,0.00098042] m/s^2

for h in [0.004,0.006], far below the declared accelerometer residual
envelopes.  A successful lower action must therefore come quantitatively from
the integrated S-chain together with the axial-bg/BA chronology.

At the weak regularizer corner tau=12, sigma=4, the SpectralMSE target exceeds
the shipping r_S clamp, so the applied target is capped at r_S=100 m*s and
T_S at 0.15 s (before carried smoothing/commit lag).  Hence this corner must
be included explicitly in any uniform K7 proof; assuming a strong S update
there would be invalid.

## 10. Result

The G0 circularity is removed from the statement of the nominal-AW theorem,
and the static field-axis kernel is bounded sharply by (K1--K3).  The
complete-reader/source-action formulation is exact.  But the desired strict
bound is NOT derivable from the currently proved independent envelopes:
the unclosed quantity is now the signed temporal reachability/total variation
of the shipping field-axis attitude-error loop, equivalently the constrained
minimum source action (K7).

This is a reachability obligation, not another geometric observability lemma
and not a request for a stronger MARINE MOTION assumption.  Proving (K6) or
(K7) from the literal axial-bg + BA + AW/S recursion is the next analytical
step.  A carried maximum or an independently chosen TV(M_b) box would be
fitted and must not be substituted.


## 11. Coupled axial loop: necessary pathology rate and no-go for total action

This section attacks (K6/K7) directly.  It obtains a quantitative necessary
condition for a forbidden nominal-AW word, but also proves that the proposed
TOTAL source-action comparison is not the correct closing functional under
the present deterministic sensor contract.

### 11.1 Free axial rotation cannot reach the field threshold

On the field-axis branch let M_b b=b and let theta_b be its signed angle.
Separate the part of dot(theta_b) generated without estimator corrections:
physical gyro-bias drift/residual plus the commissioned fast gyro residual.
The declared amplitude/rate channel gives the conservative free angular-rate
bound

    Omega_free = N_g + D_g = 0.02001 rad/s,

where charging D_g as a rate is conservative on a one-second normalization;
using the exact bias-history integral can only improve the word bound.

For a complete T-second window, integration by parts gives

    ||mean(M_b a_phys)||
       <= 2 Vmax/T + Vmax Omega_free
          + correction-induced term.                        (AX1)

The sampling-fidelity defect adds J h_max/4.

The static field-axis nuisance bound is

    K_static =
       2 g sigma_w sin(theta_max/2)
       + B_ba,post + N_a.                                   (AX2)

At the worst declared field fraction sigma_w=1/5 and retained
theta_max=6 deg,

    2 g sigma_w sin(3 deg) = 0.205296162115946,
    K_static = 1.130462767099900 m/s^2.

Choose a 32-s nominal-mean window; this changes no physical assumption and is
inside the recurring A21 word.  Then

    2 Vmax/T + J h_max/4 = 0.49375 m/s^2,
    Vmax Omega_free       = 0.110055 m/s^2.

Therefore every forbidden word must obtain at least

    Delta_corr =
      g/5 - K_static - 0.49375 - 0.110055
      = 0.227062232900100 m/s^2                             (AX3)

from correction-induced field-axis rotation/rectification.  Equivalently,
under the integration-by-parts relaxation, it needs average additional
variation at least

    Omega_corr,needed = Delta_corr/Vmax
                      = 0.0412840423454727 rad/s.            (AX4)

Thus physical gyro/bias/fast-noise transport alone cannot sustain the
gravity-scale nominal AW pathology.  Any such shipping trajectory must use
the estimator correction loop itself at a quantitatively nontrivial rate.

For a general field fraction retain the symbolic margin

    Delta_corr(sigma_w,T)
      = g sigma_w
        - 2 g sigma_w sin(theta_max/2)
        - B_ba,post - N_a
        - 2 Vmax/T - J h_max/4
        - Vmax Omega_free.                                  (AX5)

Only when this is positive does AX4 give a useful necessary correction rate.

### 11.2 Why total measurement action cannot close K7

The deterministic sensor contract bounds each fast residual in amplitude; it
does not impose stochastic cancellation or a finite all-time l2 budget.
Consequently

    A_acc(T)=sum r_acc,k^T R_acc,k^-1 r_acc,k

and the analogous magnetic action may grow linearly with the number of
samples even for an admitted coherent bounded residual.  Extending T therefore
does not make A_available small.  The minimum action required to replenish OU
leakage also grows linearly with T.  A comparison

    A_min(g sigma_w) > A_available

based only on TOTAL square-summed action is therefore structurally incapable
of exploiting the fact that the target is a DC/signed-mean quantity.  This
invalidates K7 as the final scalar closure, while retaining the complete-reader
factorization for finite-error supplies.

### 11.3 The required functional is low-frequency correction transport

Let delta theta_k^c be the field-axis part of the actual attitude correction
and let M_b,k be the resulting carried field-axis error.  The exact dangerous
term is not sum |delta theta_k^c| and not sum NIS_k.  It is the signed pairing

    R_corr(W) =
       (1/T) sum_k v_k^T (M_b,k^+ - M_b,k^-) + reset defects, (AX6)

or its exact SO(3) counterpart before linearization.

A sufficient shipping lemma is

    |R_corr(W)| <= C_corr < Delta_corr(sigma_w,T).           (AX7)

At sigma_w=1/5,T=32 s it is enough to prove

    C_corr < 0.227062232900100 m/s^2.                        (AX8)

Equivalently, the conservative TV version needs only
0.0412840423454727 rad/s average correction-induced variation, but AX6 is
strictly preferable because prediction/correction sawteeth that keep M_b near
zero cancel before the norm is taken.

### 11.4 Exact covariance identity for correction-induced axial motion

For an accepted correction with full gain K, innovation covariance Omega and
field-axis attitude selector q_b, put

    delta theta_b = q_b^T K r,
    Delta P_b = q_b^T K Omega K^T q_b >=0.

Then Cauchy--Schwarz in measurement space gives exactly

    |delta theta_b|^2
       <= (r^T Omega^-1 r) Delta P_b.                       (AX9)

Meanwhile Joseph form gives

    q_b^T P^+ q_b = q_b^T P^- q_b - Delta P_b

before the reset congruence.  Thus large axial corrections consume axial
covariance information.  Prediction replenishes that covariance only through
the literal gyro/gyro-bias process block and reset transport.

Equation AX9 is noncircular: it uses no G0 and no nominal-AW premise.
However, summing |delta theta_b| with Cauchy--Schwarz introduces a
sqrt(number-of-corrections) loss and is quantitatively useless.  The next
valid operation is to insert AX9 into the SIGNED pairing AX6 and telescope
Delta P_b against prediction replenishment before taking norms.

### 11.5 Exact remaining certificate

Define the chronological axial correction factor

    z_k = sqrt(Delta P_b,k) sign-compatible with q_b^T K_k,

and whiten the corresponding residual coordinate so
|u_k|^2<=NIS_k and delta theta_b,k=z_k u_k in the scalar relaxed channel.
Transport every occurrence to the physical-velocity pairing before squaring:

    Xi_b = [ signed coefficient of each fresh gyro/process/measurement factor
             in R_corr(W) ].

Then

    |R_corr(W)|^2 <= (Xi_b Xi_b^T) A_corr(W).                (AX10)

The point is that Xi_b is formed AFTER the Joseph decrements and prediction
replenishments telescope.  No per-event gain norm, NIS sum, raw TV, or G0
geometry enters.  This is the axial scalar analogue of the signed
factor-space reader already used for R_q/R_d.

The source-uniform theorem now required is

    sup_(reachable W)
      (Xi_b Xi_b^T) A_corr(W)
        < Delta_corr(sigma_w,T)^2.                           (AX11)

For the 32-s, sigma_w=1/5 domain the right side is

    Delta_corr^2 = 0.0515572575 (m/s^2)^2.

AX11 is narrower than the previous nominal-AW reader problem: it concerns one
field-axis signed correction/velocity pairing.  It uses only the literal
gyro/gyro-bias process covariance, correction Joseph decrements, reset
transport, and deterministic physical velocity bound.  S/AW/BA enter through
the actual correction factors but do not require a separate observability
floor.

What is proved here is AX1--AX10 and the numerical necessary margin AX3/AX4.
AX11 is NOT yet evaluated source-uniformly.  Therefore the shipping
impossibility theorem is not claimed closed.  A proof that merely replaces
AX11 by total NIS/action or raw correction TV would repeat a demonstrated
quantitative relaxation failure.


## 12. Persistent exact-degeneracy lemma: closed without G0

The exact branch condition itself supplies a simpler contradiction than the
minimum-action construction.  Let b be the applicable unit committed field,
sigma_w=||e_z x b||>0, and suppose on a complete interval W of duration T
the literal nominal specific force satisfies

    f_hat(t_k)=a_hat_w(t_k)-g e_z parallel b                 (DG1)

at every applied accelerometer epoch used by the normalized signed mean.
Then for every such epoch

    P_perp a_hat_w = P_perp(g e_z) =: c,                    (DG2)

where P_perp=I-bb^T and

    ||c||=g sigma_w.                                        (DG3)

Thus exact degeneracy forces a FIXED transverse nominal AW component; lambda(t)
along b is irrelevant.

Use the exact world accelerometer identity, projected by P_perp.  With a
consistent true/nominal attitude convention it can be written

    P_perp a_hat_w
      = P_perp M_b a_phys
        + P_perp g(e_z-M_b e_z)
        + P_perp e_ba
        + P_perp n_a
        - P_perp r_acc,                                     (DG4)

where M_b is the remaining field-axis attitude mismatch after the magnetically
observed component is charged in the magnetic residual/action, e_ba is the
physical-minus-estimated residual accelerometer bias, n_a the commissioned
fast accelerometer residual, and r_acc the actual innovation.  On an exact
zero-innovation compatibility branch r_acc=0.  More generally its normalized
signed mean must be retained as an explicit residual term.

For the exact compatibility branch average DG4 with the nonnegative
trapezoidal weights of the window.  The physical acceleration term obeys the
sampling-fidelity identity

    ||mean a_phys||
       <= 2 Vmax/T + J h_max/4.                             (DG5)

A rotation about b leaves b fixed, so its gravity defect has the sharp bound

    ||g(e_z-M_b e_z)||
       <= 2 g sigma_w sin(theta_max/2).                     (DG6)

Use the declared universal post-projection estimator/physical BA error
B_ba,post and fast accelerometer residual N_a.  Since projection cannot
increase the norm,

    g sigma_w
      <= 2 Vmax/T + J h_max/4
         + 2 g sigma_w sin(theta_max/2)
         + B_ba,post + N_a.                                 (DG7)

Therefore exact persistent degeneracy is impossible whenever

    T >
    2 Vmax /
    [ g sigma_w(1-2 sin(theta_max/2))
      - J h_max/4 - B_ba,post - N_a ],                      (DG8)

provided the denominator is positive.

For the declared worst field fraction sigma_w=1/5,

    g sigma_w                         = 1.96133,
    2 g sigma_w sin(3 deg)           = 0.205296162115946,
    J h_max/4                         = 0.15,
    B_ba,post                         = 0.625166604983954,
    N_a                               = 0.3.

The remaining physical-mean allowance is

    D = 1.96133 - 0.205296162115946
        -0.15 -0.625166604983954 -0.3
      = 0.680867232900100 m/s^2.

Hence

    T_crit = 11/D = 16.156012... s.                         (DG9)

So every exact zero-innovation field-axis-degenerate interval longer than
16.157 s is excluded by the declared same-history physical/bias/sensor
envelopes.  A 17-s complete window has strict margin

    g/5 - [11/17 + 0.15 + 0.205296162115946
           +0.625166604983954+0.3]
      = 0.033808527... m/s^2.                               (DG10)

This proof uses NO G0, no AW covariance ceiling, no S-gain sign, no carried
0.348/0.371 value, and no strengthened MARINE MOTION assumption.  OU and S
can only affect how the estimator attempts to remain on DG1; they cannot
alter the algebraic requirement DG2.

### Nonzero-innovation robust version

If the branch is only measurement-compatible up to an accelerometer residual,
define its signed transverse mean

    Rbar_acc =
      || sum_k alpha_k P_perp Rhat_k^T r_acc,k ||.           (DG11)

Then DG7 becomes

    g sigma_w
      <= 2 Vmax/T + J h_max/4
         +2 g sigma_w sin(theta_max/2)
         +B_ba,post+N_a+Rbar_acc.                            (DG12)

Consequently a 17-s branch is excluded whenever

    Rbar_acc < 0.033808527... m/s^2                         (DG13)

at sigma_w=1/5.  This is a SIGNED residual-mean requirement, not a pointwise
innovation or total-NIS requirement.  It is exactly the quantity the complete
covariance-weighted reader should bound; coherent high-frequency residuals
cancel if their signed mean does.

The previous claim that T_crit was about 21.44 s was arithmetic overcharging:
it included the free gyro/velocity rotation term even though DG4 already
works in the field-axis-rotated physical acceleration frame.  For exact
degeneracy the direct projected measurement identity gives the sharper
16.156-s threshold above.

### Consequence for the single proof path

The feared persistent exact trajectory is now eliminated on any 17-s
qualified complete interval, before invoking G0.  Therefore G0 may be applied
noncircularly AFTER establishing that every 17-s interval either

1. leaves the field-axis-degenerate set by a definite geometric amount, or
2. carries a signed accelerometer residual mean at least the DG13 threshold.

Case 2 is not yet impossible under the deterministic sensor contract: the
declared fast residual is amplitude-bounded but need not have zero signed
mean, and DG13 is much smaller than its 0.3 m/s^2 envelope.  Thus the exact
zero-innovation pathology is closed, while a near-degenerate persistent
trajectory can still hide in a coherent signed innovation unless the
complete reader/S-chain bounds Rbar_acc.

This is the precise remaining robustification.  It is much narrower than the
old AW-mean reachability problem: prove a source-uniform signed transverse
accelerometer-innovation mean below the geometric escape margin, OR use its
nonzero value directly as information/action that forces departure from the
degenerate set.  Do not bound it by the pointwise 0.3 residual envelope.
