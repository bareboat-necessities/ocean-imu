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


## 13. Correction: innovation is endogenous; robust DG12 is not a source bound

The zero-innovation special case DG1--DG10 is a valid incompatibility
calculation.  The proposed nonzero-innovation extension DG11--DG13 is NOT a
closure theorem: r_acc is the filter innovation

    r_acc = f_meas - f_pred,

hence it is an endogenous function of the same nominal AW/attitude/BA state.
Moving its signed mean to the physical-source side and then attempting to
bound it independently is tautological.  A persistent signed innovation is
precisely how the filter corrects a bad nominal state.

Retain DG1--DG10 only as the exact zero-innovation subcase.  Replace the
robustification by the following posterior-loop identity.

Let c=P_perp(g e_z), |c|=g sigma_w.  If the POSTERIOR nominal force remains
exactly degenerate after each accepted correction, then

    P_perp a_hat_k^+ = c.                                  (PL1)

Let a_k^- be the predicted AW immediately before the accelerometer correction,
and collect every non-accelerometer AW mean correction in xi_k.  The actual
accelerometer AW increment is

    Delta a_k^acc = K_aw,k r_k.

Projecting PL1 and summing gives the exact required correction balance

    sum_k P_perp K_aw,k r_k
      = sum_k [ c-P_perp a_k^- ]                            (PL2)

over the accelerometer epochs, with the prediction/S/magnetic chronology
inside a_k^- rather than replaced by a free input.  Using the literal AW-loop
telescoping identity,

    sum_acc Gamma_k(e_k-eta_k)
      = e_0-e_N + sum xi_k
        -sum_pred[(1-phi_k)a_hat_k+Delta a_phys,k],          (PL3)

shows that physical acceleration increments telescope.  On PL1 the OU part
contains the fixed term

    sum_pred (1-phi_k)c,                                    (PL4)

while the S corrections contribute their actual restoring increments through
xi_k.  Thus maintaining the branch requires the accelerometer corrections to
replenish at least the OU+S loss of c, modulo endpoints and magnetic/reset
terms.

Crucially PL2/PL3 must now be combined with the SAME innovations
r_k=f_meas-f_pred.  There is no independent Rbar_acc source.

### Required gain-weighted incompatibility lemma

A sufficient robust bridge is:

For every qualified carried chronology and every interval W on which

    ||P_perp f_hat_k|| <= delta_f                            (PL5)

at all accepted accelerometer epochs,

    || sum_k P_perp K_aw,k r_k
       - required_replenishment_W(c) ||
       <= E_phys+E_BA+E_sensor+E_impl,                       (PL6)

and the right side is strictly smaller than the OU+S replenishment required
when |c|-delta_f >= g sigma_w-delta_f.

Equivalently, after substituting r_k, prove a positive lower bound on

    sum_k <c_hat, (I-Gamma_k)(c-P_perp y_k)>
      + S-dissipation                                      (PL7)

where Gamma_k=K_aw,k Rhat_k and y_k is the physical world acceleration plus
the exact attitude/BA/sensor term.  This is a closed-loop passivity inequality,
not an innovation bound.

The S update already has exact covariance-metric dissipation and zero
deterministic source.  The remaining quantitative issue is the effective
accelerometer-to-AW operator Gamma_k on the field-transverse direction under
the reachable Riccati/tuner chronology.  An arbitrary PSD covariance cannot
be used; doing so reintroduces the original algebraic cancellation loophole.

Therefore the controlling lemma is now precisely:

    inf_(reachable k, field-transverse u)
       u^T Sym(Gamma_eff,k) u >= gamma_acc > 0              (PL8)

in the COMPLETE OU+S cycle metric, or, more generally, the word version

    sum_k <u, Gamma_k u> + D_S(W,u)
       >= gamma_W sum_k |u|^2,                              (PL9)

with gamma_W>0 source-uniformly.  D_S is the exact nonnegative S-chain
dissipation after eliminating its internal state.  PL9 is allowed to hold
only over a complete scheduler cycle/word; pointwise gain signs are not
required.

If PL9 is proved, bounded physical velocity makes the signed physical
acceleration mean O(1/T), BA has the declared amplitude/rate bounds, and the
fast deterministic sensor residual contributes through the same positive
closed-loop operator rather than as an independent 0.3 box.  The forbidden
posterior c then cannot be an invariant mean.

This is the original Riccati-reachability/passivity question in its correct
minimal form.  The exact-degeneracy calculation usefully identifies c and
the needed word length, but it does not remove PL9.


## 14. Quantitative cycle target after covariance synchronization

The shipping wrapper performs posterior AW covariance maintenance every
ADAPT_EVERY_SECS=0.1 s in Live, independently of tuner enablement.  The default
PSD-floor path stages

    Delta = Pi_+(Sigma_aw-P_aw)

inside the next prediction.  With S_factor=1 and the applied sigma floor,

    Sigma_aw >= sigma_min^2 I,   sigma_min=0.05 m/s^2.       (CY1)

Thus every 0.1-s cycle contains a prediction whose AW marginal is restored in
the deficient eigendirections before subsequent measurements.  This prevents
permanent collapse of AW measurement leverage.

In a scalar decoupled comparison, an accelerometer row H_aw=1 with
P_aw>=sigma_min^2 and R_acc<=r_acc,max^2 would give

    Gamma >= sigma_min^2/(sigma_min^2+r_acc,max^2).          (CY2)

With r_acc,max=0.3010398645 this is about 0.02684.  CY2 is NOT valid for the
full shipping gain because attitude/BA/LIN cross covariance can alter the AW
row of P H^T.  It is only the scale indicating that a 17-s word contains
roughly 170 covariance-maintenance opportunities, so a modest complete-cycle
floor would suffice.

The correct source-uniform cycle statement is

    sum_(k in cycle) <u, Gamma_k u>
       + D_S(cycle,u)
       >= gamma_cycle |u|^2,                                (CY3)

for every field-transverse unit u and every covariance reachable immediately
after the scheduled AW floor, where D_S is the exact nonnegative homogeneous
S-chain dissipation expressed in the same Schur/completed-square metric.
Prediction/process covariance and the PSD floor are included in the cycle
start; tuner tau,sigma,r_S,T_S are the coupled applied tuple.

A proof of CY3 may use the Schur complement of nuisance coordinates rather
than P_aw alone.  Let X=P_aw, C=P_aw,n and N=P_nn.  The conditional AW
covariance

    X_c = X-C N^dagger C^T >=0                              (CY4)

is the part of AW uncertainty not explainable by nuisance coordinates.
For the accelerometer observation, after conditioning on nuisance, the
effective scalar/vector gain is generated by X_c with effective noise
R_eff>=R_acc.  Therefore a positive cycle floor follows if the scheduled
PSD AW floor supplies a uniform positive lower bound

    X_c >= p_c I                                            (CY5)

at at least one accepted accelerometer update per 0.1-s cycle.

CY5 is the exact covariance-reachability sublemma.  The marginal floor
P_aw>=sigma_min^2 I does NOT imply CY5: pre-existing AW/nuisance correlation
can in principle make the Schur complement zero.  The scheduled PSD increment
helps because it is added ONLY to the AW block while preserving cross
covariances.  If immediately before the floor the joint covariance is
[[X,C],[C',N]], then after adding Delta>=0 to AW,

    X_c^new = X_c^old + Delta.                              (CY6)

Hence every strictly positive eigenvalue of Delta is added one-for-one to the
conditional AW covariance.  The remaining case is an eigendirection u in
which Delta u=0, meaning u'X u is already at/above the stationary target in
the spectral-positive-part sense.  Such a direction needs a separate
conditional-correlation argument; marginal largeness alone is insufficient.

This yields the precise dichotomy for CY5:

1. floor-active direction: u'Delta u>=d0 -> u'X_c^new u>=d0;
2. floor-inactive direction: P_aw is already large in u, and either its
   conditional variance is positive, or AW is almost deterministically
   encoded by nuisance coordinates.  In the latter case the S-chain and
   nuisance measurements observe that encoding, so the complete-cycle
   dissipation D_S must be used instead of the direct accelerometer gain.

Thus the desired cycle floor is naturally a sum
direct-conditional-accelerometer information + S/nuisance dissipation, not a
pointwise K_aw sign statement.

What remains unproved is a quantitative lower bound gamma_cycle in CY3 over
the compact reachable covariance/tuner set.  The existing coarse nuisance
upper box proves compactness after 17 s but is too large for a useful
Schur-complement number.  A sharp proof must exploit CY6 and the exact S-chain
cancellation before taking block norms.


### 14.1 Floor-active direct gain bound

Condition on every nuisance coordinate at the post-floor prior.  In a
field-transverse scalar direction u let

    p = u^T X_c u.

The conditioned accelerometer observation has unit AW coefficient (rotation
preserves norm) and effective noise variance r_eff>=R_acc; nuisance
conditioning removes, rather than boxes, the correlated cancellation.  Its
posterior conditional AW update has scalar information increment 1/r_eff and
mean gain

    gamma_c = p/(p+r_eff).                                  (CY7)

If u^T Delta u>=d0, CY6 gives p>=d0.  With the shipping accelerometer
effective variance bounded above by

    r_eff <= r_acc,max^2 + r_unconditioned,

a useful numerical lower gain requires an UPPER bound on the remaining
conditioned measurement noise, not merely R_eff>=R_acc.  If all nuisance
coordinates are conditioned exactly, the only physical measurement-noise
term is R_acc, so in the formal conditional filter

    gamma_c >= d0/(d0+r_acc,max^2).                          (CY8)

This is an information decomposition, not the literal marginal K_aw entry.
The complete-reader proof may allocate the fresh Delta component to this
conditional channel and the nuisance reconstruction to the complementary
reader.

### 14.2 Floor-inactive directions carry stored conditional information

If Delta has zero quadratic form in u, the sync supplies no fresh conditional
variance.  There are two possibilities:

(a) u^T X_c u >= p0.  Then the same conditional accelerometer information
argument applies with p0.

(b) u^T X_c u < p0.  Then AW in direction u is already determined to accuracy
p0 by the nuisance state.  In covariance language there exists the linear
conditional predictor

    a_u = C_u N^dagger n + epsilon_u,
    Var(epsilon_u)<p0.                                      (CY9)

The dangerous DC AW component is therefore carried almost entirely by the
nuisance coordinates.  But the S-chain identity annihilates AW/root/sync
columns and leaves a signed observation of the neutral integrator/noise
coordinates.  A persistent DC value encoded in nuisance must consequently
appear either in the S=0 innovation or in the bounded v/p/S root/endpoints.
Thus case (b) converts the AW problem into the already exact S-chain
readout, with residual proportional to sqrt(p0).

This proves a qualitative cycle dichotomy for every p0>0:

    direct conditional AW information >= p0/(p0+r_acc,max^2)
    OR
    AW is p0-close (in conditional variance) to an S-chain-observed nuisance
    predictor.                                                (CY10)

No arbitrary covariance cancellation survives this dichotomy.

What is still needed for a numerical gamma_cycle is the quantitative
coefficient mapping the nuisance predictor C_u N^dagger n in CY9 into the
S-chain readout.  The existing S-chain formula gives that coefficient exactly
for a supplied chronology; a source-uniform lower singular value on the
relevant predictor subspace has not yet been derived.  This is now the only
covariance-geometric quantity in the robust cycle lemma.


## 15. Four-row S information: explicit coupled lower value

For four equally spaced S rows at h=T_S and fixed applied tau over the four-row
microblock, eliminating the quadratic nuisance columns (v,p,S) is exact by a
third finite difference.  With

    psi(t)=tau^3[(t/tau)^2/2-t/tau-expm1(-t/tau)],

the vector orthogonal to the three polynomial columns at h,2h,3h,4h is
(-1,3,-3,1), whose squared norm is 20.  Hence the scalar Schur complement for
a unit initial AW component is exactly

    kappa_S(tau,h,rS)
      = [Delta_h^3 psi(h)]^2/(20 rS^2)
      = tau^6 exp(-2h/tau)(1-exp(-h/tau))^6/(20 rS^2).       (KS1)

This is an analytical formula, not a fitted/replay quantity.

The shipping variables are coupled:
h=clamp((0.015/1.1)tau,0.005,0.15), and the default SpectralMSE rS target is
clamped to [0.15,100] after its coupled tau/sigma law.  The minimum of KS1 on
the deployed coupled envelope occurs at the short-tau floor corner

    tau=0.02 s, h=0.005 s, rS=0.15 m*s.

There

    Delta_h^3 psi
      = tau^3 exp(-h/tau)(1-exp(-h/tau))^3
      = 6.7432167881e-8 m*s per (m/s^2),

and therefore

    kappa_S,min = 1.0104660589e-14                         (KS2)

in the corresponding weighted S-action units per squared AW amplitude.

The independent-box corner tau=.02,h=.15,rS=100 is unreachable and must not
be used.  Conversely KS2 is positive but quantitatively tiny: even hundreds
of disjoint four-row blocks do not by themselves give a useful gravity-scale
mean contradiction.

This does NOT invalidate the complete-cycle lemma.  At the corner producing
KS2, OU decay is strongest: over one 0.1-s covariance-sync interval the
homogeneous AW multiplier is exp(-0.1/.02)=exp(-5)=0.00673795.  Thus the
S-chain is the wrong mechanism to price the short-tau corner.  The uniform
cycle certificate must retain OU leakage + direct conditional accelerometer
information + S information jointly.  At long tau, where OU leakage is weak,
KS1 is much larger before the rS clamp and the pseudo cadence is slower; at
the rS=100 corner its absolute S-action is still small, so the accelerometer
channel remains essential.

Conclusion: the useful source-uniform constant is NOT kappa_S,min alone.
KS2 closes the requested four-row Schur-complement calculation and proves
that an argument of the form '170 cycles times kappa_S,min' cannot close the
near-degenerate branch.  The next scalar minimization must be

    gamma_cycle =
      inf_reachable [ D_OU + D_acc,conditional + D_S ],      (KS3)

with all three terms evaluated in one normalized covariance/information
metric.  Taking their separate global minima would again combine incompatible
corners and lose the coupled chronology.


## 16. Joint-cycle minimization: S-only floor is not the closing constant

The explicit four-row result is

    kappa_S(tau,h,rS)
      = tau^6 exp(-2h/tau)(1-exp(-h/tau))^6/(20 rS^2).

On the coupled shipping schedule its minimum is
1.0104660589e-14 at tau=.02,h=.005,rS=.15.  This is too small to close a
17-s word by multiplication.

The correct joint quantity cannot be the arithmetic sum of three unrelated
numbers called D_OU,D_acc,D_S.  OU is deterministic mean leakage, the
accelerometer is endogenous feedback driven by the physical measurement, and
S is a zero-valued pseudo-measurement.  They must be combined in the same
input/output map.

For one conditioned transverse AW direction, eliminate l=(v,p,S) over a
complete block and write the posterior scalar recursion at accelerometer
epochs as

    a_(k+1) = A_k a_k + B_k y_k,                            (JG1)

where y_k is the exact conditioned physical accelerometer source
(physical acceleration + retained AG/BA/sensor term), and A_k,B_k are
generated by the literal OU prediction, every intervening S=0 correction,
AW covariance sync and the conditioned accelerometer update.  For a word W,

    a_N = A_W a_0 + sum_j G_(W,j) y_j.                      (JG2)

The dangerous quantity is the signed/DC induced gain

    G_DC(W)
      = sup_(bounded-primitive y, nonzero)
          |mean_W a| / source_norm_W(y),                    (JG3)

with the physical-acceleration component measured by its bounded primitive
rather than L-infinity amplitude.  A useful cycle certificate is

    rho_DC = sup_reachable |A_W| < 1
    and
    G_DC * declared_source_budget < g sigma_w.              (JG4)

This is the mathematically coherent version of the requested
OU+accelerometer+S minimization.

### Endpoint behavior on the coupled manifold

Although a full source-uniform G_DC still requires the reachable conditioned
Riccati map, the coupled endpoints show why no single mechanism's global
minimum is useful.

At tau=.02, the 0.1-s homogeneous OU multiplier is

    exp(-.1/.02)=exp(-5)=0.006737946999,

so before feedback a DC AW component loses 99.3262% per sync interval.  The
S-only kappa is minimal here because the rapidly dying AW has little time to
integrate into S.

At tau=12, the OU multiplier is

    exp(-.1/12)=0.9917012926,

so OU leakage is weak.  This is precisely the region where the integrated
chain persists and the S/accelerometer feedback must carry the certificate.
The shipping cadence is clamped at .15 s for tau>=11 s, so a four-S-row block
spans at most .6 s.  A 17-s interval contains at least 28 disjoint four-row
blocks, while AW covariance maintenance occurs about every .1 s.

The endpoint mechanisms therefore complement rather than add independently.

### What can and cannot be concluded numerically now

The number kappa_S,min is rigorous but proves that the S-only route is
quantitatively useless.  The OU endpoint numbers are also exact.  They do NOT
supply a rigorous numerical gamma_cycle because the conditioned
accelerometer gain depends on the reachable joint covariance.  The PSD AW
sync gives the exact Schur increment

    P_(a|l,s)^new=P_(a|l,s)^old+Delta,

but Delta may vanish in a floor-inactive direction.  In that case the
conditional AW component is either already informative to the accelerometer
or encoded in l and exposed by S; quantifying that tradeoff is the remaining
Riccati minimization.

Thus the next proof object is a 1-D conditioned Riccati/mean comparison over
the compact coupled tuple, not a sum of separately minimized dissipations.
A valid explicit lower constant must be obtained from the joint map JG1/JG2.


## 17. Scalar conditional variance is not a closed Riccati state

The proposed reduction to p=P_(a|l,s) alone is false for the literal
integrated-OU process.  After conditioning on s=(theta,b_g,b_a), one axis has

    x=(v,p,S,a),   P_c in S_+^4.

Prediction is

    P^- = Phi P^+ Phi^T + Q_d(tau,h,sigma^2),               (RC1)

and shipping Q_d has nonzero q_va,q_pa,q_Sa and the corresponding integrated
cross terms.  In particular q_Sa is explicitly nonzero.  Therefore the next
Schur complement

    p_new=P_aa-P_al P_ll^dag P_la                           (RC2)

depends on the complete prior conditional covariance, not only on p_old.
An S correction

    P^+ = P^- - P^- e_S
          (e_S^T P^- e_S+r_S^2)^-1 e_S^T P^-               (RC3)

also changes P_aa through P_aS^2/(P_SS+r_S^2), again requiring the cross
covariance.  Hence two shipping-admissible conditional covariances with the
same p can have different p_new.  A one-dimensional Riccati recurrence would
reintroduce an unjustified covariance relaxation.

The smallest exact one-axis conditioned Riccati state is the symmetric 4x4
P_c for (v,p,S,a), i.e. ten scalar entries.  Its maps are explicit:

prediction RC1;
S update RC3;
accelerometer update, after conditioning on AG/BA,

    P^+ = P^- - P^- e_a
          (e_a^T P^- e_a+r_a)^-1 e_a^T P^-;                (RC4)

and AW sync

    P^+ = P^- + delta e_a e_a^T,                            (RC5)

for the conditional component supplied by the PSD floor.  Mean maps use the
same gains.

Thus the exact joint DC optimization is a 4-state/10-covariance scalar-axis
problem, not a 21-state problem but not 1-D either.  Any proof claiming a
closed p recurrence must first establish an invariant relation expressing all
P_al cross terms as functions of p; shipping provides no such identity.


### 17.1 Exact four-state mean/covariance cycle

For one conditioned world direction define x=(v,p,S,a).  At each operation
the homogeneous mean and covariance obey the SAME linear map:

prediction:
    x<-Phi x,
    P<-Phi P Phi^T+Q;

S correction:
    A_S=I-K_S e_S^T,
    x<-A_S x,
    P<-A_S P A_S^T+K_S r_S^2 K_S^T;

accelerometer correction:
    A_a=I-K_a e_a^T,
    x<-A_a x + K_a y,
    P<-A_a P A_a^T+K_a r_a K_a^T;

sync:
    x unchanged,
    P<-P+delta e_a e_a^T.

Thus for a literal cycle C the homogeneous DC retention is

    M_C = product_chronological A_i Phi_i,                  (RC6)

and the physical/source transfer is the corresponding Duhamel sum.  The
quantity relevant to a constant transverse AW mode is

    rho_a(C)=|e_a^T M_C e_a|                                (RC7)

only when the root has v=p=S=0; for arbitrary carried LIN roots use the
augmented readout norm

    rho_read(C)
      = || e_a^T M_C P0^(1/2) || /
        sqrt(e_a^T P0 e_a),                                 (RC8)

or the covariance-weighted operator norm.  This is the exact composed object
that replaces the nonclosed scalar p recurrence.

Every measurement operation is nonexpansive in its pre/post covariance
metric, while prediction/process and sync change the metric through PSD
increments.  Therefore a source-uniform strict rho_read<1 requires a
strict-information event on every nonzero homogeneous direction over the
chosen word.  The four-row S Chebyshev lemma supplies this for the LIN
subspace; the accelerometer supplies direct a information.  This recovers
qualitative strict contraction without G0.

A useful explicit numerical rho still requires lower-bounding the joint
information Gramian in a proof-scaled 4-state coordinate.  The repository's
earlier raw OU controllability Gershgorin bound failed from conditioning, so
the correct next representation is the scaled factor/LDLT already identified
in the proof ledger, now only 4x4 per axis.


## 18. Full 4x4 information floor: valid but quantitatively non-closing

A correction to the proposed target is necessary.  The quantity

    lambda_min(I_tilde)

depends on the arbitrary proof coordinate scaling.  The contraction statement

    rho <= (1+iota)^(-1/2)

is coordinate-invariant only for the PRIOR-WHITENED information

    J_C = P0^(1/2) I_C P0^(1/2),                            (IF1)

equivalently for a Loewner comparison

    I_C >= iota P0^-1.                                      (IF2)

Thus the desired constant is

    iota0=inf_reachable lambda_min(J_C).                     (IF3)

A single .1-s AW-sync cycle cannot have a uniform 4-state information floor:
at T_S=.15 it may contain no S row, while conditioned accelerometer rows
observe only AW.  Use a four-S-service superword instead (duration <=.6 s),
or a longer selected-row word for better conditioning.

The existing exact recurring root covariance certificate uses fixed LIN proof
scales

    D=diag(2.4,18,132,4)

and supplies the 4x4 lower covariance matrix whose LDL pivots are

    4.2027825893e-8,
    3.9268873816e-9,
    9.3267257619e-11,
    1.9601549308e-9.

Its actual smallest eigenvalue is approximately

    lambda_min(L_root)=4.15919906e-11.                      (IF4)

This is a proved Loewner lower matrix; the decimal eigenvalue is diagnostic,
while a rational LDL/norm bound may be used for formal promotion.

For four CONSECUTIVE S rows at the coupled short-tau corner
tau=.02,T_S=.005,r_S=.15, adding the guaranteed accelerometer rows removes
the tiny AW Schur direction.  In the same fixed proof coordinates the
combined information eigenvalues are approximately

    1.5999956e-7, 1.7999998, 1.7925380e2, 3.0976090e6.      (IF5)

The small eigenvalue is now a neutral-integrator direction, not AW.  The
coordinate-invariant crude product comparison therefore gives only

    iota0 >= lambda_min(L_root) lambda_min(I_C)
           ~= 6.65e-18,                                    (IF6)

for this corner.  The corresponding contraction bound is

    rho <= (1+iota0)^(-1/2)
         = 1-O(3.3e-18),                                    (IF7)

which is useless on a 17-s horizon.

Selecting four S rows spread through a longer superword improves raw
conditioning, but does not repair the many-orders-of-magnitude loss caused by
the generic root covariance lower certificate.  The latter was designed only
to prove coercivity of the full stability energy, not a sharp AW reader
metric.

Therefore the full 4x4 minimum-eigenvalue route is classified as a
QUANTITATIVE RELAXATION FAILURE, not a physical counterexample.  It proves
strict information but cannot close the nominal-AW bridge.

### Reader-specific Schur information is much stronger

For the AW readout, eliminate (v,p,S) from the information matrix instead of
taking its smallest eigenvalue.  Accelerometer rows add only positive AW
information, so their contribution survives the Schur complement unchanged.
With proof AW scale d_a=4, worst accepted-sample spacing h_acc=.006,
tau>=.02 and R_acc<=0.3010398645^2, the root-AW information from the
conditioned accelerometer sequence has the analytical lower bound

    I_a,acc
      >= (d_a^2/R_acc,max)
          sum_(j>=1) exp(-2 j h_acc/tau_min)
      = (16/R_acc,max)
          exp(-.6)/(1-exp(-.6))
      ~= 214.7520821.                                       (IF8)

Using a finite superword only truncates the geometric series by its explicit
positive tail.  This number is independent of S-row conditioning and shows
that the dangerous AW coordinate itself is strongly read by the conditioned
accelerometer channel.

Combining IF8 with the GENERIC root conditional-AW lower pivot
1.9601549308e-9 still gives only about

    4.21e-7                                                (IF9)

of prior-whitened AW information, because that generic covariance floor is
extremely pessimistic.  But this is eleven orders of magnitude better than
the full 4x4 eigenvalue comparison and identifies the correct next object:
a reader-specific lower bound on the conditional AW variance supplied by the
.1-s PSD synchronization/previous-cycle information, not a full LIN
covariance eigenvalue.

Conclusion: the requested full 4x4 LDL/information certificate can establish
iota0>0, but its source-uniform value from existing covariance certificates is
far too small.  Do not use it to claim closure.  The nominal-AW bridge must
use the AW Schur/readout metric, where direct accelerometer information is
O(10^2) in the fixed proof coordinates, together with the exact sync
conditional-variance increment and the S-chain only for the floor-inactive
encoded branch.


## 18. No standalone post-sync conditional-AW variance floor

The default sync computes

    Delta = Pi_+(Sigma_aw_stat - P_aa)

from the MARGINAL AW block and adds Delta only to P_aa.  Therefore

    P_(a|l,s)^new = P_(a|l,s)^old + Delta                  (CF1)

is exact, but it does NOT imply a positive uniform lower bound on
P_(a|l,s)^new.

Counterfamily (scalar AW versus one encoded state z): for any stationary
target Sigma>0 and epsilon>0,

    P = [[Sigma, sqrt(Sigma(Sigma-epsilon))],
         [sqrt(Sigma(Sigma-epsilon)), Sigma]]

is PSD, has P_aa=Sigma so Delta=0, but

    P_(a|z)=epsilon.                                       (CF2)

Let epsilon->0.  Thus

    inf P_(a|l,s)^post-sync = 0                             (CF3)

over algebraically admissible covariances even with the marginal exactly at
the stationary target.  A positive reader-specific variance floor cannot be
proved from the sync policy alone.

This is not a shipping counterexample.  It proves that the desired certificate
must include the information carried by the state that encodes AW.

### Correct reader-specific coercive quantity

After conditioning on s=(theta,bg,ba), partition one-axis LIN covariance as

    P_c = [[P_ll, c],
           [c^T, P_aa]],  l=(v,p,S).

Let

    p_c = P_aa-c^T P_ll^dag c >=0,                          (CF4)

and define the best linear predictor coefficient

    L = c^T P_ll^dag,   a = L l + epsilon,
    Var(epsilon)=p_c.                                      (CF5)

Normalize a unit AW readout.  The fraction represented by epsilon is handled
by direct conditioned accelerometer information; the fraction represented by
L l is handled by the chronological S-chain information.  Therefore the
coercive object is not p_c but the infimal decomposition energy

    C_cycle =
      inf_(epsilon,Ll: epsilon+Ll=1)
        [ I_acc |epsilon|^2 + I_S(Ll) ],                    (CF6)

with I_acc the complete conditioned accelerometer information over the word
and I_S the exact S-row action of the encoded LIN predictor.  In scalar
relaxation, if I_S >= kappa_enc |Ll|^2,

    C_cycle >= I_acc kappa_enc/(I_acc+kappa_enc).           (CF7)

This harmonic-mean form is the exact completed-square minimum over how an AW
mode splits between fresh conditional variance and encoded LIN state.

The sync strengthens CF6 by increasing only the epsilon/conditional component
whenever Delta>0; when Delta=0, the marginal is already at target and any
small p_c necessarily puts almost all of the AW mode into the encoded branch.
No arbitrary cross-covariance cancellation is left uncharged.

The remaining quantitative problem is therefore to lower-bound kappa_enc for
the SPECIFIC predictor subspace reachable from the conditional covariance,
not the global S-only kappa_S,min.  The global value 1.01e-14 permits an
arbitrary (v,p,S) combination chosen solely to hide AW at four S rows.  A
shipping covariance predictor L=c^T P_ll^dag is generated by the same OU
process Q and previous corrections; it is not arbitrary.  Bounding its
reachable coefficient cone is the next required covariance-reachability
lemma.

The large I_acc>=214.752 remains useful: once kappa_enc is established,
CF7 is essentially kappa_enc whenever kappa_enc << I_acc.  Thus the bottleneck
is now entirely the reachable encoded-LIN predictor, not direct AW
accelerometer information.


### 18.1 Fresh OU-process predictor has a sign restriction

For one prediction interval started from deterministic LIN state, the fresh
OU process covariance is the Gram covariance of a scalar white-noise input
passed through positive kernels.  With state ordering l=(v,p,S),a, all
same-axis fresh cross covariances

    Q_av, Q_ap, Q_aS

are nonnegative for h>0,tau>0: their integral representations are products of
the positive OU kernel with its positive first/second/third integrals.
Likewise Q_ll has nonnegative entries.

Thus the fresh-process regression of a on l belongs to a restricted cone; it
is not the arbitrary coefficient vector used in the global four-row Schur
minimum.  Measurement corrections can change the carried regression, but
S=0 corrections are themselves generated by the positive S column and
accelerometer corrections observe a directly.  This strongly suggests the
polynomial-cancellation predictor attaining kappa_S,min is not reachable.

A proof still needs an invariant cone or a signed total-positivity statement
for the Riccati maps.  Merely observing positive Q entries is insufficient:
matrix inversion in L=c^T P_ll^-1 can change coefficient signs.  The next
analytical target is therefore a total-positivity/variation-diminishing
invariant for the conditioned one-axis covariance factors, which would
exclude the alternating predictor needed to cancel four S rows.


## 19. Reachable regression cone: exact measurement invariances

After conditioning on s=(theta,bg,ba), let the one-axis covariance of
(l,a), l=(v,p,S), be positive definite and write its precision as

    J=P^-1=[[J_ll,J_la],[J_al,J_aa]].

The covariance regression of a on l is exactly

    L=P_al P_ll^-1 = -J_aa^-1 J_al.                         (RG1)

This precision representation makes the measurement maps transparent.

### S=0 update

A noisy S measurement has H=[e_S^T,0], so in information form

    J^+ = J^- + diag_l(e_S e_S^T/r_S^2).                   (RG2)

The a-row J_al and scalar J_aa are unchanged.  Hence

    L^+ = L^-.                                              (RG3)

Thus every S=0 correction leaves the regression of a on the FULL carried
l=(v,p,S) exactly invariant.  S changes the distribution/covariance of l and
therefore the action of the encoded component, but it does not rotate L.

### Accelerometer update

A conditioned accelerometer observation of a has H=[0,1].  Therefore

    J_aa^+ = J_aa^- + 1/r_a,
    J_al^+ = J_al^-.

Hence

    L^+ = alpha_a L^-,
    alpha_a = J_aa^-/(J_aa^-+1/r_a) in (0,1).              (RG4)

An accelerometer correction scales L toward zero by a positive scalar.  It
cannot change the sign pattern or projective direction of L.

### AW sync

Adding delta to the covariance aa block leaves P_ll and P_al unchanged, so

    L^+ = L^-                                               (RG5)

exactly.  The sync increases only the conditional residual variance
p_c=P_(a|l) by delta.

Therefore, among the correction/sync operations,

    S update:       L unchanged,
    AW sync:        L unchanged,
    accelerometer:  L -> positive scalar * L.               (RG6)

Only OU prediction/process injection can rotate the projective regression
direction.  This is a major reduction of the reachable-cone problem.

### Prediction formula

Partition the exact transition as

    l^- = A l^+ + b a^+ + w_l,
    a^- = phi a^+ + w_a,

with

    A=[[1,0,0],[h,1,0],[h^2/2,h,1]],
    b=[phi_va,phi_pa,phi_Sa]^T,
    phi=exp(-h/tau),

and fresh process covariance
Q=[[Q_ll,q_la],[q_al,q_aa]].

Given P^+ blocks (X=P_ll,c=P_la,s=P_aa),

    X^- = A X A^T + A c b^T + b c^T A^T
          + s b b^T + Q_ll,                                (RG7)

    c^- = phi(A c+b s)+q_la,                               (RG8)

    s^- = phi^2 s+q_aa.                                    (RG9)

Therefore

    L^- = c^-T (X^-)^-1.                                   (RG10)

RG7--RG10 are the ONLY projective-rotation map that needs a cone invariant.

The fresh-process vectors b and q_la have strictly positive entries for
h,tau>0.  A and Q_ll are generated by nested positive integration kernels.
This makes prediction a positive-kernel covariance transformation, but
positivity of c alone does not immediately imply positivity of L because
X^-1 can have alternating signs.

The correct object for a cone proof is the PRECISION cross row
q=-J_al=J_aa L.  Measurements either leave q fixed (S), leave its direction
fixed while changing J_aa (accelerometer), or covariance-sync changes q only
through the scalar conditional variance update while preserving L.  Prediction
is the sole nontrivial map.

### Immediate consequence

Any alternating regression direction capable of realizing the unrestricted
four-row polynomial cancellation must be CREATED by OU prediction.  It cannot
be manufactured by repeated S corrections, accelerometer corrections, or AW
sync.  Thus a source-uniform cone proof may be reduced to one prediction step
acting on the post-accelerometer cone, followed by a positive scalar
contraction.

This reduces the remaining invariant problem from the full Riccati recursion
to the sign-regularity of RG7--RG10 under the integrated-OU prediction.


## 20. Proposed three-coefficient cone induction is not closed

Write the conditioned covariance in regression form

    a=L l+epsilon,  Cov(l)=X,  Var(epsilon)=p_c,
    Cov(l,epsilon)=0.

Then c=X L^T and s=L X L^T+p_c.  Substitution into the literal prediction
gives, with T=A+bL,

    X^- = T X T^T + p_c b b^T + Q_ll,                      (RG14)
    c^- = phi[T X L^T+p_c b] + q_la,                       (RG15)
    L^- = c^-T (X^-)^-1.                                   (RG16)

Therefore L^- depends on X and p_c separately, not on L alone.  Two priors
with the same regression vector but different conditional residual variance,
or different X, generally produce different projective directions after one
prediction.  The proposed induction L in C => L^- in C is thus not a closed
shipping reachability statement.

The fresh-Q moment curve remains useful but cannot characterize carried
predictors by itself: prediction mixes transformed carried covariance,
the rank-one conditional term p_c b b^T, and fresh Q, and regression of a
matrix sum is not a scalar convex combination of component regressions.

The minimum exact regression-form state is (X,L,p_c): six+three+one=ten
scalars, algebraically equivalent to the conditioned 4x4 covariance.
A lower-dimensional proof now requires a stronger invariant relation between
X,L,p_c.

The natural candidate is not a rectangular L cone but a total-positive Gram
FACTOR invariant.  The fresh OU Q is the Gram matrix of the nested positive
kernels (k_v,k_p,k_S,k_a).  Prediction appends fresh kernel-factor columns
after transporting the carried factor; Kalman measurements are conditioning/
orthogonal-projection operations on that factor.  A variation-diminishing
factor invariant, if preserved by those projections, would exclude the
arbitrary four-row polynomial-cancellation covariance while remaining closed
under RG14--RG16.

Until that factor invariant is proved, the three-coefficient cone shortcut is
classified as a nonclosed-state relaxation and must not be promoted.


## 21. Total-positive Gram-factor invariant is not preserved by Kalman conditioning

The proposed factor invariant fails under noisy coordinate conditioning.
Let P=F F^T and let a scalar measurement have row h and variance R>0.  Put

    z=F^T h^T,  q=z^T z+R.

Then

    P^+ = F [I-z z^T/q] F^T.                               (TP1)

The bracket is SPD, but a square-root factor is a dense rank-one contraction
whose off-diagonal signs follow -sign(z_i z_j).  It is not a
total-positive multiplier in general.  Equivalently, the covariance update

    P^+=P-P h^T(hPh^T+R)^-1 hP                             (TP-two)

is a rank-one subtraction.  SPD is preserved, but non-principal minors are
differences of products and their signs are not preserved in dimension >=3.
The special S and accelerometer coordinate rows do not restore a general
minor-sign theorem for the full 4x4 covariance/factor.

Therefore strict total positivity, and likewise a global sign-regular minor
pattern, is not a closed invariant of the shipping Kalman cycle.  Prediction
and fresh OU factor appending are compatible with positive-kernel structure;
measurement conditioning is the breaking operation.

The weaker structures that DO survive are:

1. the exact sparse regression identities of section 19;
2. PSD/Loewner order;
3. information monotonicity,

    J^+=J+H^T R^-1 H.                                      (TP3)

The third property gives a better route to kappa_enc.  Do not propagate a
covariance cone.  Select a fixed subset of literal chronological measurement
rows over a superword, transport them to a common root with the exact OU
transition, and form the deterministic information design

    I_sel = sum_i Phi_i^T H_i^T R_i^-1 H_i Phi_i.           (TP4)

Partition the root as l=(v,p,S) and a.  The scalar AW information after
eliminating l is the Schur complement

    kappa_sel =
      I_aa-I_al I_ll^dag I_la.                             (TP5)

Adding any extra applied measurement row adds a PSD term to I_sel, and the
minimum residual characterization

    kappa_sel =
      min_z || R^-1/2 (H_a-H_l z) ||^2                     (TP6)

shows that adding rows cannot decrease kappa_sel.  Thus a lower bound proved
from selected S/accelerometer rows remains valid under the complete Kalman
conditioning chronology; no factor-sign invariant is required.

This selected-row Schur certificate differs from the failed global 4x4
eigenvalue bound: it eliminates the neutral polynomial root directions FIRST
and targets only the AW reader.  It also differs from the unrestricted
S-only kappa: accelerometer rows directly observe AW and cannot be fitted by
the l polynomial columns, so they regularize the near-polynomial cancellation.

The next quantitative calculation is therefore explicit: choose four
well-separated guaranteed S rows over the shortest uniform superword and a
guaranteed subset of accelerometer rows; build TP4 using the exact
piecewise-tau transitions; analytically lower-bound TP6 over the coupled
tau/rS/cadence tuple.  This scalar residual norm is the correct
conditioning-invariant kappa_enc certificate.


## 22. Finite selected accelerometer certificate: use the first guaranteed row

The previous infinite geometric sum I_acc>=214.752 treated later
accelerometer rows as open-loop root observations H Phi(t).  That is not a
rigorous lower bound on information about the same carried root after earlier
Kalman corrections: the later residual sensitivity is H M_(k:0), with the
closed-loop mean/covariance sensitivity including those earlier corrections.
Using open-loop Phi for every later row can over-count root information.

A completely safe finite certificate uses only the FIRST guaranteed
accelerometer row after the chosen root.  No earlier accelerometer correction
can have reduced its root sensitivity.  Scheduled finite accelerometer
updates are mathematically applied on the retained class by the existing
LDLT/no-invalid-input argument.

In the fixed proof AW coordinate scale d_a=4, the root-AW coefficient at the
first row is

    h_a = d_a exp(-h/tau).

With h<=h_max=.006 s, tau>=tau_min=.02 s and
R_acc<=R_acc,max=(.3010398645)^2, the scalar selected information is

    I_acc,1
      >= d_a^2 exp(-2 h_max/tau_min)/R_acc,max
       = 16 exp(-.6)/(.3010398645)^2
       = 96.89364056.                                      (FA1)

This is finite, chronological and noncircular.  It uses one guaranteed row,
not an infinite sum or a fitted carried statistic.

If the selected design also includes four guaranteed S rows, eliminate the
neutral root l=(v,p,S) by the least-squares Schur complement.  The first
accelerometer row has zero l columns in the conditioned one-axis open-loop
root model, so its residual cannot be canceled by the l fit.  Therefore

    kappa_sel >= I_acc,1 >= 96.89364056                     (FA2)

for that selected root design.

CAVEAT: FA2 is a deterministic selected-row design bound for the conditioned
root model.  To convert it into the literal FILTER closed-loop root-loss
coefficient, the selected residual must be represented with the same-history
innovation/source factor used by the complete-word information identity.
Only the first accelerometer row is immune to prior accelerometer
conditioning; later rows require closed-loop sensitivities.  Do not restore
the 214.752 infinite-series value without that calculation.

Combining FA1 with a prior conditional-AW variance p gives the one-row
dimensionless information

    iota_1 >= p * 96.89364056.                              (FA3)

The generic root p floor ~1.96e-9 still makes FA3 weak.  Thus FA1 confirms
again that measurement information is ample; the quantitative bottleneck is
the reader-specific prior decomposition between conditional AW variance and
the encoded LIN component, not the number of accelerometer rows.


## 23. Same-history conditional decomposition alone does not improve kappa_enc

For any prescribed regression L_*, choose X=X^T>0 and p_c>0 and define

    c=X L_*^T,
    s=L_* X L_*^T+p_c.

Then P=[[X,c],[c^T,s]] is positive definite with Schur complement p_c and

    P_a,l P_ll^-1=L_*.

Thus even the exact four-S cancellation regression is compatible with an SPD
same-history conditional covariance, with arbitrarily small p_c.  Covariance
algebra alone cannot improve the global kappa_S,min=1.0104660589e-14.

A stronger kappa_enc must therefore use SHIPPING REACHABILITY of P from the
literal construction through prediction/correction/sync maps.  Generic PSD,
marginal floors, upper/lower covariance boxes and Schur identities are
insufficient unless they encode that reachability.

A promising closed reachability variable is information age.  Decompose the
pre-measurement covariance factor into

    P = P_carried + P_fresh,

where P_fresh is the sum of independent OU process factors injected since a
chosen previous accelerometer correction.  Prediction transports P_carried
and appends a fresh independent factor, so this decomposition is exact before
the next measurement.  The fresh factor has the known OU positive-kernel
Gram geometry.  If repeated accelerometer corrections uniformly contract the
AW content of P_carried while every prediction injects a nonzero fresh AW
factor, then a source-uniform lower fraction

    P_fresh,aa / P_aa >= eta_fresh>0

(or its conditional/readout analogue) would force every near-deterministic
encoded AW predictor to contain a nonzero fresh-OU component.  Its selected-S
residual can then be bounded by the explicit OU kernel rather than the
arbitrary polynomial fit.

This is a genuine reachability statement and is the next required step for
any kappa_enc substantially above 1e-14.


### 23.1 Information age must be tracked in source-factor history

Do not propagate P_carried and P_fresh through the nonlinear Riccati map as
separate posterior covariances: C_H(P_a+P_b) != C_H(P_a)+C_H(P_b).

The exact age decomposition lives instead in the linear source-factor/smoother
representation.  For independent process factor B_j injected at prediction j,
its contribution to a later homogeneous/source reader is

    G_(k,j)=M_(k:j+1) B_j,

where M contains the literal subsequent closed-loop correction maps.  The
covariance/source action is the square sum of these chronological factors.
Source age is therefore an exact label on columns, preserved through every
linearized correction.

A useful fresh-information certificate must lower-bound the AW/readout norm
of columns with age <=T_age relative to the total relevant reader norm:

    sum_{j:k-j<=T_age} ||e_a^T G_(k,j)||^2
      >= eta_fresh
         sum_j ||e_a^T G_(k,j)||^2.                         (AGE1)

Unlike an additive posterior-covariance split, AGE1 is meaningful under
Kalman conditioning.  It is also exactly the representation already used by
the complete-word square-summed reader.

The remaining analytical task is to prove eta_fresh>0 from:
(i) nonzero OU process injection each step;
(ii) finite sample spacing;
(iii) the finite conditional accelerometer correction;
(iv) the .1-s AW sync factor; and
(v) the existing covariance upper bound, which prevents an unbounded ancient
factor from dominating forever.

Once AGE1 holds, the fresh columns have explicit OU kernel geometry and give
a nonzero reachable encoded-S floor.  This is the first route identified here
that both encodes shipping reachability and survives the nonlinear Riccati
conditioning exactly.


## 24. Source-age floor is a restricted joint information-ratio problem

A terminal-survival-only lower bound for a fresh OU source is impossible:
later measurements may make its terminal sensitivity arbitrarily small.
The source is then observed rather than destroyed.  Conversely a
measurement-information-only floor can be tiny in directions removed mainly
by terminal forgetting.  The proof ledger already records this exact
phenomenon for complete words: separated information-only and forgetting-only
margins failed by many orders, while the joint smoother/information-ratio
matrix remained effective.

Therefore q_age must use the SAME joint matrix mechanism, restricted to the
fresh OU source subspace.

Let B_j be a canonical full-rank factor of the one-axis fresh OU process
covariance Q_j (for example Q_j^(1/2)).  Let A_age be the prior/source energy
matrix for that source coordinate and J_age the exact finite-horizon joint
matrix containing BOTH corrected measurement loss and terminal forgetting,
constructed by the existing smoother identity over T_age<=.6 s.  Define

    kappa_age =
      lambda_max( J_age^-1 A_age )                          (AGE1)

on range(B_j), with the usual generalized-eigenvalue interpretation.  Theorem
D / the information-ratio lemma then gives the fresh-source contraction

    rho_age <= tanh( log(kappa_age)/4 ) < 1                 (AGE2)

whenever kappa_age is finite.  Equivalently the dimensionless fresh
coercivity may be written

    q_age =
      lambda_min( A_age^-1/2 J_age A_age^-1/2 ) >0.         (AGE3)

This definition is invariant to the arbitrary factorization B_j and does not
multiply by the generic root covariance floor.

Qualitative positivity follows from the conditioned four-state zero-action
argument: a nonzero fresh OU source with zero corrected accelerometer loss,
zero four-row S loss and zero terminal forgetting would lie in the complete
fresh-source nullspace, which is trivial.

But an EXPLICIT useful lower q_age cannot be obtained by multiplying the
first-row accelerometer number 96.89 with the S-only floor.  The repository
has already demonstrated that information and forgetting can act in different
eigendirections; their generalized eigenvalues must be combined before taking
a scalar minimum.

The correct numerical certificate is thus a 4x4 RESTRICTED version of the
existing Theorem-D matrix calculation:
1. inject one fresh one-axis OU Q_j;
2. retain its four source coordinates only;
3. propagate its cross-covariance through the literal closed-loop chronology
   for <=.6 s;
4. accumulate exact smoother correction losses;
5. include terminal conditional forgetting;
6. form A_age^-1/2 J_age A_age^-1/2;
7. prove its LDLT/Schur floor uniformly over the coupled tau,sigma,rS,T_S,dt
   and scheduler phase.

This is noncircular and substantially smaller than the old full-root 4x4
information certificate: no carried root covariance appears.  It is also not
a new proof architecture; it is Theorem D restricted to a fresh process
factor already present in the complete-word joint factorization.


## 24. Source-age q0: dependence audit and noncircular restriction

The proposed q0 is NOT determined by the fresh 4-state OU factor and tuner
tuple if accelerometer smoother losses are included.  For a fresh source
coordinate z,

    Delta J_k=C_k H_k^T S_k^-1 H_k C_k^T,

and S_k=H_k P_k^- H_k^T+R_k contains the carried full covariance.  An
accelerometer S_k includes attitude/gyro/BA covariance and cross terms.
The existing recurring nuisance upper theorem explicitly does not supply an
AG upper covariance.  Bounding S_k from above through G0/AG control would
reintroduce the circularity this route was intended to avoid.

Therefore a noncircular source-age certificate must initially OMIT
accelerometer smoother loss.  Keep only:
- four guaranteed S=0 corrected losses, whose innovation covariance depends
  on LIN and is bounded by the proved nuisance comparison;
- terminal LIN forgetting/conditional remainder from the same complete-word
  identity.

This restricted joint matrix is smaller than the full J_age, so any lower
bound remains valid when accelerometer/magnetic corrections are restored by
information monotonicity.

Let J_age,S+T be the exact joint information/forgetting matrix on the fresh
one-axis OU source coordinates using only S rows and terminal LIN remainder.
Define

    q0_ST =
      inf lambda_min(
        A_age^-1/2 J_age,S+T A_age^-1/2 ).                  (AGE4)

Then

    J_age,full >= J_age,S+T

in the joint minimum-action/information sense, hence q0_full>=q0_ST.

All denominators needed by AGE4 are now controlled by the source-uniform
nuisance covariance upper theorem and R_S<=10000 I.  No G0, attitude
covariance ceiling or nominal-AW premise is used.

Caution: the S-only instantaneous Schur floor 1e-14 suggests AGE4 may still
be numerically weak, but terminal forgetting acts in complementary
directions.  The proof ledger already shows that such information+forgetting
combinations can be much stronger than either scalar margin separately.
The correct next computation is therefore the restricted 4x4 Theorem-D
matrix for one fresh OU factor with S rows + terminal LIN forgetting, not an
accelerometer-inclusive Gramian.


## 25. Source-age q_ST is not the contraction/mean certificate

A final normalization audit shows that the proposed source-age coercivity is
not, by itself, the quantity needed for the nominal-AW theorem.

Whiten one freshly injected process coordinate so its prior source action is
A_age=I.  The exact smoother decomposition partitions that source's influence
between (i) information extracted by later measurements and (ii) terminal
conditional uncertainty/forgetting.  A lower bound on

    J_age = extracted_information + terminal_accounted_energy

therefore proves that the source is not in an unaccounted nullspace.  It does
NOT say that its contribution to the terminal/mean AW readout is small.
Indeed an exactly conserved source can have J_age=I with no useful mean
attenuation, while a strongly observed source can also have large J_age.

Theorem D avoids this ambiguity by comparing two terminal covariance
problems (known-root versus diffuse/root-uncertain) through

    kappa=lambda_max(Pi^-1 P_diff),

and then converting that relative diameter to contraction.  There is no
analogous implication

    lambda_min(J_age)>0  =>  nominal AW mean < threshold

without an additional readout/terminal comparison.

Thus q_ST is a valid detectability/accounting modulus but is not the missing
nominal-mean bound.  Computing a tiny or large q_ST cannot decide whether the
gravity-sized pathological mean is shipping-reachable.

This source-age branch is therefore removed from the critical path.  Retain:
- the exact zero-innovation 17-s incompatibility subcase;
- the complete covariance-weighted nominal-AW reader;
- the endogenous-innovation correction;
- the exact S passivity identities;
- the finite first-row accelerometer information as a local fact;
- the existing Theorem-D complete-word machinery.

The actual unresolved theorem remains a SIGNED readout/source bound for the
nominal AW mean (or an equivalent closed-loop reachability contradiction).
It must compare the terminal/multi-time AW readout directly with the allowed
physical/bias/sensor source set.  Detectability of fresh process factors is
insufficient.

Consequently no numerical q_ST is promoted as a stability margin.  Doing so
would be another formulation error rather than progress on the theorem.


## 26. Return to the signed AW balance: exact gain-weighted Abel identity

The critical path returns to the literal world-frame AW loop.  Let
e_k=a_hat_k-a_k and let Gamma_k=K_aw,k Rhat_k at accepted accelerometer
corrections.  With xi_k collecting actual S/magnetic AW mean corrections, the
exact identity already proved in the repository is

    sum_acc Gamma_k(e_k-eta_k)
      = e_0-e_N + sum xi_k
        -sum_pred[(1-phi_k)a_hat_k + Delta a_k].             (SA1)

Because sum_pred Delta a_k=a_N-a_0 on the same physical history, no
independent physical-increment source remains.  Substitute a_hat=e+a and
rearrange.  For any fixed world row u^T,

    sum_pred (1-phi_k) u^T a_k
      + sum_acc u^T Gamma_k a_k
    = endpoint/error terms + S/mag terms + eta terms.        (SA2)

Thus the physical acceleration enters through one signed coefficient sequence

    W_k = (1-phi_k) I + 1_acc(k) Gamma_k                    (SA3)

after aligning prediction and correction epochs exactly (with zero Gamma on
non-acc epochs and the literal chronology used for staggered events).

For a scalar projected physical acceleration a_k=dv/dt sampled at h_k, define

    w_k = u^T W_k / h_k.

The exact first Abel summation gives

    sum_k u^T W_k a_k
      = w_N v_N-w_1 v_0
        -sum_(k<N)(w_(k+1)-w_k)v_k
        -sum_k w_k q_k,                                    (SA4)

where the jerk sampling remainder satisfies

    |q_k| <= J_max h_k^2/2.

Hence

    |sum u^T W_k a_k|
      <= Vmax[|w_1|+|w_N|+TV(w)]
         +(Jmax/2) sum h_k ||W_k||.                         (SA5)

This is the exact noncircular place where bounded physical velocity enters.
Unlike the earlier arbitrary beta reader, W_k is LOCAL:

    W_k=(1-phi_k)I+Gamma_k,

and all S effects remain on the opposite side through the signed sum xi_k.
The problem is therefore narrower than uniform variation of a backward
multi-time beta adjoint.

### S-chain must be combined before bounding TV(W)

The S correction contribution is

    xi_k^S = -K_aw,S,k S_k^-.

Do not bound sum xi^S separately.  The same covariance chronology that makes
Gamma_k periodic generates K_aw,S,k and S_k.  Over one scheduler block, move
sum xi^S to the left of SA2 and define the effective local coefficient by
eliminating the homogeneous (v,p,S) response with the exact S-chain identity.
Call the resulting physical-acceleration coefficient W_eff,k.

Then the desired signed reader bound is

    |sum u^T W_eff,k a_k|
      <= Vmax D1(W_eff/h)
         +(Jmax/2)sum h_k||W_eff,k||,                        (SA6)

plus root, BA/attitude/sensor and magnetic defects already isolated in the
complete reader.

This identifies the genuine quantitative target:

    D1_eff =
      sup_shipping [
        |w_eff,1|+|w_eff,N|
        +sum||w_eff,k+1-w_eff,k||
      ].                                                     (SA7)

The carried 0.371 m/s2 sync-locked rectification is evidence about SA6 but is
not a bound.

The advantage over the abandoned beta-TV route is structural: W_eff is
generated by a LOCAL OU+accelerometer+S block, not by a backward adjoint over
the whole readout.  Tuner/sync jumps therefore enter only at their actual
local events, and S-chain cancellation is performed before variation is
taken.

A useful theorem is

    Vmax D1_eff
      +(Jmax/2) sup sum h||W_eff||
      + C_root+C_BA+C_att+C_sensor+C_mag
        < g sigma_w.                                       (SA8)

No covariance detectability surrogate implies SA8; it must be bounded
directly from the reachable local gain chronology.


## 27. Local W_eff block: exact form and covariance dependence

The signed AW balance gives the local pre-elimination coefficient

    W_k=(1-phi_k)I+Gamma_k,
    Gamma_k=K_aw,k Rhat_k.

To eliminate S mean feedback over a scheduler block, retain the block root
z0=(v,p,S,a) and write the literal LIN mean recursion

    z_(i+1)=A_i z_i+B_i y_i,                                (WE1)

where A_i is prediction, S=0 correction, or the LIN part of an accelerometer
correction, and y_i is the same-history physical/nuisance accelerometer
source.  For a scalar transverse readout of the AW balance, backward
elimination over the LOCAL block gives

    R_block = lambda_0^T z0 + sum_i w_eff,i y_i.            (WE2)

Thus W_eff exists as a local block coefficient only together with a single
carried block-root term lambda_0^T z0.  Hiding the root inside W_eff would
make the coefficient nonlocal again.

For an S event specifically,

    A_S=I-K_S e_S^T,
    K_aS=P_aS(P_SS+R_S)^-1,                                (WE3)

so the AW correction is -K_aS S.  Therefore the eliminated coefficients
depend not only on tau,sigma,R_S,T_S but on the carried LIN covariance
through P_aS and P_SS.  Mean S-chain elimination does NOT remove this
covariance dependence.

This blocks any claim that D1_eff is a function only of the tuner tuple.

### Noncircular LIN gain bounds

Unlike accelerometer AG coupling, the S gain is purely LIN.  The recurring
nuisance upper comparison gives source-uniform principal bounds on
P_aa,P_SS.  PSD Cauchy gives

    ||P_aS|| <= sqrt(||P_aa|| ||P_SS||).                    (WE4)

Since P_SS+R_S >= R_S,

    ||K_aS||
      <= sqrt(||P_aa|| ||P_SS||)/lambda_min(R_S).           (WE5)

This is rigorous and noncircular but likely very coarse.

A sharper identity uses the Joseph decrement:

    K_aS (P_SS+R_S) K_aS^T
       = P_aa^- - P_aa^+ |_S >=0.                           (WE6)

Hence for any direction u,

    ||u^T K_aS||^2
      <= [u^T(P_aa^- - P_aa^+|_S)u]/lambda_min(R_S).        (WE7)

Summing WE7 over S events telescopes only after adding AW process/sync
replenishment between events.  This is the appropriate way to control local
gain variation without multiplying the enormous P_SS upper box by
1/R_S,min.

### Effective variation decomposition

For aligned physical sample epochs define

    w_eff,k = u^T W_eff,k/h_k.

Its first variation splits exactly into

    Delta w_eff
      = Delta[(1-phi)/h] u^T
        + Delta[ u^T Gamma/h ]
        + Delta w_S,elim.                                  (WE8)

The OU scalar part has an analytic derivative bound.  Put
f(h,tau)=(1-exp(-h/tau))/h.  For h in [.004,.006],
tau in [.02,12], f is positive and smooth; tuner tau changes only through
the applied smoothed chronology.  Thus

    |Delta f|
      <= L_h |Delta h| + L_tau |Delta tau|,                 (WE9)

with explicit sup derivatives on the compact rectangle.

The accelerometer and eliminated-S parts must be treated jointly through
their Joseph covariance decrements and the .1-s AW sync/process
replenishment.  Bounding Delta Gamma by independent covariance boxes repeats
the failed huge-TV relaxation.

The correct local block target is therefore a covariance-budget variation
inequality

    sum_(k in block) ||Delta[Gamma/h + W_S,elim/h]||
      <= C_block sqrt( sum measurement decrements
                       + sum AW process/sync replenishment ), (WE10)

followed by Cauchy--Schwarz across blocks.  All quantities in WE10 are LIN
covariance quantities; no G0 or AG ceiling is required for the S part.
The accelerometer Gamma still contains AG/BA cross covariance, so its
variation cannot be bounded from the LIN nuisance theorem alone.  It must
remain in the full signed accelerometer identity rather than be boxed.

Conclusion: the S elimination is useful, but it does not by itself produce a
tuner-only D1_eff.  The remaining noncircular local calculation is to prove
WE10 for the S contribution and OU term, while leaving Gamma inside the
endogenous accelerometer balance.


## 28. Explicit OU variation: raw D1_OU+S is quantitatively impossible

For the OU part of the signed coefficient,

    f(h,tau)=(1-exp(-h/tau))/h.

On h in [.004,.006], tau in [.02,12],

    f_min = 0.0833125035 1/s,
    f_max = 45.31731173 1/s.                               (OV1)

The derivatives are

    partial_h f =
      [(1+h/tau)exp(-h/tau)-1]/h^2,

    partial_tau f =
      -exp(-h/tau)/tau^2.                                  (OV2)

Their absolute suprema on the rectangle occur at the fast corner
h=.004,tau=.02 and are approximately

    L_h = 1095.194 1/s^2,
    L_tau = 2046.827 1/s^2.                                (OV3)

Thus tuner/sample variation can be bounded explicitly by

    TV(f) <= L_h sum|Delta h| + L_tau sum|Delta tau|,

with the actual smoother chronology used for the second term.

However the first-Abel charge includes endpoint terms even when f is
CONSTANT:

    D1_OU = |f_1|+|f_N|+TV(f) >= 2 f.                      (OV4)

At a representative h=.005,tau=.02,

    f=44.23984339 1/s,

so bounded physical velocity alone gives

    Vmax * 2f = 5.5*88.47968677
              = 486.6382772 m/s^2.                         (OV5)

This exceeds g sigma_w by over two orders of magnitude.  At the exact
rectangle maximum, 2 Vmax f_max is about 498.49 m/s^2.

Therefore NO refinement of the S contribution can make the raw
D1_OU+S certificate useful uniformly.  The failure is already present with
zero S gain and constant tuner parameters.

This is a structural relaxation failure, not evidence of an unstable
shipping trajectory.  Abel has separated

    sum_pred (1-phi_k) a_k

from the error endpoints in SA1.  For fast OU, (1-phi)/h is O(1/tau), so the
velocity-endpoint charge is huge; but the SAME OU prediction also contracts
the carried AW/error state by phi.  Bounding those two effects separately
destroys their cancellation.

### Correct regrouping

Return to SA1 before moving the OU physical term alone.  At a prediction,

    e^+ = phi e^- -(1-phi)a^- - Delta a
        = phi a_hat^- - a^+.                               (OV6)

Hence the combination

    endpoint error + sum_pred[(1-phi)a_hat+Delta a]

must be telescoped at the level of a_hat/e together.  For one prediction
with no intervening correction, the identity is exact and has no
1/tau-amplified physical-velocity endpoint.

Over a block, define the OU-propagated error endpoint

    E_OU(block)
      = e_end - Phi_block e_start

and retain all accelerometer/S corrections as signed injections.  Variation
of constants gives

    E_OU
      = -sum_j Phi_(end<-j)
          [(1-phi_j)a_j+Delta a_j]
        + signed correction transport.                     (OV7)

The physical coefficient now contains the DECAYING future multiplier
Phi_(end<-j).  Summation by parts acts on

    beta_OU,j =
      Phi_(end<-j)(1-phi_j),

not on (1-phi_j) alone.  For constant h,tau,

    beta_OU,j=(1-phi) phi^(N-j),

whose total first variation is bounded independently of 1/tau:

    endpoint+TV <= 2(1-phi) <=2.                            (OV8)

More directly, the geometric coefficients telescope and their l1 mass is

    sum_j beta_OU,j = 1-phi^N <=1.                          (OV9)

This is the missing cancellation.

The S corrections must be transported by the SAME future OU multipliers
before their Joseph-budget/variation bound is taken.  Thus the useful local
coefficient is a forward-decayed block kernel, not W=(1-phi)+Gamma at the
same epoch.

Conclusion: an explicit finite D1_OU+S exists, but the requested raw local
D1 is provably useless (OV5).  The next valid signed calculation is the
forward-decayed OU+S block kernel beta_OUS, for which the OU component has
uniform l1 mass <=1 and first variation <=2 in the constant-parameter case.
Tuner variation can then be charged as a perturbation of a probability-like
decay kernel rather than as O(1/tau) endpoint variation.


## 29. Variable-phi forward OU kernel

Let T_j=product_(m=j..N) phi_m and T_(N+1)=1.  The forward-decayed physical OU coefficient is

    beta_j=(1-phi_j) product_(m=j+1..N)phi_m=T_(j+1)-T_j.

For arbitrary shipping 0<phi_j<1, T_j is nondecreasing, hence beta_j>=0 and

    sum beta_j=1-T_1<=1.

For augmented first variation D1=|beta_1|+|beta_N|+sum|Delta beta|, every nonnegative sequence satisfies D1<=2 sum beta, therefore

    D1(beta_OU)<=2(1-T_1)<=2.

This is exact for arbitrary time-varying tau/dt; no tuner smoothness bound is needed.

At an S event s, future OU transport gives Z_s=T_(s+1)K_aS,s. Joseph gives

    ||u'Z_s||^2 <= T_(s+1)^2 [u' DeltaP_aa,s^S u]/lambda_min(R_S,s).

Thus

    sum_s ||u'Z_s|| <= sqrt(H_S) sqrt(sum_s u'DeltaP_aa,s^S u),
    H_S=sum_s T_(s+1)^2/lambda_min(R_S,s).

Also D1(Z_S)<=2 sum_s||Z_s||.  A naive independent-extrema bound on H_S is too loose; the remaining scheduler/tuner problem is exactly to bound H_S under the coupled applied tau,R_S,T_S chronology.  The covariance factor is the S-induced AW Joseph decrement and can be paired with forward-decayed OU-process/sync replenishment.  No AG covariance enters.


## 30. Uniform coupled bound on the transported S chronology sum

Recall

    H_S=sum_s T_(s+1)^2/lambda_min(R_S,s),

where T_(s+1) is the future OU tail product.  The applied tuner tuple is
committed only at the adaptation cadence (~.1 s); inside each activation cell
tau,T_S,r_S are fixed.  The pseudo scheduler preserves elapsed phase when T_S
is retargeted.

Use cells of length L=.1 s.  In a cell with fixed applied tuple, a conservative
event-count bound is

    N_S <= 1+ceil(L/T_S),                                   (HS1)

where the extra one covers a phase-preserved event at the cell boundary.
The smallest S standard deviation is the Y factor .50 times the base r_S, so

    lambda_min(R_S) >= (.5 r_S)^2,   r_S>=.15.             (HS2)

Future squared OU transport across a complete cell is

    d(tau)=exp(-2L/tau).                                    (HS3)

Define the local worst charge

    A(tau,r_S)
      =[1+ceil(L/T_S(tau))]/(.5 r_S)^2.                    (HS4)

The backward chronology sum satisfies the scalar comparison

    H_i <= A_i+d_i H_(i+1).                                (HS5)

Therefore any

    M >= sup_reachable A/(1-d)                             (HS6)

is invariant under arbitrary sequences of applied cells: if H_(i+1)<=M then
H_i<=A_i+d_i M<=M.

Using the shipping coupled cadence

    T_S(tau)=clamp((.015/1.1)tau,.005,.15),

tau in [.02,12], and retaining the rigorous possibility r_S=.15 even at
large tau (SpectralMSE can remain on the base floor as physical wave RMS
approaches its 1e-6 guard), the maximum of HS6 is at

    tau=12 s, T_S=.15 s, r_S=.15 m*s.

There

    N_S<=2,
    A=2/.075^2=355.5555556,
    d=exp(-.2/12)=0.9834714538,

hence

    H_S <= 21511.60494.                                    (HS7)

This is a source-uniform analytical comparison for arbitrary coupled tuner
cell sequences.  It is conservative relative to the frozen-parameter
geometric sum (about 7200.37 at the same corner) because HS1 permits one
retarget/phase-boundary event in every .1-s cell.

### Quantitative consequence

Combining HS7 with the Joseph/Cauchy estimate gives

    D1(Z_S)
      <=2 sqrt(21511.605)
         sqrt(sum_s u^T DeltaP_aa,s^S u)
      ~=293.337 sqrt(S-decrement budget).                   (HS8)

Even with a sharp O(10) AW covariance/decrement budget, HS8 is far too large
for a g sigma_w mean certificate.  Thus the chronology sum H_S is now
rigorously finite and explicitly bounded, but Joseph + l1/Cauchy is a
demonstrated quantitative relaxation failure.

The failure is the same structural one seen earlier: S corrections are signed
and correlated through the integrated chain, while HS8 replaces them by the
sum of their norms.  A useful proof must keep the signed future-transported
S-chain combination before Cauchy/variation, just as the OU kernel had to be
kept in telescoping form before Abel.

Consequently:
- variable-phi OU is closed sharply: mass<=1, D1<=2;
- coupled S chronology is closed as H_S<=21511.605;
- but the generic Joseph-to-l1 conversion is unusable.

The remaining S calculation is to derive the SIGNED transported S-chain
kernel itself and seek a telescoping/divided-difference identity for it,
rather than bounding individual K_aS events.


## 31. Signed S jump identity and regression-ratio nonclosure

At an S=0 update,

    Delta a^S = -K_aS S^-,
    Delta S   = -k_SS S^-,

with

    K_aS=P_aS/(P_SS+R_S),
    k_SS=P_SS/(P_SS+R_S).

Hence, whenever P_SS>0,

    Delta a^S = C_S Delta S,
    C_S=P_aS/P_SS.                                         (SK1)

The pseudo-measurement denominator cancels exactly.  Future OU transport gives

    Xi_S=sum_s T_(s+1) C_S,s Delta S_s.                    (SK2)

Thus the signed S contribution is a discrete Stieltjes sum.  Summation by
parts moves differences onto A_s=T_(s+1)C_S,s and leaves S endpoints plus
the exact propagation gaps between S events.  Those gaps are generated by
the same (v,p,S,a) integrated chain and must be combined with the physical/OU
balance before norms.

An S update itself leaves C_S invariant because both P_aS and P_SS are
multiplied by R_S/(P_SS+R_S).

However C_S is NOT a closed scalar under the rest of shipping chronology.
For one OU prediction (one axis, ordering v,p,S,a), let
r_Srow=[h^2/2,h,1,phi_Sa].  Then

    P_aS^- =
      phi[ (h^2/2)P_av + h P_ap + P_aS + phi_Sa P_aa ]
      + q_aS,                                               (SK3)

while

    P_SS^- = r_Srow P r_Srow^T + q_SS.                     (SK4)

Therefore C_S^-=P_aS^-/P_SS^- depends on the full 4x4 LIN covariance.
An accelerometer correction also changes numerator and denominator by
different Schur products, so it does not preserve C_S.

Consequently the proposed scalar lemma

    TV[T_future P_aS/P_SS] <= C_reg

cannot be proved from a recurrence for C_S alone.  The exact signed jump
identity SK1 remains valuable, but a useful variation theorem must either:
(a) keep the full conditioned 4x4 LIN covariance in the coefficient history,
or
(b) avoid TV(C_S) entirely by combining SK2 with the propagation equations
before introducing C_S as a separate coefficient.

Route (b) is preferable.  Substitute Delta S_s=S_s^+-S_s^- directly into
the chronological S-state recursion and telescope the ACTUAL S jumps before
dividing by P_SS.  This keeps the regular product K_aS S^- intact and avoids
the artificial singular ratio when P_SS is small.

Thus the next signed derivation should use the pair
    (future-transported AW equation, future-transported S equation)
as a 2-row block and eliminate S jumps by block Gaussian elimination.
The elimination coefficient is computed at the block-matrix level, where
P_SS+R_S remains regular; do not form C_S eventwise.


## 32. Two-row AW/S service-block elimination

Take one complete S-service block from immediately after S service b to
immediately before service b+1.  Let x=(v,p,S,a).  Collapse every literal
operation inside the open block (OU predictions, accelerometer corrections,
syncs and their same-history source inputs) into

    x_(b+1)^- = F_b x_b^+ + G_b y_b,                        (BL1)

where y_b denotes the stacked physical/accelerometer/implementation sources.
No S pseudo-update is hidden in F_b.

At the terminal S service, with h_S=e_S^T,

    nu_b = 0-h_S x_(b+1)^- = -S_(b+1)^-,                   (BL2)

    x_(b+1)^+ = x_(b+1)^- + K_b nu_b
              = (I-K_b h_S)x_(b+1)^-.                      (BL3)

The two rows of interest are therefore

    a_(b+1)^+
      = e_a^T F_b x_b^+ + e_a^T G_b y_b + K_aS,b nu_b,     (BL4)

    S_(b+1)^+
      = e_S^T F_b x_b^+ + e_S^T G_b y_b + k_SS,b nu_b.     (BL5)

Together with BL2, this is a 2-row block with ONE internal scalar nu_b.
Do not divide BL4 by BL5 or form P_aS/P_SS.

### Regular Schur elimination

In the covariance-weighted/minimum-action reader the internal pseudo
innovation has quadratic cost

    nu_b^2 / Omega_S,b,
    Omega_S,b=P_SS,b^-+R_S,b >0.                            (BL6)

The gain column is

    k_b=[K_aS,b, k_SS,b]^T
       =[P_aS,b^-,P_SS,b^-]^T/Omega_S,b.                    (BL7)

Eliminating nu_b at BLOCK level is therefore a regular one-dimensional Schur
complement with denominator Omega_S,b, never P_SS alone.  Equivalently, for
any two-row adjoint lambda=[lambda_a,lambda_S],

    min_nu {
       (nu^2/Omega_S)
       +2 nu (lambda_a K_aS+lambda_S k_SS)
    }

has minimizer

    nu_*=-Omega_S(lambda_a K_aS+lambda_S k_SS)              (BL8)

and decrement

    -Omega_S(lambda_a K_aS+lambda_S k_SS)^2
     =-[lambda_a P_aS+lambda_S P_SS]^2/Omega_S.             (BL9)

This is exactly the Joseph/information cancellation but performed BEFORE
separating AW and S.  It remains finite as P_SS->0 and automatically retains
the sign/correlation of P_aS.

### Chronological composition

Let L_b be the two-row reader map at the end of block b.  Pull it backward
through the terminal S update:

    lambda_b^- =
      (I-h_S^T K_b^T) lambda_b^+
      = lambda_b^+ - h_S^T(K_b^T lambda_b^+).               (BL10)

Then pull through F_b:

    lambda_b^root = F_b^T lambda_b^-.                       (BL11)

The source coefficient is

    z_b = G_b^T lambda_b^-.                                 (BL12)

Thus an S service changes only the S component of the backward adjoint by the
scalar K_b^T lambda.  Its covariance-weighted action is charged exactly by
BL9.  No eventwise |K_aS|, no P_aS/P_SS ratio, and no sqrt(N_S) appears.

Across multiple S blocks the complete signed reader is obtained by repeating
BL10--BL12.  This is the desired block Gaussian elimination.

### What remains after the elimination

BL10 shows that S cannot be summarized by a tuner-only scalar kernel: the
elimination coefficient depends on the CURRENT two-row adjoint through
K_b^T lambda.  But this dependence is favorable: BL9 gives an exact negative
square in the same variable.  Therefore the correct bound is a block
completed-square inequality, not TV of an S gain.

For each block, combine:
1. the forward-decayed OU physical coefficient (mass<=1, D1<=2);
2. the endogenous accelerometer source coefficient inside G_b;
3. the terminal S square BL9.

The remaining local inequality has schematic form

    signed physical/source contribution
      - [lambda_a P_aS+lambda_S P_SS]^2/Omega_S
      <= block supply.                                     (BL13)

A large S-to-AW coupling increases the negative square rather than the
variation charge.  This is exactly the cancellation lost by Joseph+Cauchy.

The next quantitative target is therefore to complete the square between the
physical OU/source coefficient z_b and BL9 over one service block, deriving a
source-uniform block supply constant C_b.  Summing C_b over blocks preserves
sign and does not incur sqrt(N_S).


## 33. Block completed square: project source reader onto terminal S innovation

For one S-service block, let the scalar signed source reader after backward
transport be

    r_b = z_b^T y_b.

Let the terminal pseudo innovation be

    nu_b = -S_(b+1)^-
         = n_b^T y_b + nu_root,b,                           (BC1)

where n_b^T=-e_S^T G_b and nu_root,b=-e_S^T F_b x_b^+.
Work in the same covariance/action inner product used by the complete-word
reader.  Decompose the SOURCE part of r_b into its projection on the source
part of nu_b and an orthogonal residual:

    r_b = alpha_b nu_b + r_b^perp,                          (BC2)

with

    alpha_b = <r_b,nu_b>/<nu_b,nu_b>
            = Cov(r_b,nu_b)/Omega_b                         (BC3)

when the complete innovation variance Omega_b is used, and
<r_b^perp,nu_b>=0.  Root pieces are retained separately rather than hidden in
the source projection.

The terminal S Schur elimination contributes the negative quadratic action

    -nu_b^2/Omega_b

in normalized innovation coordinates (equivalently BL9 in adjoint
coordinates).  Therefore the correlated source component completes exactly:

    alpha_b nu_b - nu_b^2/Omega_b
      <= Omega_b alpha_b^2/4
       = Cov(r_b,nu_b)^2/(4 Omega_b).                       (BC4)

No triangle inequality is used.

The unexplained source is only r_b^perp.  Its variance/action is the Schur
residual

    Var(r_b^perp)
      = Var(r_b)-Cov(r_b,nu_b)^2/Omega_b.                   (BC5)

Thus the pseudo measurement splits the block source reader into:
1. an S-correlated component, charged by the completed square BC4;
2. an orthogonal component BC5 that the S update cannot control.

This is the exact same-history cancellation sought in the block argument.

### Covariance form without explicit projection coefficient

Let the joint source covariance of (r_b,nu_b) be

    Sigma_b=[[V_r,C_rS],[C_rS,Omega_b]].

Then the conditional/source residual is

    V_perp = V_r-C_rS^2/Omega_b >=0.                        (BC6)

The block bound can be written

    r_b - nu_b^2/Omega_b
      <= r_b^perp + C_rS^2/(4 Omega_b),                    (BC7)

with r_b^perp handled by its declared physical/bias/sensor source class.
Large correlation C_rS REDUCES V_perp; it is not an independent defect.

### Important limitation

BC4 is an action/quadratic completion.  The nominal-AW theorem is a signed
linear mean bound.  To turn BC5 into a deterministic mean supply C_b, each
source class still needs its own deterministic constraint:
- physical acceleration: bounded velocity/jerk, handled by the already
  forward-decayed OU kernel and summation by parts;
- commissioned sensor residual: amplitude box;
- BA/attitude/field terms: retained deterministic envelopes;
- endogenous accelerometer innovation: must remain in the literal feedback
  identity, not be assigned a stochastic covariance norm.

Therefore the S block completion does NOT by itself bound the whole z_b^T y_b
by a covariance variance.  It removes exactly the component aligned with the
S innovation and leaves a smaller signed reader on the deterministic source
classes.

### Resulting block supply structure

After OU telescoping and S completion, one block has the form

    R_b
      <= C_phys,b + C_BA,b + C_sensor,b + C_impl,b
         + C_Scorr,b
         + R_acc,endog,b,                                   (BC8)

where

    C_Scorr,b = C_rS^2/(4 Omega_b),                         (BC9)

and the physical coefficient retains the forward OU mass/variation bounds.
R_acc,endog,b is the ONLY term that cannot be bounded without using the
accelerometer feedback equation itself.

Summing blocks does not incur sqrt(N_S): every S innovation is eliminated
locally by its own negative square before block supplies are added.

This closes the S-chain structurally.  The remaining critical calculation is
now to eliminate R_acc,endog,b with the accelerometer measurement/update
identity in the same way, but its innovation covariance contains AG/BA.
For the SIGNED mean theorem, use the deterministic identity
r_acc=y_phys-h(x_hat), not an innovation-covariance lower bound.


## 34. Deterministic accelerometer-innovation elimination

The literal accelerometer mean model is

    f_pred = R_wb (a_hat_w-g) + lever + b_a(T),
    r_acc  = f_meas-f_pred.

Write the physical body measurement on the same history as

    f_meas =
      R_true (a_phys-g) + lever_true + b_a,true(T) + eta_a.

Transport the residual to the nominal world frame.  Then exactly

    R_wb^T r_acc
      = a_phys-a_hat_w
        + d_att + d_BA + d_lever + eta_w,                  (AC1)

where d_att=(R_wb^T R_true-I)(a_phys-g), d_BA is the physical-minus-nominal
temperature-dependent accelerometer bias, d_lever is the lever-model defect,
and eta_w=R_wb^T eta_a.  No innovation covariance appears.

Let

    Gamma_k=K_aw,k R_wb,k.

The AW mean correction is therefore

    Delta a_hat_k^acc
      = Gamma_k(a_phys,k-a_hat_k)
        + K_aw,k[d_body,k],                                 (AC2)

or, collecting declared world defects,

    a_hat_k^+
      =(I-Gamma_k)a_hat_k^-
        +Gamma_k a_phys,k
        +d_acc,k.                                           (AC3)

This is the central deterministic cancellation: Gamma multiplies the physical
acceleration and the negative nominal AW with the SAME matrix.  They must not
be bounded separately.

### Combine with OU prediction

For the simple literal ordering prediction -> (optional S) -> accelerometer,
temporarily denote the post-S pre-accelerometer AW by a_tilde_k.  Prediction
from the previous posterior is

    a_pred,k=phi_k a_hat_(k-1)^+.

The S update contributes xi_S,k, so

    a_tilde_k=phi_k a_hat_(k-1)^+ + xi_S,k.                 (AC4)

Substitute AC4 into AC3:

    a_hat_k^+
      =(I-Gamma_k)phi_k a_hat_(k-1)^+
       +Gamma_k a_phys,k
       +(I-Gamma_k)xi_S,k
       +d_acc,k.                                            (AC5)

Thus the homogeneous closed-loop AW multiplier is

    A_k=(I-Gamma_k)phi_k,                                   (AC6)

and the physical input coefficient is Gamma_k.  The same-history S correction
is attenuated by I-Gamma_k if it precedes the accelerometer update.

For arbitrary scheduler ordering, AC5 generalizes by composing the literal
rank-one S maps and accelerometer map in actual order; the key pair
(I-Gamma),Gamma remains exact at every accepted accelerometer correction.

### Error form

Define e_k=a_hat_k-a_phys,k.  From AC3,

    e_k^+
      =(I-Gamma_k)e_k^-
        +d_acc,k,                                           (AC7)

at the measurement instant (same physical sample).  This is stronger than
the earlier AW-loop rearrangement: physical acceleration CANCELS from the
measurement error recursion exactly.  Across prediction,

    e_pred,k
      =phi_k e_(k-1)^+
       -(1-phi_k)a_phys,k-1
       -Delta a_phys,k.                                     (AC8)

Hence the only physical forcing of the error occurs through the OU prediction
identity already handled by the forward-decayed beta_OU kernel.  The
accelerometer update does not introduce an independent physical-acceleration
reader at all; it applies I-Gamma to the existing error.

### Complete deterministic word

Iterating AC7--AC8 with the exact S block maps gives

    e_N
      = M_cl e_0
        + sum_j M_(N:j)
            [ -(1-phi_j)a_phys,j - Delta a_phys,j ]
        + D_S
        + D_att+BA+lever+sensor+impl.                       (AC9)

Here M_cl is the literal closed-loop product of OU, S and (I-Gamma)
accelerometer maps.  Crucially there is NO separate sum Gamma_j a_phys,j.
That term was an artifact of moving Gamma a_hat to the other side before
substitution.

This removes the remaining endogenous accelerometer source term
R_acc,endog structurally.

### What remains quantitatively

The nominal AW mean is

    a_hat = a_phys+e.

On a long window, the physical mean is bounded by velocity endpoints and
jerk as already used.  The error reader AC9 has:
- OU physical kernel with exact variable-phi mass<=1 and D1<=2 BEFORE the
  intervening correction maps;
- homogeneous root term M_cl e_0;
- S corrections, already eliminable by the two-row block Schur completion;
- declared deterministic attitude/BA/lever/sensor/implementation defects.

The remaining mathematical issue is now the effect of the matrices
(I-Gamma_k) on the forward OU kernel.  If they are nonexpansive in the
relevant signed/readout metric, the beta_OU mass/variation bounds survive and
the physical contribution closes.  Euclidean nonexpansiveness is NOT
automatic for an arbitrary Kalman gain.  The correct metric is the
covariance/information metric in which a Kalman correction is contractive.

Therefore the last feedback lemma is:

    the forward-decayed OU physical reader transported through literal
    accelerometer corrections has mass/variation no larger than its
    covariance-weighted complete-word reader bound, without converting back
    to Euclidean gain norms.                                (AC10)

This is now a pure correction-transport lemma; there is no endogenous
innovation source left.


## 35. Full-metric correction transport and final-reader constant

For every literal accepted Kalman correction with A=I-KH, Joseph gives

    P^+ = A P^- A^T + K R K^T >= A P^- A^T.

Hence

    A^T (P^+)^-1 A <= (P^-)^-1,                            (MT1)

so full-state covariance-weighted error is nonexpansive.  This is exactly the
actual-gain prefix inequality already proved in the corrected-word theorem.

The AW-only block I-Gamma does NOT inherit MT1 after projection: accelerometer
corrections exchange energy among attitude, BG/lever, AW and BA coordinates.
Therefore forward OU physical defects must be embedded in the full state,
transported through the full A matrices, and projected to the AW reader only
at the end.

For a unit transverse AW row u and a transported AW column
x=M_(N<-j) E_aw u,

    |u^T E_aw^T x|
      <= sqrt(u^T P_aa,N u) ||x||_(P_N^-1)
      <= sqrt[
           (u^T P_aa,N u)
           (u^T E_aw^T P_j^-1 E_aw u)
         ].                                                 (MT2)

Thus a valid reader constant is the SAME-HISTORY product

    c_read^2 =
      sup_(shipping,j,u)
       (u^T P_aa,N u)
       (u^T E_aw^T P_j^-1 E_aw u).                         (MT3)

The second factor is the inverse conditional AW covariance at epoch j.

### Independent covariance extrema are useless

The recurring nuisance upper comparison permits

    P_aa <=156^2=24336.

The recurring post-prediction lower covariance certificate has a very small
generic conditional-AW floor (order 1e-9 in the fixed proof coordinates).
Combining these independently gives c_read of order 10^6, far too large for
the signed OU Abel bound.  This is another demonstrated relaxation failure:
the two extrema come from incompatible histories/directions.

Therefore c_read must be bounded as the product MT3 on ONE realized
covariance chronology.  The required theorem is a same-word covariance ratio,
not separate upper/lower boxes.

### Riccati-order route

Let R_(N<-j) denote the literal Riccati map from epoch j to N.  For the actual
P_j,

    P_N=R_(N<-j)(P_j).

The desired scalar is

    F(P_j)=
      [u^T R(P_j)_aa u]
      [u^T (P_j)_(a|rest)^-1 u].                           (MT4)

Riccati monotonicity alone does not make F monotone because the two factors
move in opposite directions.  But this is precisely a projective/diameter
quantity: it compares a terminal marginal to an initial conditional
precision along the same Riccati trajectory.

The existing Theorem-D interval

    Pi <= P_N <= P_diff

and its root-to-terminal information ratio can potentially bound MT4 without
a covariance ceiling.  The next calculation should express MT4 through the
joint Gaussian (x_j,y_word,x_N) Schur complements and reduce it to the
word diameter kappa_W, rather than to independent covariance extrema.

If one can prove

    c_read^2 <= kappa_W                                    (MT5)

(or a modest fixed multiple), then

    c_read <= sqrt(kappa_W),

and the already established word contraction/diameter certificate directly
supplies the final reader conversion.  This would connect the signed nominal
AW proof to Theorem D at exactly one final point, without reintroducing gain
norms or G0.


## 36. The proposed c_read^2 <= kappa_W inequality is false

The same-history product

    (terminal marginal variance)*(initial conditional precision)

is NOT controlled by the Riccati diameter kappa_W in general.

Scalar counterexample: let

    y=x_0+v,   Var(v)=R=1,
    x_N=x_0+w, Var(w)=Q=1.

For a known root, terminal covariance is

    Pi=Q=1.

For a diffuse root, the measurement leaves variance R=1 before prediction,
so

    P_diff=Q+R=2,

and Theorem D gives

    kappa_W=P_diff/Pi=2.                                   (CR1)

For an actual root prior Var(x_0)=p>0, the posterior after y is p/(p+1), so

    P_N=1+p/(p+1).

The proposed reader product is

    P_N * p^-1
      = 1/p + 1/(p+1),                                     (CR2)

which tends to infinity as p->0 while kappa_W remains 2.  Hence no universal
inequality

    c_read^2 <= kappa_W

(or any fixed multiple independent of the root covariance) can hold.

The failure is conceptual: kappa_W measures the DIAMETER of terminal Riccati
covariances as the root prior varies from known to diffuse.  The reader
product multiplies terminal variance by INITIAL precision; an arbitrarily
well-known root has arbitrarily large initial precision even though the
terminal diameter remains finite.

### Correct normalization

The transported deterministic column generated at epoch j must be normalized
by the covariance/source channel that actually CREATES that column, not by
the total initial conditional precision.  For OU physical forcing, the
coefficient beta_j is deterministic and its sharp bounded-velocity estimate
already supplies the source normalization.  Using P_j^-1 on E_aw beta_j
penalizes a perfectly known AW root even though the physical forcing occurs
AFTER that root and is unrelated to its uncertainty.

Therefore the final-reader conversion should start at the injection epoch
AFTER the OU physical defect is added.  If the defect is d_j=E_aw q_j, let
P_j^def be the covariance metric immediately after the corresponding
prediction/process channel.  Since prediction adds Q_j,

    P_j^def >= Q_j                                         (CR3)

on the process-supported subspace, and subsequent corrections are
nonexpansive.  A terminal AW readout obeys

    |u^T E_aw^T M d_j|
      <= sqrt(u^T P_aa,N u)
         sqrt(d_j^T (P_j^def)^-1 d_j)
      <= sqrt(u^T P_aa,N u)
         sqrt(d_j^T Q_j^dag d_j).                           (CR4)

This is source/action normalization, not root-precision normalization.

CR4 alone may still be quantitatively loose for deterministic physical
acceleration because Q_j can be tiny.  The sharp beta_OU/Abel bound and the
metric correction transport therefore need a mixed argument: keep the
physical sequence in its signed primitive norm through prediction, and use
Kalman metric contraction only for correction maps.  There is no scalar
c_read depending only on kappa_W that automatically converts between these
two norms.

Conclusion: Theorem D does not close the final reader bridge.  The remaining
lemma is a two-norm transport problem: show that interposed Kalman correction
maps do not amplify the OU signed-primitive operator from bounded physical
velocity to final AW readout.  This requires exploiting the special
accelerometer correction structure, not only generic Riccati diameter.


## 37. Accelerometer passivity: false for isolated AW, exact at the full measurement port

The deterministic AW error update is

    e_aw^+=(I-Gamma)e_aw^- + cross/defect terms,
    Gamma=K_aw R_wb.

From the literal gain construction,

    K_aw =
      [ P_a,theta J_att^T
        +P_aa R_wb^T
        +P_a,ba
        +P_a,bg J_bg^T ] S_acc^-1,                          (PA1)

with the BA/BG terms present according to the shipping mode.  Therefore

    Gamma=K_aw R_wb                                         (PA2)

contains cross-covariance terms of unrestricted sign.  PSD of P and SPD of
S_acc do NOT imply sym(Gamma)>=0, nor ||I-Gamma||_2<=1.
Consequently there is no generic Euclidean positive-real/passivity theorem
for the isolated AW correction.  A zero-mean physical input could be
rectified by an arbitrary algebraic AW gain if the coupled attitude/BA
coordinates are discarded.

The full accelerometer correction DOES have an exact passive-port identity.
Let the full linearized error be e, measurement row H, physical/model defect
d, and innovation

    r=-H e + d                                               (PA3)

(up to the fixed sign convention).  The update is

    e^+=A e + K d,   A=I-KH.                               (PA4)

Joseph gives

    P^+=A P^- A^T+K R K^T.                                 (PA5)

Equivalently in the complete-word action/dual formulation, the measurement
port contributes the nonnegative square associated with S_acc=H P^- H^T+R.
Completing that square keeps the SAME combination H e that contains attitude,
AW, BA and lever/BG effects.  Thus accelerometer feedback is passive in the
full measurement port, not in the AW coordinate alone.

### Consequence for the signed physical-input theorem

The remaining two-norm lemma cannot be

    isolated AW correction preserves OU beta variation.

That statement is false without additional cross-covariance restrictions.

Instead combine the forward OU physical forcing with the FULL accelerometer
port before projection, exactly as section 32 combined AW and S.  Over an
accelerometer-service block retain two objects:
1. the transverse AW readout adjoint;
2. the accelerometer predicted-specific-force row H e.

Eliminate the internal accelerometer innovation by a block Schur/completed
square with denominator S_acc, but substitute the deterministic physical
measurement identity first.  The physical acceleration appears in both the
OU forcing and H e with fixed opposite signs.  Their signed combination is
the candidate positive-real supply rate.

This block-level port formulation has the required properties:
- no sym(Gamma) assumption;
- no gain norm;
- no AG covariance ceiling;
- all attitude/BA cross covariance retained;
- same-history physical acceleration appears once;
- the negative measurement square grows when cross coupling is large.

The next analytical object is therefore a FULL accelerometer-port block
matrix, analogous to the successful two-row AW/S block:
    [final AW reader ; accelerometer force-error port].
Its Schur complement in the innovation variable is regular because
S_acc>=R_acc>0.  The question becomes whether the resulting deterministic
physical supply has nonpositive DC gain (up to declared attitude/BA/lever
defects).  This is the correct positive-real formulation.


## 38. Full AW-reader / accelerometer-port Schur elimination

Consider one accepted accelerometer correction.  Let e^- be the full
linearized state error immediately before correction, H the literal full
accelerometer Jacobian, d the deterministic physical/model/sensor defect after
putting the true physical acceleration on the same history, and

    r = -H e^- + d.                                        (AP1)

The correction is

    e^+ = e^- + K r,                                       (AP-two)

with

    S_acc = H P^- H^T + R_acc >0,
    K=P^- H^T S_acc^-1.                                    (AP3)

Let lambda^+ be an arbitrary backward reader at the post-correction state.
Then

    lambda^T e^+
      = lambda^T e^- + (K^T lambda)^T r.                   (AP4)

Put

    q=K^T lambda=S_acc^-1 H P^- lambda.                     (AP5)

The innovation is therefore one internal 3-vector entering the reader through
q^T r.

### Regular 3x3 Schur completion

In the complete covariance-weighted action, the accelerometer innovation has
quadratic cost r^T S_acc^-1 r (up to the fixed action convention).  Completing
the square,

    q^T r - r^T S_acc^-1 r
      <= (1/4) q^T S_acc q                                (AP6)

for the convention with unit quadratic coefficient.  Equivalently,

    (1/4) q^T S_acc q
      =(1/4) lambda^T P^- H^T S_acc^-1 H P^- lambda.        (AP7)

Thus the full accelerometer coupling is charged by an exact PSD Schur square.
No gain norm, AW-only Gamma sign, AG covariance ceiling, or cross-covariance
box is required.

The backward reader itself pulls through as

    lambda^-=(I-H^T K^T)lambda^+,                           (AP8)

the full-state analogue of the S-block formula.

### Deterministic physical substitution

For the literal shipping model,

    r =
      R_wb(a_phys-a_hat_w)
      + d_att_body+d_BA+d_lever+eta_a.                     (AP9)

Hence

    q^T r =
      (R_wb^T q)^T a_phys
      -(R_wb^T q)^T a_hat_w
      +q^T d_decl.                                         (AP10)

The physical acceleration and nominal AW have EXACTLY the same coefficient
with opposite signs.  Define

    gamma_lambda = R_wb^T q
                 = R_wb^T S_acc^-1 H P^- lambda.            (AP11)

Then the accelerometer port contribution is

    gamma_lambda^T(a_phys-a_hat_w)+q^T d_decl.              (AP12)

This is the full-state version of the cancellation previously seen in the AW
row, but now all attitude/BA/BG cross covariance remains inside q and the
negative Schur square AP7.

### Combine with OU prediction before bounding

At the preceding OU prediction, the AW error receives

    d_OU=-(1-phi)a_phys-Delta a_phys.                       (AP13)

Transport its AW reader to the accelerometer port.  The combined signed
physical supply over prediction+correction is therefore of the form

    c_OU^T[-(1-phi)a_phys-Delta a_phys]
      +gamma_lambda^T(a_phys-a_hat_w)
      - measurement_square
      +declared defects.                                   (AP14)

Do NOT bound gamma_lambda^T a_phys and
-gamma_lambda^T a_hat_w separately.  Substitute
a_hat_w=a_phys+e_aw at the same correction epoch:

    gamma_lambda^T(a_phys-a_hat_w)
      =-gamma_lambda^T e_aw.                               (AP15)

Thus the accelerometer port contributes NO independent physical-acceleration
forcing.  Physical acceleration remains only in AP13, the OU prediction
forcing.  This recovers the deterministic error recursion AC7 at the full
reader/action level while retaining the exact negative measurement square.

### Consequence

The accelerometer Schur block therefore closes structurally:
- endogenous innovation eliminated;
- physical acceleration cancels from the correction port;
- attitude/BA/lever/sensor enter only through declared d_decl;
- arbitrary cross covariance strengthens/changes the PSD Schur square but
  cannot create a free physical source.

After every accepted accelerometer correction, the backward reader is AP8.
The only physical forcing over the whole word is the sequence of OU
prediction defects AP13.

The unresolved two-norm issue is correspondingly narrower: transport of
those OU defects through AP8 may rotate the final AW reader into
attitude/BA coordinates, but it cannot create another a_phys source.  To
bound the signed physical primitive, perform Abel summation on the FULL
backward reader's AW component at the prediction epochs, not on an isolated
AW gain.

Let

    b_j = E_aw^T lambda_j^- (1-phi_j)                       (AP16)

be the vector coefficient of a_phys at prediction j after all future full
corrections have been pulled backward.  The exact final physical reader is

    sum_j b_j^T a_phys,j + corresponding Delta-a terms.     (AP17)

The final remaining quantitative lemma is now simply a variation bound on the
AW COMPONENT of the full backward reader:

    |b_1|+|b_N|+sum|b_(j+1)-b_j| <= C_full-reader.          (AP18)

Unlike the abandoned beta-TV problem, AP18 contains no independent
accelerometer innovation or S gain; both have already been Schur-eliminated.
It is the literal complete-word backward reader variation and is the correct
object for bounded-velocity Abel summation.


## 39. Exact 0.1-s block reader and two-norm coarse graining

Partition a regular word at the applied-tuner/covariance-sync activation
boundaries t_B, with block length L_B about .1 s.  Freeze the realized
coefficients only for the auxiliary linear comparison.  Collapse all literal
homogeneous operations in block B into

    e_(B+1)=M_B e_B + sum_(j in B) M_(B+1<-j) d_j.          (CB1)

All accepted accelerometer/S/magnetic corrections are inside M_B.  Their
measurement/source defects retain the existing square-summed factor
representation; do not convert them to eventwise l1 norms.

Let Lambda_(B+1) be the full backward reader at the block end.  The block-root
reader is

    Lambda_B=M_B^T Lambda_(B+1).                            (CB2)

For a physical OU prediction defect at sample j,

    d_j^phys=E_aw[-(1-phi_j)a_j-Delta a_j].

Its exact signed coefficient in the block reader is

    c_(B,j)^T
      =-Lambda_(B+1)^T M_(B+1<-j) E_aw.                    (CB3)

Thus the block physical contribution is

    R_phys,B =
      sum_(j in B) c_(B,j)^T
        [(1-phi_j)a_j+Delta a_j].                           (CB4)

Because the future transport M_(B+1<-j) includes internal correction maps and
LIN prediction mixing, c_(B,j) is NOT generally a common boundary AW reader
times the scalar pure-OU kernel.  Claiming
R_phys,B=(1-Phi_B) r_B a would discard the very correction rotations being
proved about.

### Internal summation by parts

Write a_j=(v_(j+1)-v_j)/h_j plus the bounded jerk remainder.  Apply discrete
summation by parts only inside B.  This gives

    R_phys,B =
      B_B^R v_(B+1)-B_B^L v_B
      + R_var,B + R_jerk,B,                                (CB5)

where B_B^{L,R} are exact block endpoint reader coefficients and R_var,B
contains only differences of c_(B,j)/h_j INSIDE the block.

Do not bound R_var,B by raw total variation.  Split each difference into:
1. OU prediction transport, whose scalar forward kernel has mass<=1 and
   D1<=2;
2. correction-induced reader jumps.

For an accelerometer correction the AW-reader jump is
-R_wb^T q_acc and its squared norm is bounded by the exact corrected loss
q_acc^T S_acc q_acc / lambda_min(R_acc).  S, magnetic and covariance sync
have no direct AW-reader jump; their indirect effect enters through later
prediction mixing and is retained in the full block metric.

Therefore the internal correction part obeys a block l2 estimate

    ||R_var,B^corr||
      <= C_B sqrt(Loss_B) sqrt(PrimitiveEnergy_B),          (CB6)

with C_B depending only on the finite block horizon/scales, not on the number
of corrections individually.  The exact sharp C_B still needs derivation;
using eventwise Cauchy would give sqrt(n_B) and is not promoted.

### Global block composition

Sum CB5 over blocks.  Adjacent velocity-boundary terms combine as

    sum_B [B_B^R v_(B+1)-B_B^L v_B]
      = endpoint terms
        +sum_internal (B_(B-1)^R-B_B^L) v_B.               (CB7)

Thus bounded physical velocity is charged only by BLOCK-TO-BLOCK reader
variation, not sample-level variation.

The residual correction terms should be accumulated with the existing
complete-word square-summed input identity:

    sqrt(sum_B Loss_B)                                     (CB8)

rather than sum_B sqrt(Loss_B).  This avoids a sqrt(number of blocks) loss.

### Remaining quantitative lemma

The useful block theorem is now:

    D_block =
      |B_0^L|+|B_last^R|
      +sum_B |B_(B-1)^R-B_B^L|
      <= C_boundary,                                       (CB9)

and

    sum_B ||R_var,B^corr||
      <= C_corr sqrt(total corrected loss)
                   * physical-primitive budget,             (CB10)

with C_boundary,C_corr source-uniform and modest.

CB9 is a variation bound on only ~170 boundary readers.  CB10 keeps all
~3400 internal corrections in their natural l2 action.  This is the exact
two-norm coarse graining needed for the signed physical theorem.

The next calculation is to express B_B^L,B_B^R in terms of Lambda_B and the
block OU tail products, then test whether CB9 telescopes across the covariance
sync boundary (sync has identity mean map) so that only tuner changes and
block corrected-loss remainders contribute.


## 40. Exact block endpoint coefficients: physical OU forcing is a coboundary

The physical prediction forcing must be kept as the pair

    (1-phi_j)a_j + Delta a_j
      = a_(j+1)-phi_j a_j.                                 (BE1)

Let r_j=E_aw^T lambda_j denote the full backward reader's AW component at the
appropriate prediction boundaries.  Ignoring corrections for one moment, the
physical contribution of prediction j is

    -r_(j+1)^T [a_(j+1)-phi_j a_j].                        (BE2)

Summing over a prediction-only block j=m,...,n gives

    R_phys,B =
       phi_m r_(m+1)^T a_m - r_(n+1)^T a_(n+1)
       +sum_(j=m+1..n)
          [phi_j r_(j+1)-r_j]^T a_j,                        (BE3)

where r_j is the AW reader immediately after pulling through prediction
j-1.  Thus the internal coefficient is the FAILURE of the pure prediction
adjoint relation.

For the full LIN prediction, BR3 gives

    r_j =
      phi_j r_(j+1)
      +phi_va,j lambda_v,j+1
      +phi_pa,j lambda_p,j+1
      +phi_Sa,j lambda_S,j+1.                              (BE4)

Hence

    phi_j r_(j+1)-r_j
      =-[phi_va lambda_v
         +phi_pa lambda_p
         +phi_Sa lambda_S]_(j+1).                           (BE5)

So even before measurement corrections, the only internal physical
coefficient is the integrated-LIN reader mixing; there is no raw OU
variation term at all.

### Include corrections

At a correction between predictions, the full reader changes by

    lambda^- - lambda^+ = -H^T q.                           (BE6)

For S and magnetic updates H_aw=0, so they do not directly alter r.  For an
accelerometer update,

    Delta r=-R_wb^T q_acc.                                 (BE7)

Therefore the exact internal coefficient in BE3 is

    phi_j r_(j+1)-r_j
      = -c_L,j^T lambda_L,j+1
        + correction_jump_terms,                            (BE8)

with c_L=[phi_va,phi_pa,phi_Sa].  The correction jump terms are nonzero only
for accelerometer ports in the AW component and are square-summed by the
corrected-loss identity.

### Explicit block endpoint coefficients

Comparing BE3 with

    R_phys,B =
      B_B^R a_(n+1)-B_B^L a_m + R_internal,B,

the exact endpoint coefficients are

    B_B^L = -phi_m r_(m+1),                                 (BE9)
    B_B^R = -r_(n+1),                                       (BE10)

up to the fixed sign convention of R_phys.  Equivalently, if the block
boundary is chosen immediately before the first prediction, pull r_(m+1)
through that prediction and write

    B_B^L = -[r_m-c_L,m^T lambda_L,m+1].                    (BE11)

Thus B_L and B_R are boundary AW readers plus a single prediction-mixing
term; they do NOT contain a sum of internal OU coefficients.

### Adjacent block cancellation

Choose blocks at identity-mean covariance-sync boundaries.  The final reader
of block B and initial reader of B+1 are the SAME full reader at the common
boundary because sync has identity mean map.  Therefore the leading AW-reader
parts of

    B_B^R-B_(B+1)^L

cancel exactly.  The mismatch is only:
1. the first-prediction LIN mixing term of the new block;
2. any accelerometer correction whose literal event ordering lies exactly at
   the boundary.

There is NO tuner-jump term by itself: phi of the new block appears only in
the local first-prediction identity BE11, and the physical forcing pair BE1
remains exact for arbitrary phi.

Hence the block-boundary variation is controlled by the integrated-LIN
reader at ONE prediction per block plus boundary correction ports, not by
variation of tau.

### Two-norm global bound

The integrated-LIN boundary charge is

    g_B =
      phi_va lambda_v
      +phi_pa lambda_p
      +phi_Sa lambda_S.                                    (BE12)

The accelerometer boundary jump is R_wb^T q_acc.  Both have natural
covariance-weighted l2 controls:
- g_B through the LIN reader metric/process covariance;
- q_acc through q_acc^T S_acc q_acc.

Therefore

    D_block
      <= endpoint_reader_charge
         + sum_B ||g_B||
         + sum_boundary_acc ||q_acc||.                      (BE13)

Do NOT apply Cauchy separately to the two sums.  Stack all g_B and q_acc as
columns of the complete-word factor/reader operator and use its single
operator-norm <=1 identity.  The physical boundary velocities are then the
input coefficients.  The desired inequality has the form

    |sum_B g_B^T v_B + acc-boundary terms|
      <= C_2 sqrt(sum_B ||v_B||^2_weighted),                (BE14)

which is still l2 in block velocities.  To exploit only |v_B|<=Vmax without
sqrt(N_B), one needs additional sign/variation structure of g_B.  Thus exact
block cancellation removes all OU/tuner variation, but the LIN-reader
boundary sequence remains the final l1-vs-l2 obstacle.

This is much narrower than CB9: C_boundary is not yet modest from existing
loss identities alone.  The remaining sequence is specifically the
integrated-LIN mixing g_B at one prediction per .1-s block.


## 41. Second Abel step is admissible: use the declared physical displacement bound

The authoritative MARINE MOTION constants include

    P_max=8.1 m,
    V_max=5.5 m/s,
    A_max=8.8 m/s^2,
    J_max=100 m/s^3.

Thus a second summation-by-parts step from physical velocity to physical
displacement does NOT strengthen the assumptions.

After the first exact block coboundary reduction, the unresolved physical
boundary term has the form

    R_g = sum_(B=1..M) g_B^T v_B,                           (A2-1)

where, at one prediction per ~.1-s block,

    g_B =
      phi_va,B lambda_v,B
      +phi_pa,B lambda_p,B
      +phi_Sa,B lambda_S,B.                                (A2-2)

Let L_B=t_(B+1)-t_B and p_B be the physical displacement primitive at the
block boundary.  Taylor/integral remainder gives

    v_B = (p_(B+1)-p_B)/L_B + eps_v,B,                     (A2-3)

with

    ||eps_v,B|| <= A_max L_B/2                              (A2-4)

when v_B is the boundary velocity at the chosen end; the corresponding
one-sided convention changes only the sign of the remainder.

Define

    h_B = g_B/L_B.                                          (A2-5)

Then

    R_g =
      sum_B h_B^T(p_(B+1)-p_B)
      +sum_B g_B^T eps_v,B.                                (A2-6)

Discrete Abel gives exactly

    sum_B h_B^T(p_(B+1)-p_B)
      = h_M^T p_(M+1)-h_1^T p_1
        -sum_(B=1..M-1)(h_(B+1)-h_B)^T p_(B+1).             (A2-7)

Therefore

    |R_g|
      <= P_max[
           ||h_1||+||h_M||
           +sum||Delta h_B||
         ]
         +(A_max/2) sum_B L_B ||g_B||.                     (A2-8)

This replaces the previous V_max * sum||g_B|| charge by:
- P_max times first variation of the much smoother h_B=g_B/L_B;
- an acceleration remainder weighted by L_B.

### Divided-difference structure of g_B

The integrated-OU coefficients satisfy

    phi_va = tau(1-phi),
    phi_pa = tau^2(x+exp(-x)-1),
    phi_Sa = tau^3(.5x^2-x-exp(-x)+1),  x=h/tau.           (A2-9)

For small h these are respectively

    h+O(h^2/tau),
    h^2/2+O(h^3/tau),
    h^3/6+O(h^4/tau).                                      (A2-10)

Hence at one prediction,

    g_B/L_B
      ~ (h/L_B) lambda_v
        +(h^2/(2L_B)) lambda_p
        +(h^3/(6L_B)) lambda_S.                            (A2-11)

Since h~.005 and L_B~.1, the p and S reader contributions receive additional
small factors ~.00125 and ~2.1e-5 relative to their raw components.  The
leading v-reader factor is ~.05.

This is the quantitative reason the second Abel step is promising.

### What remains to prove

The exact first variation is

    Delta h_B =
      Delta[(phi_va/L_B)lambda_v]
      +Delta[(phi_pa/L_B)lambda_p]
      +Delta[(phi_Sa/L_B)lambda_S].                         (A2-12)

Reader jumps inside the block have already been Schur-eliminated/square-summed.
At sync boundaries the mean reader is continuous.  Therefore Delta h_B is
generated by:
1. one prediction/tuner coefficient change;
2. accumulated full-reader evolution over the preceding block, controlled in
   the complete-word metric.

A useful theorem now needs a source-uniform bound on the SECOND divided
difference/variation of these three scaled LIN reader components.  If direct
l1 variation is still too loose, the p and S terms admit third/fourth Abel
steps because physical p is already bounded but no additional physical
primitive beyond displacement is declared; therefore only one more Abel step
is legally available on the physical side.  Any further smoothing must come
from the reader dynamics themselves, not a new physical assumption.

The target is

    D2_reader =
      ||h_1||+||h_M||+sum||Delta h_B||,                    (A2-13)

with

    |R_g| <= 8.1 D2_reader
             +4.4 sum_B L_B||g_B||.                        (A2-14)

This is the exact second-Abel physical bound under the existing MARINE MOTION
contract.


## 42. Quantitative second-Abel coefficient geometry and reader-normalization gap

For

    h_B=(1/L_B)[phi_va lambda_v+phi_pa lambda_p+phi_Sa lambda_S],

the integrated-OU coefficients have exact integral representations

    phi_va = int_0^h exp(-s/tau) ds <= h,
    phi_pa = int_0^h (h-s) exp(-s/tau) ds <= h^2/2,
    phi_Sa = int_0^h (h-s)^2/2 exp(-s/tau) ds <= h^3/6.    (Q2-1)

These bounds are uniform for arbitrary tau>0 and avoid small-x expansion
remainders.

For h<=.006 and nominal activation-cell L=.1,

    phi_va/L <= .06,
    phi_pa/L <= 1.8e-4,
    phi_Sa/L <= 3.6e-7.                                    (Q2-2)

If actual cell lengths vary, replace .1 by the proved minimum activation-cell
length before using these numbers.  The earlier rough p/S factors in section
41 were too large because they missed powers of h; Q2-2 is the correct
integrated-chain scaling.

Thus, in Euclidean row norm,

    ||h_B||
      <= .06 ||lambda_v||
         +1.8e-4 ||lambda_p||
         +3.6e-7 ||lambda_S||.                             (Q2-3)

The acceleration remainder satisfies

    (A_max/2) L_B ||g_B||
      <=4.4[
          h ||lambda_v||
          +(h^2/2)||lambda_p||
          +(h^3/6)||lambda_S||
        ] L_B,                                              (Q2-4)

so per block its raw coefficient scales are at most
.0264, 7.92e-5 and 1.584e-7 times the respective reader norms for
L_B=.1,h=.006.

### What the existing LIN action certificate does and does not provide

The 16-s LIN matrix certificate proves a covariance lower comparison

    P_LIN >= D A^-1 D /2,

equivalently a precision ceiling on LIN endpoint errors.  It does NOT by
itself bound an arbitrary backward adjoint lambda_L: reader scale is set by
the terminal readout normalization and the complete-word joint reader.

Therefore substituting the LIN precision matrix directly as a bound on
||lambda_v||,||lambda_p||,||lambda_S|| would be invalid.

For the present reader the terminal normalization IS fixed (unit transverse
AW nominal-mean row).  The required quantity is the actual complete-word
minimum-action reader restricted to the LIN boundary coordinates.  If its
action is J_reader, then dual Cauchy gives

    |c^T lambda_L|
      <= sqrt(c^T P_LIN c) sqrt(lambda_L^T P_LIN^-1 lambda_L), (Q2-5)

and the second factor is part of the normalized reader action.  A useful
numeric D2 bound therefore requires a source-uniform upper bound on this
NORMALIZED reader action at the ~.1-s boundaries, not another covariance
lower bound.

### Two-norm formulation

Let a_B be the three-vector of scaled coefficients

    a_B=[phi_va/L_B,phi_pa/L_B,phi_Sa/L_B].

Let W_B be the exact 3x3 covariance/action metric induced on the LIN reader at
the boundary by the complete-word minimum-action construction.  Then

    |h_B|^2 <= (a_B W_B^-1 a_B^T)(lambda_L^T W_B lambda_L). (Q2-6)

The first factor is pure coefficient geometry and is strongly suppressed by
Q2-2.  The second is square-summed reader action.

The same representation applies to Delta h_B using the difference of two
coefficient rows and the chronological reader transport.  If the stacked
operator of all boundary rows has norm C_stack in the complete-word action
metric, then

    D2_reader <= sqrt(M+1) C_stack                           (Q2-7)

by generic Cauchy; this still costs sqrt(170) and is not enough a priori.
To avoid it, one must bound the l1 operator norm of the STACKED divided-
difference rows directly.  This is now a finite deterministic reader matrix
problem, not a covariance theorem.

### Remaining exact calculation

Construct the normalized complete-word reader L_min already used in the
repository, sample its LIN boundary rows at the ~.1-s activation boundaries,
and form the deterministic divided-difference operator

    D2 L_LIN = [h_1; h_2-h_1; ...; h_M-h_(M-1); h_M].       (Q2-8)

The required constant is the induced action-to-l1 norm

    C_D2 = sup_(||z||_action<=1) ||D2 L_LIN z||_(2,1).      (Q2-9)

A generic spectral bound reintroduces sqrt(M); a useful proof needs the
special banded/Volterra structure of L_min.  The coefficient scales Q2-2 make
this plausible, but the existing LIN covariance certificate alone does not
supply C_D2.

Thus no honest explicit D2_reader number follows yet from the current
certificate.  What HAS closed analytically is the coefficient geometry; the
remaining obstruction is precisely the normalized complete-word reader's
divided-difference l1 norm.


## 44. Dual divided-difference factorization: local slab term plus transport commutator

Let Lambda_B be the full residual-functional adjoint of the normalized
terminal-AW minimum-action reader at consecutive actual AW-sync boundaries.
Collapse the literal homogeneous mean operations in slab B into M_B and stack
the observation rows applied inside that slab into O_B, with the corresponding
minimum-action reader weights ell_B.  Backward chronology gives exactly

    Lambda_B = M_B^T Lambda_(B+1) - O_B^T ell_B.            (VD1)

(Here O_B^T ell_B denotes the sum of the individually transported correction
rows in their literal order; it is notation for the exact slab observation
reader, not a commuted measurement stack.)

Define the scaled LIN mixing row A_B by

    h_B = A_B Lambda_B,

where A_B is zero outside LIN (v,p,S) and on those coordinates contains the
literal first-prediction coefficients divided by slab duration:

    A_B|LIN = [phi_va/L_B, phi_pa/L_B, phi_Sa/L_B].         (VD2)

Then an adjacent divided difference is

    h_(B+1)-h_B
      = [A_(B+1)-A_B M_B^T] Lambda_(B+1)
        + A_B O_B^T ell_B.                                 (VD3)

This is the exact common-future-tail cancellation formula.

The second term is LOCAL: it uses only observation-reader weights in slab B.
The first term is the only surviving future-tail dependence and is multiplied
by the transport commutator

    C_B := A_(B+1)-A_B M_B^T.                              (VD4)

Thus the hoped-for statement "all common future columns cancel" is too
strong.  The correct statement is that common future tails survive only
through C_B.

### Dual operator

Let D2 h=[h_1,h_2-h_1,...,h_M-h_(M-1),h_M].  For dual block
vectors y_i, substitute VD3 and interchange slab/source sums.  The dual
functional splits exactly into

    (D2 L)^T y = Local(y) + Tail(y),                        (VD5)

where Local is block-banded in the slab observation/source columns and Tail
is a Volterra sum of C_B^T y_B pulled through future boundary adjoints.

No sqrt(M) is intrinsic to Local: disjoint unit-action source slabs can be
combined by the complete-word square-sum.  The only possible long-horizon
loss is Tail, controlled by the sequence C_B rather than by h_B itself.

### Size of the commutator

A_B already has the uniform coefficient geometry

    |A_v|<=.06,
    |A_p|<=1.8e-4,
    |A_S|<=3.6e-7

for h<=.006,L=.1.  Moreover M_B is the literal near-identity 0.1-s
closed-loop transition.  Expanding VD4,

    C_B =
      (A_(B+1)-A_B)
      - A_B(M_B^T-I).                                      (VD6)

The first term is tuner/cadence coefficient variation.  The second is a
small row A_B multiplying the full slab state transition defect.  This is a
much smaller object than the full reader variation, but it is not zero and
must be enclosed.

Importantly, using ||M_B-I|| as an arbitrary full-state Euclidean norm would
again be disastrous because corrections can rotate AG/BA coordinates.  Only
the columns seen by A_B matter.  Therefore compute/enclose the THREE pulled
rows

    A_B M_B^T                                               (VD7)

directly from the literal slab chronology.  This retains all Kalman
cancellations.

### Source-slab norm reduction

If the local term has per-slab action operator L_B^loc and the commutator
tail has operator C_B Phi_(future), the exact mixed-norm target can be bounded
without event-count Cauchy provided one proves

    sup_s sum_(B<=s) || C_B Phi_(s<-B) ||_* <= C_tail       (VD8)

and

    sup_s ||L_s^loc||_* <= C_local.                         (VD9)

Then

    C_D2 <= C_endpoint + C_local + C_tail.                 (VD10)

These are Volterra row-sum bounds; there is no sqrt(170).  The future
transport in VD8 is the SAME full correction/prediction transport already
present in the normalized reader, not a product of gain norms.

### Status

The infrastructure now exports every ingredient needed to evaluate VD3 on a
frozen shipping word: exact sync boundaries, Lambda_B, slab operation
intervals and durations.  The next source-uniform theorem is therefore not a
generic mixed-norm estimate.  It is an enclosure of the local rows
A_B O_B^T ell_B and the commutator rows C_B=A_(B+1)-A_B M_B^T over the
coupled shipping slab class.

This is narrower than the previous D2 obstruction and identifies precisely
what must be small for the Volterra cancellation to remove the sqrt(170)
factor.


## 45. Local slab leverage bound and commutator decomposition

For slab B write the exact backward recurrence

    Lambda_B=M_B^T Lambda_(B+1)-O_B^T ell_B,

and h_B=A_B Lambda_B.  The local divided-difference contribution is

    l_B^loc=A_B O_B^T ell_B.                               (LC1)

Let Sigma_B be the covariance/action Gram of the slab observation-source rows
in the SAME augmented design after marginalizing the slab's process/source
columns consistently.  Then weighted Cauchy gives

    |l_B^loc|^2
      <= chi_B^2 (ell_B Sigma_B ell_B^T),                  (LC2)

where the local leverage is

    chi_B^2 =
      A_B O_B^T Sigma_B^-1 O_B A_B^T.                     (LC3)

This is the correct source-uniform local constant.  It applies A_B before
Cauchy, retaining the coefficient ceilings
(.06,1.8e-4,3.6e-7), and introduces no count of corrections in the slab.

The complete reader action is the square-sum of the slab actions plus process
and root residual action, so the local terms can be stacked in the action
metric.  A useful theorem needs only

    sup_B chi_B <= chi_0.                                  (LC4)

Because Sigma_B includes the actual effective measurement covariance, a
lower noise floor gives a finite chi_0.  Computing chi_0 by replacing O_B
with an arbitrary row norm would be too coarse; use the literal slab rows.

### Commutator

The nonlocal term is

    C_B Lambda_(B+1),
    C_B=A_(B+1)-A_B M_B^T.                                 (LC5)

Factor the slab homogeneous map chronologically into predictions/resets and
Kalman corrections.  Pull A_B forward through this factorization.  Across a
correction A_k=I-K_k H_k,

    a A_k^T = a - (a H_k^T) K_k^T.                         (LC6)

Thus the correction part of the commutator is a sum of measurement-port rows,
not a free Euclidean matrix defect.  For S and magnetic H_aw=0 but A_B lives
on v,p,S, so S can couple directly through its S column; accelerometer can
couple through the AW/AG/BA columns after prediction mixing.  Each term has
the form

    scalar/vector port coefficient * K_k^T Lambda_future,   (LC7)

and the latter is exactly the q_k variable charged by the corrected-loss
square q_k^T S_k q_k.

Across a prediction F_k, the row evolves deterministically:

    a F_k^T.                                                (LC8)

For the pure LIN integrated chain this can be computed analytically from
phi_va,phi_pa,phi_Sa and the block duration.  Therefore split

    C_B = C_B^pred + C_B^port,                              (LC9)

where C_B^pred is tuner/integration geometry and C_B^port is retained as
corrected-loss ports.

This is essential: there is no source-uniform useful raw Euclidean bound on
C_B from tuner ranges alone because M_B contains covariance-dependent gains.
Trying to bound ||C_B|| directly repeats the failed arbitrary-gain
relaxation.

### Volterra bound in the correct norm

The dual tail is controlled by

    Tail = sum_B C_B^pred Lambda_(B+1)
           + sum_ports alpha_k q_k.                         (LC10)

The port sum is square-summed globally by the existing corrected loss.  Only
C_B^pred requires an l1/Volterra row-sum bound.  It is covariance independent.

Hence the source-uniform constants should be

    chi_0 = sup_B sqrt(A_B O_B^T Sigma_B^-1 O_B A_B^T),    (LC11)

    C_pred = sup_s sum_(B<=s)
       || C_B^pred Phi_pred(s<-B) ||,                       (LC12)

with all correction ports removed from Phi_pred and charged separately in
their action metric.

This is narrower than bounding the full C_B.  The remaining deterministic
matrix calculation is the pure prediction/reset integrated-chain
commutator C_B^pred; the stochastic/Kalman part is already covered by exact
loss squares.


## 46. Pure prediction commutator: exact neutral-integrator formula

Work on one LIN axis in state order (v,p,S,a).  The literal prediction is

    F(H,tau)=
      [1      0   0   phi_va
       H      1   0   phi_pa
       H^2/2  H   1   phi_Sa
       0      0   0   phi].                                (PC1)

The second-Abel boundary row has no AW component:

    A_B=[a_B,b_B,c_B,0],

with, for the selected first prediction of slab B divided by slab duration L_B,

    a_B=phi_va(h_B,tau_B)/L_B,
    b_B=phi_pa(h_B,tau_B)/L_B,
    c_B=phi_Sa(h_B,tau_B)/L_B.                             (PC2)

For a PURE prediction interval of total duration H, direct multiplication
gives the exact identity

    A_B F(H,tau)^T
      =[ a_B,
         H a_B+b_B,
         H^2 a_B/2+H b_B+c_B,
         0 ].                                               (PC3)

All OU-column coefficients cancel from PC3.  In particular the transported
row has zero AW component and is independent of phi,phi_va,phi_pa,phi_Sa of
the slab transport.  This is because A_B sees only the neutral v->p->S
integrator rows.

Therefore the pure-prediction commutator between adjacent slabs is

    C_B^pred =
      [ a_(B+1)-a_B,
        b_(B+1)-b_B-H_B a_B,
        c_(B+1)-c_B-H_B b_B-H_B^2 a_B/2,
        0 ],                                                (PC4)

axiswise, before attitude/reset coordinate transport.  No covariance, Kalman
gain, sigma_aw or R_S enters.

### Uniform coefficient bounds

The integral representations give

    0<=a_B<=h_B/L_B,
    0<=b_B<=h_B^2/(2L_B),
    0<=c_B<=h_B^3/(6L_B).                                  (PC5)

For h_B<=.006 and a nominal L_B=.1 these are

    a<=.06, b<=1.8e-4, c<=3.6e-7.

If H_B is the full ~.1-s slab duration, the polynomial transport terms in PC4
have ceilings

    H_B a_B <= .006,
    H_B b_B <=1.8e-5,
    H_B^2 a_B/2 <=3e-4                                     (PC6)

for H_B=.1.  Thus the p-component commutator is dominated by the coefficient
change Delta b and a .006 neutral-integrator term; the S component by
Delta c, 1.8e-5 and 3e-4 terms.

### Important consequence

Even with CONSTANT tuner and cadence, PC4 is not zero because A_B is the
first-prediction divided-difference row while F(H)^T transports across the
whole slab.  The common-tail cancellation therefore leaves a deterministic
neutral-integrator commutator of order H_B a_B.

This is not a long-horizon instability: the triangular neutral integrator
has exact polynomial structure.  Its Volterra sum must be combined across
slabs before norms.  Repeatedly summing the -H a and -H b-H^2 a/2 terms
telescopes into first/second divided differences of the boundary polynomial
reader.  Bounding PC4 slab-by-slab in l1 would unnecessarily pay O(M).

### Tuner variation

Only the coefficient differences

    Delta a_B, Delta b_B, Delta c_B                         (PC7)

carry tau/h/L variation.  Since
a=phi_va(h,tau)/L, b=phi_pa/L, c=phi_Sa/L, their source-uniform variation
must use the applied tuner chronology and sample/slab timing.  The large
neutral polynomial terms in PC4 should be telescoped exactly; only PC7 needs
a variation bound.

Thus C_pred has decomposed again:

    C_B^pred = C_B^neutral + Delta A_B,                     (PC8)

where C_neutral is an exactly summable triangular-integrator coboundary and
Delta A_B is the small tuner/cadence coefficient variation.

The next calculation is to telescope C_neutral over all slabs in the dual
Volterra sum, leaving endpoint polynomial terms, and separately bound the
total variation of (a_B,b_B,c_B) under the applied tau/h/L chronology.


## 47. Neutral Volterra telescoping closes; first-sample coefficient TV does not

Define the neutral three-state translation semigroup

    N(H)=
      [1 0 0
       H 1 0
       H^2/2 H 1],

so N(H1)N(H2)=N(H1+H2).  For the row
A_B=[a_B,b_B,c_B], section 46 gives exactly

    A_B N(H_B)^T
      =[a_B,H_B a_B+b_B,H_B^2 a_B/2+H_B b_B+c_B].          (NT1)

Hence the pure prediction commutator is the covariant difference

    C_B^pred=A_(B+1)-A_B N(H_B)^T.                         (NT2)

Let X_B be cumulative neutral time from boundary B to a fixed terminal
boundary, so X_B=H_B+X_(B+1).  Right-transport NT2 to that common terminal
frame:

    C_B^pred N(X_(B+1))^T
      = A_(B+1)N(X_(B+1))^T
        -A_B N(X_B)^T.                                     (NT3)

Therefore the Volterra sum telescopes EXACTLY:

    sum_(B=r..s) C_B^pred N(X_(B+1))^T
      = A_(s+1)N(X_(s+1))^T-A_r N(X_r)^T.                 (NT4)

The apparent slab terms -H a and -H b-H^2 a/2 do not accumulate at all.
For constant intrinsic A_B they reduce entirely to endpoint polynomial rows.
More generally, after placing every A_B in a common neutral frame, only
changes of the intrinsic A_B remain.  This closes item (1) without any
sqrt(number of slabs) or O(M) neutral charge.

### First-sample coefficient variation is not source-uniform

The current definition uses the first physical prediction of each slab:

    a_B=phi_va(h_B,tau_B)/L_B,
    b_B=phi_pa(h_B,tau_B)/L_B,
    c_B=phi_Sa(h_B,tau_B)/L_B.                             (NT5)

The declared sampling contract bounds h_B in [.004,.006] but does not impose
a total-variation bound on the sequence of sample intervals.  An admissible
sequence can alternate .004,.006 at every slab.  Since phi_va is strictly
increasing in h,

    TV(a_B)

then grows linearly with the number of slabs even with constant tau and L.
The same issue affects b_B,c_B.  Thus no horizon-independent useful uniform
TV(a,b,c) follows from the current chronology.

This is not a physical obstruction; it is an artifact of selecting ONE
sample to represent a .1-s slab.

### Correct block coefficient: use the whole slab

The first-sample row must be replaced by the exact block-integrated OU
coefficient.  For a slab with predictions k=1..n and no correction terms
(the latter remain separate ports), define the intrinsic physical leakage
row by summing each prediction's contribution transported through the later
PURE prediction maps.  Semigroup composition gives exactly the same
coefficient as one integrated-OU prediction over the total slab duration L_B
when tau is constant inside the applied-tuner cell:

    Abar_B =
      [ Phi_va(L_B,tau_B)/L_B,
        Phi_pa(L_B,tau_B)/L_B,
        Phi_Sa(L_B,tau_B)/L_B ].                            (NT6)

This coefficient is independent of the subdivision h_k and therefore immune
to sample-jitter TV.  It uses the entire slab rather than a representative
first sample.

If tau changes only at the slab boundary (the applied tuner chronology),
NT6 is exact.  If a commit can occur inside a chosen sync slab, split the
slab at that commit; no new physical assumption is needed.

The uniform sizes for L<=.1 are

    abar<=1,
    bbar<=L/2<=.05,
    cbar<=L^2/6<=.0016667,                                 (NT7)

which are larger than the first-sample coefficients but have controlled
chronology.  More importantly, their neutral covariant differences telescope
by NT3--NT4; only tuner/boundary variation remains.

### Tuner variation

For fixed L, Phi_va(L,tau)/L is monotone in tau and lies in [0,1].
The applied tau smoother remains in [.02,12].  A total-variation bound on
tau_applied over a 17-s word is NOT implied merely by this range: a bounded
sequence can oscillate.  The actual exponential adaptation law does constrain
per-commit motion, but without a bound on target-frequency variation its
total variation can still scale with the number of commits.

Therefore item (2), as originally phrased as TV(a_B,b_B,c_B), does not close
source-uniformly from current assumptions.

The correct use of NT4 is stronger: do not bound coefficient TV separately.
Keep the covariant differences C_B^pred inside the telescoping identity.
Then arbitrary intrinsic coefficient changes contribute only through endpoint
rows when transported by the neutral semigroup; correction-induced departures
are already separated as ports.  Any non-neutral effect of changing tau lies
only in the AW column of F, which A_B annihilates in PC3.

Thus for the PURE (v,p,S) prediction commutator, tuner variation also drops
out after covariant telescoping.  No TV(tau) theorem is required.

Conclusion:
- neutral prediction Volterra tail: exact endpoint telescoping CLOSED;
- raw first-sample TV(a,b,c): cannot be uniformly bounded and should be
  discarded;
- replace first-sample A_B by the exact whole-slab intrinsic row or retain
  the covariant-difference identity directly;
- the remaining nonlocal tail is no longer a prediction/tuner term.  It is
  only the correction/reset port departure from the neutral semigroup, which
  belongs in the already square-summed action framework.


## 48. Assembly audit: global neutral endpoint collapse is not quantitatively admissible

The remaining physical bound was to combine:
(i) the neutral endpoint polynomial charge,
(ii) local leverage chi_0,
(iii) globally square-summed correction/reset ports.

Before assigning numbers, the neutral endpoint normalization must be audited.

Section 47 transports every covariant difference to one common future frame:
A_B N(X_B)^T.  Algebraically this telescopes exactly.  But for a 17-s window,

    A N(X)^T =
      [a, X a+b, X^2 a/2+X b+c].                           (AS1)

Even with a<=1 for the whole-slab intrinsic row, the S coefficient at X=17
can be O(144).  Multiplying such a far-frame endpoint row by the declared
physical displacement bound P_max=8.1 is enormous.  Thus "telescope globally,
then apply the displacement norm" destroys the local primitive scaling.

The physical displacement bound is translation-local: it controls p(t) in
the chosen physical frame, while N(X)^T is the estimator integrated-chain
reader translation.  These cannot be paired after an arbitrary 17-s neutral
reader transport without carrying the corresponding physical polynomial
primitive transformation.

Therefore the correct deterministic estimate must keep each covariant
difference paired with the LOCAL displacement increment/endpoints before the
long neutral transport.  The Volterra identity is still useful for the
ACTION/source operator, but not as a single far-terminal physical endpoint
bound.

### Local endpoint charge

On one slab of duration L<=.1, the whole-slab intrinsic row satisfies

    abar=phi_va(L,tau)/L <=1,
    bbar<=L/2<=.05,
    cbar<=L^2/6<=.001667.                                  (AS2)

The second-Abel local endpoint term is therefore bounded with the local
physical displacement/velocity primitives without X^2 growth.  Adjacent
local endpoints cancel algebraically in the signed sum; correction-port
departures are handled by action.  Do not transport the physical endpoint
row through N(X) before applying P_max.

### Leverage and port terms cannot yet be assigned a source-uniform number

The local leverage is

    chi_B^2=A_B O_B^T Sigma_B^-1 O_B A_B^T.                (AS3)

This is finite on every frozen word.  A source-uniform chi_0 requires a
uniform enclosure of the literal slab observation design and its marginalized
source Gram Sigma_B.  Measurement-noise floors alone are insufficient if O_B
contains carried attitude/BA geometry with no source-uniform AG covariance
ceiling; using a raw H norm would repeat the G0/covariance circularity.

Likewise the globally square-summed correction/reset port action is useful
only after normalizing the terminal-AW minimum-action reader by a
SOURCE-UNIFORM action ceiling.  The existing joint-reader theorem explicitly
records that its historical action ceiling is still open source-uniformly.
The newly constructed AW-target reader inherits that same issue.

Hence there is presently no honest source-uniform numeric value for either
chi_0 times reader action or the global port-action contribution.

### What is closed and what remains

Closed analytically:
- physical acceleration appears only in OU prediction defects;
- two Abel steps are legal (V_max,P_max declared);
- neutral v->p->S prediction geometry telescopes exactly;
- sample-jitter/tuner TV need not be bounded;
- S and accelerometer innovations are exact Schur ports;
- correction ports are square-summed, not l1 accumulated.

Still open:
- a source-uniform normalized terminal-AW reader action on the retained
  shipping trace class;
- equivalently, a source-uniform local leverage/port bound for that reader.

This is the SAME historical-reader action obstruction already recorded for
the AG joint reader, now in the terminal-AW target.  It cannot be bypassed by
the neutral telescoping.

The next genuine theorem is therefore to prove a source-uniform action ceiling
for the terminal-AW minimum-action reader on one short (~.1-s) slab or on the
17-s qualified word, using the deterministic physical/BA/field constraints
already established here.  A frozen carried-word value can diagnose margin
but cannot certify it.


## 49. Terminal-AW action ceiling: reduction and causal-reader route

For an arbitrary terminal target row q, the joint minimum-action theorem gives

    B_q^* = Pi_q + Ttilde_q I_eff^-1 Ttilde_q^T,            (AWC1)

with the SAME

    I_eff=O_h^T Sigma^-1 O_h                               (AWC2)

as the AG reader.  I_eff depends only on the six unknown AG root columns and
the complete source/observation design, not on the terminal target.

Therefore a source-uniform floor I_eff>=mu I would immediately imply a
terminal-AW ceiling:

    B_AW^*
      <= (1+1/g) T_AW T_AW^T
         +(1+g) T_h,AW I_eff^-1 T_h,AW^T.                  (AWC3)

The source-driven AW terminal covariance T_AW T_AW^T is uniformly bounded by
the OU stationary/process construction plus the established nuisance/source
envelopes on a finite word, and T_h,AW is a finite product of bounded literal
mean maps.  Thus the only coercivity issue in AWC3 is mu>0.

But the current AG reader proof explicitly records that a source-uniform
moving-window I_eff floor remains OPEN.  Hence AWC3 does not create a new
proof of the requested AW ceiling; it reduces to the existing six-column
historical-reader obstruction.

### Why T_h,AW is not zero

Shipping prediction has no AG->AW mean block, but an accelerometer correction

    A=I-KH

does.  Its AW/AG block is

    A_AW,AG = -K_AW H_AG.                                  (AWC4)

Therefore an AG root perturbation can enter terminal AW through accepted
accelerometer corrections, and T_h,AW is generally nonzero.

### Causal cancellation identity

The same structure provides a more promising explicit trial reader.  At the
correction that creates the AG->AW transfer,

    delta a_AW^+ =
      delta a_AW^- - K_AW H_AG delta h^- - ...             (AWC5)

while the raw auxiliary observation contains

    y_i = H_AG delta h^- + H_n delta n^- + V_i s_i.         (AWC6)

Choose the local trial-reader block

    L_i^AW = - (future AW transport) K_AW                  (AWC7)

with sign according to the residual convention.  Its contribution L_i y_i
cancels exactly the AG-root term created by -K_AW H_AG at that event.  Repeat
this chronologically for every accelerometer correction.  Magnetic/S
corrections have no direct H_AW measurement row but can alter future
transport; their effects are included in the future AW transport in AWC7.

Backward induction then gives exact AG-root cancellation without solving a
six-column inverse: it is simply the Duhamel expansion of the corrected mean
recursion.

The action of this CAUSAL reader consists of:
1. transported AW process/sync factors;
2. gain-weighted accelerometer noise factors
       (future transport) K_AW V_i;
3. nuisance-root residual after AG cancellation;
4. S/mag correction effects only through their future transport.

This is exactly the actual-gain trial estimator of terminal AW.

### Source-uniform ceiling target for the causal reader

Because the minimum-action reader is no worse than any feasible reader,

    B_AW^* <= B_AW^causal.                                 (AWC8)

A source-uniform AW ceiling therefore follows if one can bound the causal
action directly.  Joseph gives a telescoping covariance budget for the
gain-weighted noise terms:

    K_i R_i K_i^T <= P_i^- - A_i P_i^- A_i^T,              (AWC9)

with subsequent future transport.  Summed chronologically, these measurement
terms plus process factors reconstruct the terminal AW covariance generated
from ZERO AG root and the bounded nuisance root.  Consequently

    B_AW^causal

is bounded by the diffuse-AG auxiliary Riccati terminal AW marginal, but this
statement alone is circular unless that marginal is bounded independently.

For AW specifically, OU prediction supplies phi<1 and fresh stationary
process covariance; S/accelerometer corrections are covariance-decreasing
globally.  The remaining source-uniform theorem can therefore be reduced to
an AW-marginal covariance upper bound with nuisance root bounded, WITHOUT an
AG root ceiling, by showing that AG-root covariance injected into AW by
corrections is canceled in the causal reader action rather than carried as
state covariance.

This is narrower than I_eff coercivity: construct B_AW^causal backward and
bound its residual nuisance/process/noise action directly using the already
proved nuisance upper comparison and OU stationary AW process budget.

The next calculation should write the exact backward action recursion for
this causal AW reader and test whether its action obeys a scalar/matrix
Lyapunov inequality driven only by the OU AW process and bounded nuisance
sector.  If yes, it proves B_AW,* without solving the open AG information
floor.


## 50. Causal terminal-AW reader: exact full-root action recursion

Let Y^+ be the terminal-AW trial-reader residual row immediately AFTER an
accepted correction. Choose the causal observation weight

    L_i = Y^+ K_i.                                         (CA1)

Then

    Y^- = Y^+(I-K_i H_i),                                  (CA2)

and Joseph gives exactly

    Y^- P_i^- Y^{-T} + L_i R_i L_i^T
      = Y^+ P_i^+ Y^{+T}.                                  (CA3)

Across a prediction x^+=F x^-+U w,

    Y^-=Y^+F,                                               (CA4)
    B^-=B^+ +(Y^+U)(Y^+U)^T,                               (CA5)

and a deterministic reset only transports Y.  Therefore, for ANY PSD root
covariance P_0 carried by the same literal covariance chronology,

    Y_0 P_0 Y_0^T
      + sum_process ||Y_k U_k||^2
      + sum_corrections L_i R_i L_i^T
      = u^T P_N u.                                         (CA6)

This is the useful causal-reader identity.  It includes the complete root
row; no root block is discarded.

### Important correction: the AG root does not cancel

The former section 51 claimed Y_root,AG=0 from the fact that the
correction-free skeleton has no AG->LIN/AW mean block.  That inference is
false.  A single correction is already a counterexample.  If q is a terminal
AW row and the pre-correction AG root injection is E_h, then

    Y^- E_h = q(I-KH)E_h = -q K H E_h                      (CA7)

whenever q E_h=0.  This is generically nonzero.  The observation term
L H E_h is part of the residual functional q x_N-L y; it must not be counted
a second time as a cancellation inside Y^-.

Equivalently, the chronological adjoint implemented by
`chronological_reader_adjoints` pulls a correction as q<-q-LH.  With the
causal choice L=qK this is exactly q(I-KH), not the correction-free skeleton
row.  Thus no historical six-column information inverse is avoided by a
fictitious zero AG row.

This correction does NOT destroy the sharp action ceiling below, because CA6
uses the full carried root covariance and telescopes to the terminal AW
principal covariance.

## 51. Sharp causal-reader action ceiling without AG-root cancellation

The literal shipping AW covariance synchronization assigns exactly

    P_aa <- Sigma_aw,stat,    Sigma_aw,stat <= 16 I.        (SC1)

Between syncs,

    P_aa^- = phi^2 P_aa^+
             +(1-phi^2) Sigma_aw,stat.                     (SC2)

Accepted Kalman corrections satisfy P^+<=P^- in Loewner order, hence cannot
increase the AW principal block.  Attitude resets leave AW unchanged.  A
later sync again assigns SC1.  Consequently every regular post-sync suffix
obeys

    P_aa(k) <=16 I.                                        (SC3)

Apply CA6 from an actual post-sync boundary with the FULL carried P_0,
including AG/LIN/BA blocks and all cross covariance.  For every unit terminal
AW direction u,

    B_AW,W^causal
      := Y_0 P_0 Y_0^T
         +sum_process ||Y_k U_k||^2
         +sum_corrections L_i R_i L_i^T
       = u^T P_N u
       <=16.                                                (SC4)

Therefore

    B_AW,W^min <= B_AW,W^causal <=16,
    sqrt(B_AW,*) <=4.                                      (SC5)

No AG covariance ceiling and no AG-root cancellation are needed.  The
arbitrary AG root is paid inside the single full-root action in SC4, and the
literal covariance chronology itself proves that the total paid action cannot
exceed the terminal AW marginal.

For the physical-transfer estimate this is harmless: the deterministic
physical primitive map has zero columns into root uncertainty.  Cauchy is
applied only after embedding that map into the SAME complete action space as
SC4.  Root action can consume part of the reader budget but cannot create an
extra physical charge.

## 52. Exact object required for the physical port enclosure

Let S_W be the complete normalized action-source space of the causal reader:
the full root factor, every prediction/process or covariance-sync factor, and
every accepted correction-noise factor, all in literal chronology.  Let
ell_W be its reader coefficient row.  SC4 states

    ||ell_W||_2 <=4.                                       (GP1)

After deterministic accelerometer elimination, keep the physical OU forcing
as the exact coboundary

    a_(k+1)-phi_k a_k,                                     (GP2)

perform the first Abel reduction with physical velocity and the second with
the declared displacement primitive, and telescope the neutral v->p->S
semigroup before taking norms.  This produces a linear map

    G_phys,W : P_W -> S_W,                                 (GP3)

where P_W is the finite-dimensional physical primitive coordinate space
(endpoint velocity/displacement plus the exact one-sided sampling remainder).
The root-source rows of G_phys,W are identically zero.  Accelerometer and S
columns are not independent controls: their entries are the Schur-completed
ports from the same-history correction chronology.

Define

    C_port = sup_W ||G_phys,W||_(P_W -> l2).               (GP4)

Then, with the dual norm on the bounded primitive coordinates,

    |R_ports,W|
       = |ell_W G_phys,W p_W|
       <=4 C_port ||p_W||_P.                               (GP5)

This is the correct single-budget formulation.  It neither requires nor
permits a second factor four for accelerometer, S, process/sync or root
families.

The old quantity

    ||D2 h||_(2,1) / sqrt(reader action)                    (GP6)

is NOT GP4.  It is a reader-dependent divided-difference diagnostic.  It
omits the physical primitive scaling, the complete normalized source columns,
and the Schur-completed port embedding.  It must not be called C_port,W or
compared with the 8.45e-3 target.

## 53. What is analytically closed and what remains

Closed, source-uniformly:
- accelerometer innovation is endogenous and is eliminated on the same
  physical history;
- GP2 is an exact OU coboundary for arbitrary applied phi;
- two Abel steps are legal under the already-declared V_max and P_max;
- the neutral LIN polynomial semigroup telescopes exactly, so no tuner or
  sample-jitter total-variation premise is needed;
- accelerometer and S corrections are retained through their complete Schur
  ports;
- the causal terminal-AW reader has the single full-root action ceiling SC5.

Still open: an explicit source-uniform numerical enclosure of GP4 over the
coupled shipping chronology.  Raw coefficient geometry is insufficient:
replacing the normalized source Gram by independent measurement floors gives
an O(1) bound, far above the available gravity-scale margin.  The enclosure
must therefore be performed on the complete Schur-normalized source operator
before norms, preserving common columns.

The next admissible analytical calculation is to write GP3 blockwise as the
composition

    bounded physical primitive
      -> exact Abel boundary/coboundary rows
      -> Schur-completed chronological source columns
      -> common normalized action space,                  (GP7)

then bound the largest singular value of the symbolic block operator using
only declared interval envelopes and exact covariance identities.  A carried
word may test this construction, but cannot supply any interval endpoint or
theorem constant.

## 54. Insert B_AW,*=16: joint action budget and remaining physical-to-port constant

The sharp reader ceiling gives

    ||reader||_action <= sqrt(B_AW,*) = 4.                 (JB1)

Do not allocate a separate factor 4 to local leverage and another factor 4
to correction/reset ports.  They are orthogonal/source blocks of the SAME
causal-reader action.  If the complete normalized source coordinates are
partitioned into slab-local observation ports, process/sync ports and root
nuisance ports, their squared reader coefficients satisfy

    sum_blocks ||ell_block||_source^2 <=16.                (JB2)

Let G_phys be the deterministic linear map from the bounded physical
primitive charges produced by the two-Abel construction into these normalized
source/port coordinates.  Then the complete port contribution obeys

    |R_ports| <= ||G_phys||_(primitive -> action-dual) * 4. (JB3)

Thus the former chi_0 and global correction-port constants should be combined
into ONE physical-to-port operator norm

    C_port := ||G_phys||.                                  (JB4)

The desired numerical physical charge is

    C_phys,total
      <= C_endpoint + 4 C_port + C_sampling/jerk.          (JB5)

No event-count or block-count factor occurs.

### Threshold comparison

For the worst field fraction sigma_w=1/5 and theta_max=6 deg, the previously
derived exact-degeneracy budget left approximately

    0.0338 m/s^2                                           (JB6)

at T=17 s after physical velocity endpoint, attitude, BA and sensor charges.
Therefore the robustified reader route closes that particular 17-s margin if

    C_endpoint,new + 4 C_port + C_sampling,new
       < 0.0338 m/s^2,                                     (JB7)

after subtracting only charges not already included in the earlier budget.

This inequality must avoid double counting: the two-Abel endpoint/sampling
terms replace the earlier crude physical-mean term; BA/attitude/sensor defects
already present in the 0.0338 calculation are not charged again through
C_port unless the port map represents an additional residual.

### What B_AW,*=16 does NOT supply

JB1 is an action normalization, not by itself a deterministic mean charge.
A numerical C_port is still required.  The earlier local leverage

    chi_B^2=A_B O_B^T Sigma_B^-1 O_B A_B^T

is one representation of the local part of C_port, but the complete
physical-to-port map should be formed before taking norms so that common
source columns and correction Schur cancellations are retained.

The infrastructure now has exactly the needed frozen-word matrices:
- terminal-AW causal/minimum-action reader;
- actual sync slabs;
- chronological adjoints;
- D2 rows;
- source/action factors.

For a source-uniform theorem, C_port must be enclosed over the admitted
shipping slab class.  Its coefficient geometry is independent of reader
normalization and carries the exact integrated-chain factors established
earlier.

### Immediate analytic bound available from coefficient geometry

For a <=.1-s local slab, the second-Abel map uses

    phi_va/L <=1,
    phi_pa/L <=.05,
    phi_Sa/L <=.001667                              (whole-slab form)

or the much smaller first-sample factors when used inside the exact local
coboundary.  Therefore a crude Euclidean bound on C_port is O(1), which would
give O(4 m/s^2) and fail JB7.  The proof requires the actual normalized
measurement/source Gram to obtain the expected much smaller leverage; raw
coefficient geometry alone cannot certify the 0.0338 margin.

Conclusion: B_AW,*=16 closes the reader normalization problem, but the final
numerical gravity-scale contradiction is NOT yet proved.  The single
remaining quantitative constant is C_port, the source-uniform induced norm of
the complete physical primitive -> normalized causal-reader port map.  A
frozen carried-word C_port may diagnose whether the 0.0338 margin is
realistic, but the theorem needs an enclosure over the shipping slab class.


## 55. Full action-dual norm is too coarse: restrict to reachable causal readers

There is a second structural correction to JB3--JB4.  If the action source
space is the orthogonal direct sum of root, process/sync and measurement-noise
factors, the unrestricted induced norm

    ||G_phys,W||_(primitive -> full action source)

cannot preserve a cancellation between two different source families.
The scalar two-port identity is enough to see this.  For one physical scalar
p entering two orthogonal ports with opposite sign,

    G p = (p,-p),

the full source norm is sqrt(2)|p|.  Yet a reachable reader with equal port
weights ell=(c,c) gives ell G p=0 exactly.  Taking ||G|| before imposing the
chronological relation between the reader weights destroys the cancellation.

The shipping signed-AW chronology has precisely this structure: physical
acceleration enters the OU prediction defect and the accelerometer correction
with linked opposite signs.  Process factors and accelerometer-noise factors
are orthogonal in the covariance action, so their cancellation is not visible
in the Euclidean norm of a reader-independent source vector.  Therefore the
phrase "complete map into normalized action-source coordinates" is
insufficient if action* means the unrestricted dual of that orthogonal direct
sum.

The correct object is the bilinear operator restricted to the REACHABLE
causal-reader subspace R_W generated by the literal backward recursion:

    C_port,W =
      sup_{ell in R_W, ||ell||_action<=1}
      sup_{||p||_P<=1} | ell G0_phys,W p |.                 (RR1)

Equivalently, if Pi_R,W denotes the action-orthogonal projection onto the
closure of reachable causal-reader rows,

    C_port,W = || Pi_R,W G0_phys,W ||_(P -> action).        (RR2)

This is not an arbitrary regression cone: R_W is defined exactly by the
shipping causal recursion L_i=Y^+K_i, prediction/reset pullback and the
terminal AW selector.  RR1 preserves the same-history process/correction
correlation by construction.  The global action ceiling still supplies
||ell||_action<=4 after rescaling.

The two-Abel reductions are now applied to the bilinear form BEFORE the
supremum in RR1.  In particular, the OU coboundary and neutral semigroup
telescope in ell G0 p; they are not replaced by the l2 norm of independent
process and measurement columns.

This identifies the actual source-uniform target:

    C_port = sup_W C_port,W,                                (RR3)

with W ranging over the retained coupled shipping chronology.  A proof of a
small RR3 must enclose the reachable backward-row recurrence jointly with the
Schur-completed correction ports.  Bounding the full orthogonal source norm
of G0 is a dead end because it discards the cancellation before the norm.



## 56. Exact one-sync-slab physical bilinear and joint accelerometer/S square

The previous source-space formulation can now be made exact without assigning
independent physical controls to corrections.

### Physical LIN lift

On one axis let

    chi=(v,p,S,a)^T

be the SAME physical kinematic history, with v'=a, p'=v, S'=p.  At prediction
j let F_j be the literal real-arithmetic shipping LIN transition and define

    delta_j = chi_(j+1)-F_j chi_j.                          (SB1)

For the ideal integrated OU transition this is exactly the Duhamel lift

    delta_j =
      int_0^h exp(A_L(h-s)) E_a
        [a_dot(t_j+s)+a(t_j+s)/tau_j] ds,                  (SB2)

and in particular

    E_a^T delta_j = a_(j+1)-phi_j a_j.                     (SB3)

Thus SB3 is the OU coboundary, while the v,p,S components of SB1 are its
integrated neutral companions.  Literal coefficient/rounding defects are
kept as the already-declared implementation/sampling remainder; they are not
folded into a new physical assumption.

Use physical error e=hat{x}-x_true.  Prediction contributes

    e_(j+1)^- = F_j e_j^+ - E_L delta_j + d_j^decl.         (SB4)

After deterministic accelerometer elimination, an accepted accelerometer
correction has NO free physical-acceleration input:

    e^+ = A_a e^- + d_a^decl,   A_a=I-K_a H_a.             (SB5)

At an S=0 pseudo update, however, the target is not the physical S.  With the
single fixed physical capture origin,

    S_true(t)=q(t)-q(T_c),

and shipping r_S=-S_hat=e_S-S_true (up to the fixed error-sign convention).
Hence

    e^+ = A_S e^- - K_S S_true + d_S^decl,
    A_S=I-K_S H_S.                                         (SB6)

The -K_S S_true term is real physical wave distortion and may not be treated
as fictitious pseudo-measurement noise.

### Exact slab bilinear

Take one actual covariance-sync slab B=[t_B,t_(B+1)] and let lambda denote
the reachable backward causal row for a terminal AW selector.  For every
accepted correction i put

    q_i = K_i^T lambda_i^+.                                 (SB7)

Ignoring only the separately declared d^decl terms, the exact physical
bilinear in the error recursion is

    B_B =
      -sum_(j in pred(B)) lambda_(j+1)^T E_L delta_j
      -sum_(s in S(B)) q_s^T S_true,s.                     (SB8)

There is no accelerometer a_phys term in SB8; it was eliminated in SB5.

Now telescope SB8 operation by operation against the SAME physical chi.
Across a prediction,

    lambda_j^T E_L chi_j
      -lambda_(j+1)^T E_L chi_(j+1)
      =-lambda_(j+1)^T E_L delta_j.                        (SB9)

Across an S update, the reachable-reader jump gives

    (lambda_s^- - lambda_s^+)^T E_L chi_s
      =-q_s^T S_true,s.                                    (SB10)

Across an accelerometer correction the physical error forcing is zero, but
the reader jumps.  Therefore its missing telescoping term must be restored
explicitly:

    (lambda_a^- - lambda_a^+)^T E_L chi_a
      =-q_a^T H_(a,L) chi_a.                               (SB11)

Magnetic/reset/sync operations have no direct LIN physical target; sync has
identity mean map.  Summing SB9--SB11 gives the EXACT one-slab identity

    boxed{
    B_B =
      lambda_B^T E_L chi_B
      -lambda_(B+1)^T E_L chi_(B+1)
      +sum_(a in acc(B)) q_a^T H_(a,L) chi_a
    }.                                                      (SB12)

For the literal accelerometer row H_(a,L) has only the AW LIN column, so

    H_(a,L) chi_a = R_wb,a a_phys,a                         (SB13)

in the fixed linearization convention.  Equations SB8 and SB12 are two
representations of the SAME bilinear.  SB8 is the correct representation for
the OU two-Abel reduction; SB12 proves that the S_true terms generated by the
regularizer are exactly the missing neutral-boundary terms.  In particular,
the declared P_AC bound must not be charged independently on every S event.

### First and second Abel inside the slab

Apply summation by parts to the AW component SB3 before norms.  The first step
moves the coboundary coefficient onto physical velocity.  The second moves
the resulting covariant neutral difference onto physical displacement.  If
the v,p,S components of SB1 are retained simultaneously, the neutral
semigroup identity

    N(H1)N(H2)=N(H1+H2)

makes the interior polynomial terms exactly the SB9 boundary difference.
The terminal S forcing SB10 cancels the S-boundary part.  Hence the complete
two-Abel result is not

    P_max * ||D2 h||_(2,1) + sum_s |q_s S_true,s|.

It is SB12 plus the one-sided jerk/implementation remainder.  This is the
precise reason the old D2 diagnostic cannot represent C_port.

### Joint accelerometer + S Schur quadratic

At each accepted correction i in {a,S}, let

    Omega_i = H_i P_i^- H_i^T + R_i >0.                    (SQ1)

For accelerometer Omega_i is 3x3; for each S axis it is the corresponding
pseudo innovation covariance (or retain the full 3x3 S block).  The exact
linear correction port is q_i^T r_i.  Completing the innovation square gives

    q_i^T r_i-r_i^T Omega_i^-1 r_i
      = -||Omega_i^-1/2(r_i-.5 Omega_i q_i)||^2
        +.25 q_i^T Omega_i q_i.                            (SQ2)

Stack the ACTUAL sequential innovations in slab chronology,

    q_B=(q_i)_(i in a,S),   r_B=(r_i)_(i in a,S).           (SQ3)

Because each Omega_i is the conditional innovation covariance after all
preceding operations, the sequential innovation factorization is block
diagonal in these coordinates:

    Omega_B = diag_chronological(Omega_i).                  (SQ4)

Therefore the exact JOINT slab square is

    Q_B(q_B,r_B)
      =-||Omega_B^-1/2(r_B-.5 Omega_B q_B)||^2
       +.25 q_B^T Omega_B q_B.                             (SQ5)

SQ5 is one quadratic form; there is no separate accelerometer factor and S
factor and no sqrt(event count).

The positive term in SQ5 must NOT be bounded independently by the reader
action.  The Kalman covariance identity gives, eventwise,

    K_i Omega_i K_i^T = P_i^- - P_i^+,                     (SQ6)

so

    q_i^T Omega_i q_i
      =lambda_i^{+T}(P_i^- - P_i^+)lambda_i^+.             (SQ7)

Thus the positive pieces of SQ5 are covariance-storage decrements.  They are
kept with the causal Joseph storage and telescope through the same chronology.
Replacing Omega_i by R_i, or bounding SQ7 separately, loses this cancellation.

### Telescope all sync slabs over 17 s BEFORE a norm

Let B=0,...,M-1 be the actual sync slabs in a 17-s regular word.  Summing
SB12 makes every internal physical boundary cancel because sync has identity
mean map and the backward row is the same row on the two sides:

    sum_B B_B =
      lambda_0^T E_L chi_0
      -lambda_N^T E_L chi_N
      +sum_(a in acc(W)) q_a^T R_wb,a a_phys,a
      +R_samp/impl.                                        (GW1)

Equivalently, summing the SB8 representation retains the OU coboundaries and
S_true forcings jointly; applying the two Abel steps and neutral-semigroup
identity gives GW1.  No slabwise norm is taken in either derivation.

Likewise the correction quadratic is formed globally first:

    Q_W = sum_B Q_B
        =-sum_i ||Omega_i^-1/2(r_i-.5 Omega_i q_i)||^2
         +.25 sum_i q_i^T Omega_i q_i.                     (GW2)

Use SQ7 together with prediction/process storage and the full-root Joseph
identity.  Only AFTER this telescoping is it legitimate to invoke the sharp
terminal-AW storage ceiling u^T P_N u<=16.  This is the single action budget;
the .25 sum q Omega q term is not a second budget.

### Remaining 17-s operator

GW1 is still not a numerical bound.  Its accelerometer sum and the two
physical endpoints are strongly correlated through the reachable q_i and
lambda_0.  Taking Cauchy on the accelerometer sum alone would recreate a
sqrt(N) loss, while charging |S_true|<=P_AC,max separately would undo SB10.

The exact remaining operator is therefore the signed map

    T_W(ell,chi)
      = lambda_0^T E_L chi_0
        -lambda_N^T E_L chi_N
        +sum_acc q_a^T R_wb,a a_a,                         (GW3)

with ell constrained by the literal causal backward recursion and with GW2
and the Joseph storage retained jointly.  The desired C_port is the induced
norm of GW3 AFTER the two-Abel primitive normalization, restricted to these
reachable rows.

This derivation closes the requested algebraic composition:
- one-sync-slab physical bilinear: SB8/SB12;
- exact joint accelerometer+S Schur form: SQ5;
- 17-s pre-norm telescope: GW1/GW2.

It does NOT yet prove C_port<8.45e-3.  The next quantitative lemma is a
source-uniform bound on the single signed endpoint+accelerometer functional
GW3 relative to the telescoped storage GW2, with V_max, P_max and the
one-sided jerk/implementation remainder inserted only once.  If that bound
cannot fit the remaining 0.0338 m/s^2 budget, GW3 is also the exact object
from which to construct an admissible analytical counterexample.


## 57. Exact signed field-axis rotation x physical-velocity cell functional

This section derives the low-frequency rotation term without replacing the
shipping attitude error by an arbitrary TV signal.

Fix the applicable unit world field direction b and put B=[b]_x.  On the
retained local chart factor the field-axis part of the relative attitude as

    M_b(theta)=Exp(theta B).                                (FV1)

The complementary magnetically observed attitude factor is retained in its
existing magnetic residual/action; FV1 is not an extra physical assumption.
For every actual accepted correction i define the EXACT field-axis increment

    kappa_i = theta_b(E_i^+) - theta_b(E_i^-),              (FV2)

where E_i^+ is obtained by the literal shipping quaternion injection
d_i=E_theta K_i r_i and theta_b(.) is the chosen smooth local twist
coordinate.  Thus

    Delta M_b,i
      = M_b,i^- [Exp(kappa_i B)-I].                         (FV3)

No replacement kappa_i=b^T d_i is made.  Instead define the exact nonlinear
twist/reset remainder

    rho_i^tw = kappa_i - b^T d_i.                           (FV4)

It is zero to first order and is charged by the existing retained-chart/reset
remainder machinery.

Consider one correction cell C=[t_L,t_R], including all sequential S,
accelerometer and magnetic corrections at their literal epochs.  Physical
velocity is continuous through estimator corrections.  Stieltjes integration
by parts gives EXACTLY

    int_C M_b(t) a(t) dt
      = M_b,R v_R - M_b,L v_L
        - int_C M_b B v d theta_b^cont
        - sum_(i in C) Delta M_b,i v_i.                    (FV5)

Hence the signed field-axis rotation x velocity functional is

    R_C(u) =
      int_C u^T M_b B v d theta_b^cont
      + sum_(i in C) u^T Delta M_b,i v_i,                  (FV6)

for the fixed terminal transverse direction u.  FV5 is the exact SO(3)
version of AX6.

### Literal correction decomposition

For the linear part of FV3 put

    c_i = u^T M_b,i^- B v_i,                               (FV7)
    d_i = E_theta K_i r_i.                                 (FV8)

Then

    u^T Delta M_b,i v_i
      = c_i b^T d_i + epsilon_i^SO3,                       (FV9)

where epsilon_i^SO3 is DEFINED by FV3/FV9 and contains both rho_i^tw and the
exact exponential remainder.  It is not silently dropped.

Partition each literal innovation by the actual measurement row.

Accelerometer:
    r_a =
      H_a,theta e_theta
      + H_a,AW e_AW
      + H_a,BA e_BA
      + nu_a,phys.                                         (FV10)

Here nu_a,phys is the same-history physical/sensor/model remainder after the
chosen true/nominal rotation convention.  Thus

    c_a b^T E_theta K_a r_a
      = J_a,theta + J_a,AW + J_a,BA + J_a,phys.            (FV11)

This is the requested AW/BA split.  AW and BA do NOT constitute separate
attitude corrections: they enter through the SAME accelerometer innovation
and gain.

Integral pseudo measurement:
    r_S = H_S e - S_true,                                  (FV12)

so

    c_S b^T E_theta K_S r_S
      = J_S,state - c_S b^T E_theta K_S S_true.            (FV13)

Even though H_S has no attitude column, K_S may have attitude rows through
carried cross covariance.  Therefore an S correction may rotate the mean and
must remain in the chronology.

Magnetic:
    r_m = H_m e + nu_m,                                    (FV14)

and

    c_m b^T E_theta K_m r_m = J_m,state + J_m,phys.        (FV15)

All nuisance cross covariance is inside the literal K_i.  No H_n deletion or
independent gain norm is used.

The continuous part of FV6 is generated by the literal gyro/gyro-bias
prediction.  Writing the exact field-axis rate in the chosen chart as

    d theta_b^cont = omega_b^err dt + d rho_pred,           (FV16)

retains physical gyro-bias history, estimator b_g, fast gyro residual and the
already-declared prediction/chart defect in one term.  Thus BA/AW/S enter the
rotation functional through their actual correction chronology, while gyro
bias enters both FV16 and future correction residuals.

### Joint Joseph/Schur square with the physical coefficient retained

The dangerous mistake would now be to bound |c_i| by V_max and then sum
correction action.  Keep c_i inside the correction port.  Define

    q_i = c_i E_theta^T b,                                 (FV17)
    z_i = K_i^T q_i.                                       (FV18)

The complete LINEAR correction part of FV6 is

    R_C^lin = sum_(i in C) z_i^T r_i.                      (FV19)

For the actual conditional innovation covariance Omega_i,

    z_i^T Omega_i z_i
      = q_i^T K_i Omega_i K_i^T q_i
      = q_i^T (P_i^- - P_i^+) q_i.                         (FV20)

The exact square completion is

    z_i^T r_i - r_i^T Omega_i^-1 r_i
      = -||Omega_i^-1/2(r_i-.5 Omega_i z_i)||^2
        +.25 z_i^T Omega_i z_i.                            (FV21)

Equations FV20--FV21 apply unchanged to S, accelerometer and magnetic events
and retain all BA/AW/nuisance cross covariance.  They are ONE chronological
Schur/Joseph construction.  The coefficient q_i depends on the same physical
v_i and M_b,i as FV5, so FV20 may not be replaced by V_max^2 times an
unweighted covariance decrement.

### Exact 17-s telescope before a norm

Let C_0,...,C_(N-1) be the actual correction cells of a 17-s word.  Summing
FV5 first gives

    int_W M_b a dt
      = M_b,N v_N - M_b,0 v_0
        - R_W^cont
        - R_W^lin
        - R_W^SO3,                                        (FV22)

where

    R_W^cont =
      sum_C int_C u^T M_b B v d theta_b^cont,              (FV23)

    R_W^lin =
      sum_(i in W) c_i b^T E_theta K_i r_i,                (FV24)

    R_W^SO3 =
      sum_(i in W) epsilon_i^SO3.                          (FV25)

Every internal physical endpoint cancels exactly.  No cellwise norm and no
raw TV(M_b) has been introduced.

Now augment the full-root causal factor chronology by ONE scalar accumulator
r whose update at correction i is

    r^+ = r^- + q_i^T K_i r_i,                             (FV26)

and whose deterministic coefficient q_i is frozen from the SAME physical
history.  Prediction carries the full root/process factors and also adds the
continuous signed supply FV23; reset/sync operations use their literal maps.
Because FV26 has no independent fresh covariance source, the covariance of
this augmented readout is obtained by the SAME root/process/sync/correction
factors as the shipping word.  In particular the correction block is the
single stacked row

    Z_corr,W = [ z_i^T Omega_i^(1/2) ]_(i chronological),  (FV27)

not separate accelerometer/S/magnetic norms.  Its Gram is

    Z_corr,W Z_corr,W^T
      = sum_i z_i^T Omega_i z_i                            (FV28)

only after every occurrence of a common upstream factor has already been
transported into the augmented full-root factorization.  The full factor
representation is therefore the safe object; FV28 alone is not promoted as
a telescope of weighted marginal decrements.

The complete signed 17-s functional is

    F_W(u) =
      u^T[M_b,N v_N-M_b,0 v_0]
      - R_W^cont
      - R_W^lin
      - R_W^SO3
      + R_W^Abel/samp
      + R_W^decl,                                          (FV29)

with R_W^Abel/samp the already-declared two-Abel quadrature remainder and
R_W^decl the existing BA/sensor/field/implementation terms, each included
once.

Let xi_W be the coefficient row of FV29 in the COMMON full-root Joseph factor
space, after substituting FV10--FV16 and the two Abel primitives.  Then the
first legitimate norm is

    |F_W(u)| <= ||xi_W||_2 ||s_W||_2
                + |endpoint/primitive remainders|.          (FV30)

The source-uniform quantity actually required is therefore

    C_rot,* =
      sup_(shipping-reachable W)
      ||xi_W||_2,                                          (FV31)

with q_i=c_i E_theta^T b generated by the same physical history.  This is a
joint root+accelerometer+S+magnetic+process operator.  It is NOT a total-NIS
bound, a TV bound, or a product V_max sqrt(sum Delta P_b).

### What this closes and what remains

Closed analytically:
- exact one-cell SO(3) Stieltjes identity FV5;
- literal S/accelerometer/magnetic chronology FV10--FV15;
- explicit AW and BA appearance through the SAME accelerometer port;
- exact nonlinear twist/reset remainder separation FV3--FV4/FV9;
- joint physical-coefficient Joseph square FV17--FV21;
- cancellation of all internal 17-s cell endpoints before norms FV22;
- one augmented full-root factor representation FV26--FV30.

Still OPEN: a useful source-uniform numerical enclosure of FV31, including
the continuous gyro/gyro-bias term FV23 and the nonlinear FV25 remainder,
tight enough that the COMPLETE nominal-AW inequality is below g sigma_w.
No 0.0338 m/s^2 residual margin is assumed available before those terms are
inserted.  Failure of a coarse enclosure is not a shipping counterexample.


## 58. Source-uniform enclosure of the continuous field-axis term and explicit envelope obstruction

Write the continuous field-axis error rate on the retained chart as

    dot(theta_b)^cont =
      b^T e_bg + b^T n_g + rho_chart,                       (CE1)

where e_bg is the SAME-HISTORY total residual gyro-bias channel from IMU BIAS,
n_g is the commissioned fast gyro residual, and rho_chart contains the
already-declared chart/prediction transport defect.  Physical angular rate
does not appear as an independent error source: it is common to truth and
nominal propagation.

The continuous contribution to the normalized physical mean is

    C_cont =
      (1/T) int_W c(t) [b^T e_bg(t)+b^T n_g(t)] dt
      +R_chart,
    c(t)=u^T M_b(t)[b]_x v(t).                              (CE2)

Since |c|<=|v|, the fast channel has the sharp envelope consequence

    |C_fast| <= V_max N_g.                                 (CE3)

For the slow residual bias, use the displacement primitive rather than the
velocity envelope.  Ignoring only the separately charged derivative of the
bounded rotation coefficient, scalar integration by parts gives

    (1/T) int v^T w e_bg dt
      = [p^T w e_bg]_0^T/T
        -(1/T) int p^T w dot(e_bg) dt
        -(1/T) int p^T dot(w) e_bg dt,                     (CE4)

with w=-[b]_x M_b^T u and |w|<=1.  Therefore the part not already assigned to
field-axis rotation/chart coupling obeys

    |C_bg,primitive|
      <= 2 P_max B_g/T + P_max D_g.                        (CE5)

At T=17 s, using the declared source envelopes,

    2 P_max B_g/T + P_max D_g
      = 2(8.1)(0.02)/17 + 8.1e-5
      = 0.0191399235294118 m/s^2,                          (CE6)

and

    V_max N_g = (5.5)(0.02)
              = 0.11 m/s^2.                               (CE7)

Thus even before R_chart and the nonlinear reset term,

    C_cont,envelope <= 0.129139923529412 m/s^2 + R_chart.  (CE8)

CE8 is a VALID source-uniform upper enclosure, but it is not useful for the
old 0.0338084 m/s^2 arithmetic remainder.

### Explicit admitted-envelope witness for the fast channel

The failure is not merely the looseness of CE3.  Consider on an integer
number of 2*pi-second cycles

    v(t)=5.5 sin(t) e_1,
    p(t)=-5.5 cos(t) e_1,
    a(t)=5.5 cos(t) e_1,
    jerk(t)=-5.5 sin(t) e_1,                               (CE9)

and choose the field axis/terminal direction so that
u^T M_b[b]_x e_1=1+O(theta_b), with commissioned fast gyro residual

    n_g(t)=0.02 sin(t) b.                                  (CE10)

Take e_bg=0 and no correction jump for this scalar envelope calculation.
Then

    theta_b(t)=-0.02 cos(t)+const,                          (CE11)

so the field-axis error amplitude is only 0.02 rad < 6 degrees, while all
physical primitive bounds are satisfied:

    |p|=5.5<8.1,
    |v|=5.5,
    |a|=5.5<8.8,
    |jerk|=5.5<100.                                        (CE12)

The normalized signed supply is

    (1/T) int v^T [b]_x^T M_b^T u (b^T n_g) dt
      = 5.5(0.02)/2 + O(0.02^2)
      = 0.055 + O(0.0004) m/s^2.                           (CE13)

This already exceeds 0.0338084.  CE9--CE13 are an explicit counterexample to
ANY proof that attempts to bound the continuous fast-gyro channel by the
declared independent physical/sensor envelopes and retained 6-degree tube
alone.

It is NOT yet a counterexample to the shipping theorem: the complete
accelerometer/S/magnetic correction chronology has not been solved for this
history.  A theorem may still close if those corrections cancel CE13 in the
same signed augmented factor row.  Such cancellation must be proved from
shipping reachability; it cannot be assumed from MAGNETIC SERVICE because a
rotation about the instantaneous field direction is precisely the
magnetically weak coordinate.

### Nonlinear SO(3) jump remainder

For an exact field-axis jump kappa,

    Exp(kappa B)-I-kappa B

has spectral norm

    r_exp(kappa)
      <= kappa^2/2                                         (CE14)

for the retained local branch.  Hence the exact jump remainder from FV9 obeys

    |epsilon_i^SO3|
      <= |v_i| kappa_i^2/2
         + |v_i| |rho_i^tw|,                               (CE15)

with a harmless refinement replacing kappa^2/2 by the exact trigonometric
remainder if desired.  Consequently

    |R_W^SO3|/T
      <= V_max/(2T) sum_i kappa_i^2
         +V_max/T sum_i |rho_i^tw|.                        (CE16)

The existing literal injection lemma gives eventwise

    |d_i|^2 <= NIS_i q_theta^T K_i Omega_i K_i^T q_theta,  (CE17)

but the current deterministic contract supplies neither a source-uniform
all-word sum of NIS-weighted decrements nor a signed bound on
sum |rho_i^tw|.  Therefore CE16 is finite on each carried word but has NO
currently proved useful source-uniform numerical constant.  Replacing it by
number-of-events times a pointwise maximum would resurrect the failed
sqrt(N)/TV architecture.

The correct escape is to keep the exact jump FV3 in the augmented nonlinear
functional, or prove a complete-word quadratic injection-action bound from
the literal accepted-correction chronology.  Until one of those is closed,

    C_rot,* <= 0.129139923529412
               + C_chart + C_SO3 + C_joint-correction      (CE18)

is the strongest simple source-uniform analytical enclosure supplied by the
declared envelopes, and it is quantitatively insufficient.

### Margin consequence

The previously quoted 17-s remainder

    Delta_old = 0.0338084 m/s^2

cannot be certified as spare margin for the corrected proof.  Even the
admitted-envelope fast-channel witness CE13 contributes about 0.055 m/s^2.
Therefore the test

    C_endpoint,new + 4 C_* + C_sampling,new < 0.0338

FAILS as a source-envelope argument before any positive endpoint, sampling or
nonlinear-reset charge is added.

This failure identifies the exact missing lemma rather than strengthening
MARINE MOTION:

    SAME-HISTORY FAST-GYRO CANCELLATION LEMMA:
    the signed fast-gyro term CE2 plus the literal
    accelerometer/S/magnetic correction jumps FV24 must admit a joint
    source-uniform bound substantially below V_max N_g.     (CE19)

CE19 must be derived from the literal estimator/measurement chronology with
the nominal AW/BA/S coupling retained.  If CE9--CE13 can be extended to make
those literal corrections satisfy every retained service/gate contract while
preserving a >available-margin signed residual, it becomes a genuine
shipping-reachable analytical counterexample.  Otherwise the mechanism that
prevents that extension is exactly the lemma needed to continue the proof.


## 59. Propagating the sinusoidal fast-gyro witness through the literal correction equations

The CE9--CE13 envelope witness can be embedded much further into the literal
sensor chronology.  This section separates what is exact from the remaining
self-consistent Riccati problem.

### Exact gyro/magnetic completion

Choose the committed unit field axis b and a smooth physical body-to-world
attitude

    R_true(t)=Exp(theta(t)[b]_x),
    dot(theta)(t)=-0.02 sin(t).                             (WC1)

Choose zero physical gyro-bias residual and the commissioned fast gyro
residual

    n_g(t)=+0.02 sin(t) b.                                  (WC2)

Then the measured gyro supplied to the nominal estimator is exactly zero:

    omega_meas = omega_true+n_g =0.                         (WC3)

A level nominal attitude therefore has zero prediction rotation.  Moreover a
rotation about the committed field leaves the magnetic sample invariant:

    R_true(t)^T B = B.                                      (WC4)

Thus with nominal qref=I and v2ref=B,

    r_mag(t)=0                                               (WC5)

at every callback.  Zero magnetic residual does NOT imply absence of magnetic
service.  The service information depends on H_m, the covariance and the
accepted cadence.  The historical exact field-axis obstruction already proves
that a 25-Hz zero-residual magnetic word can satisfy the one-second
MAGNETIC SERVICE floor; WC4 uses the same symmetry.  Therefore magnetic
corrections do not necessarily cancel the fast-gyro field-axis ambiguity.

### Exact accelerometer compatibility map

Let w(t) denote the nominal AW trajectory that the estimator would need on
this branch and take BA=0 first.  With lever arm zero, choose physical
acceleration

    a_phys(t)
      = g e_z + R_true(t)[w(t)-g e_z].                      (WC6)

Then the literal physical accelerometer sample is

    f_meas = R_true^T(a_phys-g e_z)
           = w-g e_z,                                      (WC7)

which is EXACTLY the nominal level accelerometer prediction for AW=w.
Hence

    r_acc=0                                                 (WC8)

whenever the nominal AW state equals w.  A nonzero BA trajectory simply
replaces w by w+ba in WC6/WC7.

For the desired large translation take the principal component

    w_1(t)=5.5 cos(t).                                      (WC9)

Because rotation about b can be chosen with b=e_1, this component is unchanged
by R_true.  The gravity compensation in WC6 has norm at most

    2 g sin(0.02/2)=0.1961297312 m/s^2.                    (WC10)

Therefore the instantaneous acceleration magnitude is below

    sqrt(5.5^2+0.19613^2)=5.50350 <8.8,                    (WC11)

before a tiny DC centering adjustment.  Its derivative is likewise far below
J_max=100.  The mean O(theta^2) vertical gravity defect in WC6 must be removed
by the corresponding O(theta^2) DC component of w (or BA); otherwise physical
velocity would acquire a secular drift.  This centering is below 0.001
m/s^2 and lies inside the declared envelopes.  The resulting p,v primitives
remain within the CE9 bounds plus O(0.2) transverse corrections.

Thus the physical kinematics and the literal gyro, accelerometer and magnetic
MEASUREMENT EQUATIONS do not exclude the witness.

### Why zero accelerometer innovation is not an invariant shipping trajectory

The estimator AW mean is not free.  Between corrections,

    w^-_(k+1)=phi_k w_k^+                                  (WC12)

in its AW component, with the corresponding exact v,p,S lift.  If r_acc=0 at
all epochs, accelerometer corrections cannot replenish the loss
(1-phi_k)w_k.  The S=0 correction can change AW through P_AW,S, but its input
is fixed by

    r_S=-S_hat.                                             (WC13)

Therefore WC8 can persist only if the homogeneous OU+S corrected map has the
required unit-frequency orbit.  There is no architectural identity asserting
this.

Allow the literal accelerometer innovation u_k.  Over one scheduled S
interval the exact frozen lifted scalar recurrence from the shipping code is

    x_(j+1)=A_S,j x_j + sum_l B_(j,l) u_(j,l),              (WC14)

where x=(vhat,phat,Shat,what), A_S contains the exact analytic OU propagation
and the actual S gain, and every B_(j,l) is the actual chronological
accelerometer gain transported through later operations.  Magnetic
corrections and BA coupling enlarge WC14 but do not change its affine form.

For a period-m lifted word write

    x_(j+m)=A_per x_j + B_per u_[j,j+m).                    (WC15)

A periodic AW target w_req sampled from WC6 is self-consistent iff

    (I-A_per)x_j = B_per u_per,                             (WC16)
    e_a^T x_k = w_req,k                                     (WC17)

at every accelerometer epoch, together with the literal innovation identity

    u_k =
      R_true,k^T(a_phys,k-g e_z)
      -Rhat_k^T(what_k-g e_z)
      -bhat_a,k                                             (WC18)

and the BA/S/magnetic state recurrences.

Equations WC14--WC18 are the exact finite-dimensional compatibility system
for the proposed periodic counterexample.  They show immediately that the
correction chronology does NOT NECESSARILY cancel the dangerous fast-gyro
term: the earlier frozen-word DC calculation has generic nonzero
accelerometer-to-AW gain, and no source identity forces B_per or the
unit-frequency transfer in WC16 to vanish.

Conversely, WC14--WC18 also show why CE9--CE13 is not yet a complete shipping
counterexample.  One must solve the actual periodic covariance/gain orbit,
because K_acc, K_S and K_mag are generated by that same orbit.

### Joint signed balance on a compatible periodic orbit

Suppose WC14--WC18 have a period-m solution.  Sum the exact attitude-error
balance over one period.  Since theta_b returns to its initial value,

    0 =
      int_period b^T n_g dt
      +sum_acc b^T E_theta K_a r_a
      +sum_S   b^T E_theta K_S r_S
      +sum_mag b^T E_theta K_m r_m
      +R_twist.                                             (WC19)

For WC4, r_mag=0, but the magnetic covariance still changes future gains.
Multiply the operationwise balance BEFORE summation by the physical signed
coefficient c_i.  The desired cancellation lemma would require the weighted
version

    int c(t)b^T n_g dt
      +sum_acc c_i b^T E_theta K_a r_a
      +sum_S c_i b^T E_theta K_S r_S
      +R_twist,c
      = small.                                              (WC20)

WC19 does NOT imply WC20 because c_i varies with physical velocity.  This is
the exact mathematical reason magnetic service plus bounded attitude error
does not by itself cancel the 0.055 m/s^2 supply.

Hence there is presently NO analytical necessity lemma forcing cancellation.
The only remaining discriminator is the self-consistent periodic
Riccati/mean system WC14--WC18.

### Status of the candidate

The witness now satisfies analytically:
- MARINE primitive amplitude and jerk limits, after the stated O(theta^2)
  centering;
- IMU BIAS with zero slow physical gyro bias;
- commissioned fast gyro residual exactly at its 0.02-rad/s envelope;
- gyro sample compatibility WC3;
- magnetometer sample compatibility WC4--WC5;
- accelerometer sample compatibility map WC6--WC8;
- field-axis error amplitude 0.02 rad, inside the retained 6-degree domain;
- a 25-Hz accepted zero-residual magnetic chronology of the same symmetry
  class for which the historical proof supplies one-second service.

Not yet proved:
- existence of the literal periodic covariance/gain orbit satisfying
  WC14--WC18 with the actual coupled tau,sigma,R_S adaptation chronology;
- all applied tuner/gate states on that orbit;
- the resulting exact weighted residual in WC20.

Therefore the next proof calculation is no longer another norm inequality.
It is a finite-dimensional PERIODIC RICCATI COMPATIBILITY problem: prove
WC14--WC18 has no solution uniformly over the retained coupled tuner
chronology, which would be the missing same-history cancellation/exclusion
lemma; or construct one exact/interval-enclosed solution and evaluate WC20,
which would complete the shipping-reachable counterexample.


## 60. Periodic Riccati/mean compatibility: analytic reduction and exclusion criterion

The remaining WC14--WC18 problem is not an arbitrary nonlinear fixed point.
For a prescribed periodic physical input history, the deployed front-end
variance/frequency channels are deterministic stable filters.  After their
transients, the tuner state satisfies a unique periodic forced recurrence.
The applied tau, sigma_aw, S cadence and R_S therefore form a fixed periodic
coefficient word U_*; they are not independent controls.

For this fixed word, covariance propagation is independent of the innovation
VALUES.  Let R_j denote one complete literal covariance step (prediction,
pending AW sync, scheduled S correction, accelerometer correction, magnetic
correction and resets in actual order).  Then

    P_(j+1)=R_j(P_j),    R_(j+m)=R_j.                       (PR1)

On the retained regular class all process/measurement covariances are
positive on their declared channels and the periodic word contains recurring
magnetic/accelerometer/S observations.  Hence any stabilizing periodic
covariance solution P_j^* is determined by U_* alone.  Existence/uniqueness
may be established by the standard finite-horizon Riccati monotonicity once
the already-required periodic detectability/stabilizability conditions are
inserted; no physical innovation can be chosen to alter P_j^*.

Freeze those literal periodic gains K_j^*=K(P_j^*,U_*).  The MEAN dynamics are
then an affine periodic linear system

    x_(j+1)=A_j x_j+B_j u_j+d_j,                            (PR2)

where u_j is the accelerometer measurement innovation and d_j contains the
fixed physical/BA/magnetic/S forcing not assigned to u.  Let

    Phi=A_(m-1)...A_0,                                     (PR3)

and define the one-period reachability matrix

    G=[A_(m-1)...A_1 B_0, ..., B_(m-1)].                   (PR4)

Periodic closure is exactly

    (I-Phi)x_0 = G u + d_per.                              (PR5)

This has the Fredholm criterion

    y^T(G u+d_per)=0
    for every y in ker((I-Phi)^T).                         (PR6)

If I-Phi is nonsingular there is NO Fredholm obstruction:

    x_0=(I-Phi)^-1(G u+d_per)                              (PR7)

for every periodic innovation word u.  Therefore contraction of the corrected
mean word actually favors existence of a periodic forced orbit; it does not
exclude the witness.

The physical compatibility equations close the loop.  Stack the actual
accelerometer rows over one period.  Because innovation is measured minus
predicted,

    u = z_phys - C x - c0.                                 (PR8)

Substitute PR7 into PR8.  The exact periodic compatibility equation is

    [ I + C (I-Phi)^-1 G ] u
      = z_phys-c0-C(I-Phi)^-1 d_per.                       (PR9)

Call the bracket D_per.  The witness is EXCLUDED iff either:
(a) PR6 fails in the singular case, or
(b) D_per is singular with the right side outside its range, or
(c) the unique solution violates a retained physical/gate/tuner condition.

There is no architectural reason for (a) or (b).  In fact D_per is the
finite-word closed-loop innovation sensitivity.  Positive R_acc means the
Kalman correction never imposes an exact algebraic measurement constraint;
for a finite covariance and finite word, the standard innovation map from
measurement sequence to innovation sequence is block lower triangular with
IDENTITY diagonal.  Therefore it is invertible.  Equivalently, chronological
Kalman filtering defines a bijection

    measurement word <-> innovation word                   (PR10)

for a fixed initial mean/coefficient word.  The periodic boundary condition
adds only the finite-dimensional root equation PR5.

This yields an important conclusion:

    PERIODIC RICCATI DYNAMICS ALONE CANNOT EXCLUDE
    THE FAST-GYRO WITNESS.                                 (PR11)

Any exclusion must come from the PHYSICAL admissibility of the unique
periodic closed-loop solution (amplitude/jerk/bounded primitives, BA
projection, tuner/gates, or magnetic-service qualification), not from a
missing mean fixed point.

### Constructive periodic solution map

For the proposed field-axis history the physical measurement word is an
explicit smooth function of the desired physical acceleration.  Define the
periodic closed-loop transfer from physical acceleration samples a to the
nominal AW samples by

    w_hat = T_aw,a a + t_aw,                               (PR12)

where T_aw,a is obtained by eliminating u and x_0 with PR5/PR8.  The exact
accelerometer compatibility construction WC6 requires

    a = g e_z + R_true(w_hat-g e_z)+a_free,                 (PR13)

with a_free reserved for the chosen principal physical oscillation/centering.
Substitute PR12:

    [I - R T_aw,a] a
      = g e_z - R g e_z + R t_aw + a_free.                 (PR14)

Thus the genuine counterexample/exclusion problem is one finite linear
periodic equation for the physical acceleration word, followed by deterministic
inequality checks.  If I-R T_aw,a is nonsingular, PR14 has a UNIQUE periodic
solution.  Again, nonsingularity produces the candidate rather than excludes
it.

A proof of impossibility therefore requires a quantitative statement that
EVERY PR14 solution violates at least one existing source envelope or gate.
No such statement follows from OU decay or S=0 architecture alone.

### What can be proved without numerical fitted constants

The previous frozen-word calculation already established that the
accelerometer-to-AW DC transfer is generically nonzero.  The same
chronological structure at frequency omega=1 gives a rational matrix transfer
in z=e^{i h} for a frozen subword.  Positive acceleration measurement noise
and finite gains make this transfer finite.  Hence, away from isolated
algebraic zeros of det(I-R T_aw,a), the implicit-function theorem gives a
locally unique periodic compatible solution depending continuously on
(tau,sigma,R_S,K).  The coupled tuner law restricts those coefficients to a
compact periodic path but supplies no identity pinning it to an algebraic
zero.

Consequently an ANALYTIC UNIVERSAL EXCLUSION of the periodic orbit cannot be
obtained from the present structural equations.  To promote a genuine
counterexample one still needs a rigorous interval enclosure of ONE actual
periodic coefficient/gain orbit and PR14 solution, followed by the declared
physical/service/gate checks.  Such interval evaluation is proof arithmetic,
not a fitted theorem premise.

This resolves the requested dichotomy at the structural analytical level:
the periodic Riccati/mean equations do not supply the missing cancellation
lemma; generically they admit a unique forced periodic solution.  The next
rigorous step is constructive interval certification of one shipping periodic
orbit, not another symbolic exclusion argument.


## 61. Constructive periodic certificate: first literal candidate excluded, lower-frequency refinement

To make the forcing exactly commensurate with the literal 5 ms IMU clock and
25 Hz magnetic callbacks, first choose omega=pi/3 rad/s (6 s period), with

    v_y=3.5 sin(omega t),
    n_g,x=0.02 sin(omega t),
    theta_x=(0.02/omega) cos(omega t),                      (PC1)

and physical acceleration a_y=3.5 omega cos(omega t).  The exact signed
fast-gyro supply remains

    <v_y n_g,x> = 3.5(0.02)/2 = 0.035 m/s^2.               (PC2)

The literal shipping wrapper was run from startup on this smooth exogenous
history, with measured gyro identically zero, B=(75,0,0), and the exact
world-to-body rotated accelerometer sample.  It reaches Live at step 7031 and
BA-active/refined operation at step 24016.  The committed field is exactly
(75,0,0), magnetic innovation is zero on the late word, and the late nominal
states remain finite (recorded AW norm <=4.70116, BA estimate norm <=0.235298).

However, this first candidate is EXCLUDED by an existing theorem condition:
over the late period the literal relative attitude error reaches

    0.130682 rad = 7.488 deg > pi/30.                       (PC3)

Thus it cannot serve as a counterexample inside the retained 6-degree local
domain.  This is a useful genuine exclusion: the failure is not periodic
Riccati solvability but the existing nonlinear-domain gate.

The signed supply PC2 is frequency independent at fixed velocity and gyro
residual amplitudes.  Therefore refine without changing any assumption to the
commensurate omega=pi/6 rad/s (12 s period, 2400 IMU samples, 300 magnetic
callbacks).  Then the analytic physical amplitudes are

    |p| <= 3.5/(pi/6) < 6.685 m,
    |v| = 3.5 m/s,
    |a| <= 3.5(pi/6) < 1.834 m/s^2,
    |jerk| <=3.5(pi/6)^2 <0.961 m/s^3,
    |theta_x| <=0.02/(pi/6) <0.03820 rad,                  (PC4)

while PC2 remains exactly 0.035 m/s^2.  All pre-compensation physical
amplitudes are strictly inside the retained envelopes.  This lower-frequency
candidate is the next interval-certificate target because the accelerometer
has substantially less dynamic forcing to misattribute to tilt.

No conclusion is promoted from the finite native replay.  The certificate
must still enclose the late 12-s tuner/covariance/mean orbit, prove the
6-degree bound and all gates on that orbit, and evaluate the complete weighted
functional including correction jumps.


## 62. Attempted source-uniform K_17<=250: horizon mismatch and exact obstruction

The proposed next step was to combine the existing aggregate world-frame
geometry, MAGNETIC SERVICE, jerk/sample fidelity, S-chain cancellation,
gyro-bias persistence and P_ba<=I/1600 on ONE 17-s word to prove

    sup_W kappa_nu(W) <= K_17 <=250.                        (K17-1)

That implication is NOT available from the current proved lemmas.

### Horizon audit

The relevant existing quantitative geometry has incompatible horizons.

1. Theorem G0 uses two accelerometer windows W1,W2 of length L=16 s
   separated by G=64 s, plus 2 s of endpoint room for the Lemma-T tube.
   Its certified constant s^2>=1.486786e-3 therefore belongs to an
   approximately 100-s construction, not a 17-s word.

2. The jerk/sampling theorem excludes the fixed-attitude sampled alias only
   over 32 s.

3. Its positive joint 3-D measured-vector information result is a 64-s
   statement.

4. The 17-s nuisance comparison is genuinely 17 s, and the S-chain
   cancellation and P_ba<=I/1600 are horizon-compatible, but they do not by
   themselves supply the missing slow AG quotient information.

Consequently inserting the G0/32-s/64-s constants into a 17-s diameter would
mix different words and is invalid.  Composition only allows information
actually contained in the chosen word.

### The G0 premises are not source consequences

Even on its proper long horizon, the explicit G0 number uses

    m_perp <= 2/5 m/s^2,
    u1     <= 6/5,                                         (K17-2)

for the NOMINAL force windows and an injection-free transported array.
Those are satisfied by the carried audits but remain unproved source-uniform
consequences of MARINE MOTION / IMU BIAS / MAGNETIC SERVICE.  The literal
injection-frame extension is also open.  Therefore G0 cannot currently be
promoted even on 100 s.

### What IS source-uniformly finite

For a radius-local retained class the exact variational reduction already
gives the right statement.  Eliminate nuisance coordinates by the Schur
complement and append the physical-kernel precision mu=1/c:

    G_red,mu = G_red + mu nu nu'.                           (K17-3)

On a compact same-history coefficient class, if

    x' G_red,mu(W) x >0                                    (K17-4)

for every unit slow x and every W in the class, continuity gives

    g_*(c,r)=min_(W,|x|=1) x'G_red,mu(W)x >0.              (K17-5)

Together with the exact S-chain fast elimination and bounded literal source
factors this implies a FINITE radius-local quotient diameter

    K_17(c,r)<infinity.                                    (K17-6)

This is an existence proof, not an explicit useful number.  Turning K17-5
into a numerical modulus requires quantitative same-word versions of:
magnetic transverse information after nuisance projection, accelerometer/BA
compatibility, gyro transport, and the kernel-row angle.  The current
published 17-s lemmas do not provide those constants.

### Why 250 cannot be claimed

The carried 16-s kappa_nu values near 196--198 demonstrate feasibility only.
There is no proved inequality placing every admissible word below them or
below 250.  In particular:
- source-uniform nominal signed AW-window statistics are OPEN;
- source-uniform literal injection-frame transport for G0 is OPEN;
- the complete slow Schur information floor after nuisance elimination is
  OPEN;
- MARINE MOTION T_E and theta_E are still symbolic, so no fixed 17-s
  excitation amount can be inserted.

Hence K_17<=250 is presently UNPROVED, and no valid algebraic combination of
the listed certificates yields it.

### Productive correction to the proof target

There are two legitimate paths.

A. Keep 17 s as the nuisance/covariance warm-up horizon, but prove a
RADIUS-LOCAL quotient bound K_17(c,r) jointly with retained-region invariance.
This matches the current finite-error architecture: source-only rho need not
be proved before the nonlinear radius.

B. Use a longer contraction superword whose horizon actually contains the
available geometry (at least 64 s, and 100 s for G0 as currently stated).
Then derive a same-superword K_T and compose finite-error supplies over that
horizon.  The carried diagnostics suggest the diameter improves strongly
with word length, but no carried number is promoted.

Path A is preferable if the goal is the existing 17-s recurring clock.
Its next lemma is not K_17<=250 outright but an explicit radius-local
variational floor

    G_red,mu(W) >= g_17(c,r) I                              (K17-7)

for every retained 17-s word, followed by conversion of g_17 and the
fast/S-chain factors into K_17(c,r).  This avoids horizon mixing and uses
the fact that storage itself bounds nominal AW/BA error on the candidate
region.


## 63. Adopt a 100-s contraction superword; retain 17 s only as nuisance warm-up

The proof clock is now separated from the estimator clock.  No shipping
schedule, tuner, correction cadence, assumption or quality gate is changed.

Let t_r be a regular A21 post-prediction root for which the existing 17-s
nuisance/root comparison has matured.  Define the contraction superword

    W_100=[t_r,t_r+100 s].                                  (SW1)

All literal operations inside W_100 are retained.  The 17-s theorem is used
only to establish the nuisance covariance/source class at t_r and at
intermediate roots; it is NOT asserted to contract the complete state.

### Why 100 s is the first convenient existing-proof horizon

A 100-s moving superword can contain, on the SAME history:
- a first 16-s accelerometer window W1;
- 64 s of separation;
- a second 16-s accelerometer window W2;
- the endpoint room already used by the Lemma-T tube construction, by placing
  the windows inside the available 100-s interval as in the existing G0
  proof;
- every 1-s MAGNETIC SERVICE subwindow;
- the complete 32-s jerk/alias exclusion interval;
- at least one complete 64-s joint measured-vector information interval;
- repeated 17-s nuisance warm-up blocks, exact S-chain cancellations and the
  proved P_ba<=I/1600 ceiling.

Thus no 32/64/100-s constant is imported from outside W_100.

For any prefix/suffix decomposition W_100=UV, the exact word identities give

    M_W=M_V M_U,                                            (SW2)

and the complete source covariance composes as

    S_W=M_V S_U M_V^T+S_V.                                 (SW3)

Information can only increase under adding observations to the word, while
the Riccati diameter is nonincreasing under informative prefixing.  Hence any
certified informative subword may be used inside W_100 without changing its
constant to a fictitious 17-s one.

### 100-s diameter theorem target

Define the same physical-kernel completion at the superword root,

    J_100,mu = J_100 + mu nu_0 nu_0^T,    mu=1/c_100.       (SW4)

Let Pi_100 be the known-root terminal covariance and P_nu,100 the terminal
covariance with only the kernel precision in SW4.  Put

    K_100 =
      sup_(admissible W_100)
      lambda_max(Pi_100^-1 P_nu,100).                      (SW5)

If the scalar kernel set is invariant between successive 100-s roots,

    nu_1^T P_nu,100 nu_1 <= c_100,                          (SW6)

then Corollary K gives

    rho_100 <= 1-1/K_100 <1.                               (SW7)

The recurring theorem is then stated on the 100-s roots.  Intermediate 17-s
roots are controlled by the exact every-prefix finite-error composition; they
need not contract individually.

### What existing geometry now supplies legitimately

On a complete moving W_100, the horizon issue is removed.  The existing G0
construction may be embedded as one slow AG information reader, the 32-s
sampling lemma may be used for alias exclusion, and the 64-s measured-vector
lemma may be used for physical-vector information.  Exact S-chain
cancellation removes the neutral LIN/AW root from the corresponding reduced
reader, while the 17-s nuisance comparison controls the remaining nuisance
action.  Gyro-bias persistence and MAGNETIC SERVICE act on the same W_100,
and P_ba<=I/1600 supplies the BA part of the kernel precision.

This gives the structural implication

    [G0 literal premises on W_100]
      + [same-history kernel return <1]
      => K_100<infinity
      => rho_100<1.                                        (SW8)

### Two obligations remain; changing the word does not erase them

The longer word fixes only the horizon mismatch.  It does NOT promote the
two source-open inputs of G0:

G100-1. NOMINAL FORCE PREMISES.
The explicit G0 number s^2>=1.486786e-3 uses

    m_perp<=0.4 m/s^2,    u1<=1.2                          (SW9)

on its two nominal 16-s windows.  The carried audits satisfy these bounds,
but no source/radius-uniform theorem currently derives them from the literal
AW loop.  On a retained storage ball the valid target is instead

    m_perp(r)<=m_phys,perp + 4 r + curvature(r),            (SW10)

or the stronger signed AW-loop estimate, evaluated on the same W_100.

G100-2. LITERAL INJECTION FRAME / KERNEL RETURN.
G0 is proved for the injection-free world array.  The literal reset/injection
transport and the equality case of the scalar kernel return must be handled
on W_100.  In particular the exact recurrence

    D_01(c)=d_perp+ell^2/(j0+||nu_0||^2/c)                 (SW11)

still forbids a finite invariant c if an admissible same-history exact-kernel
pair has unit persistence with d_perp>0.  A 100-s word gives substantially
more geometry with which to exclude that equality, but does not make the
exclusion automatic.

Therefore K_100 is not yet assigned a numerical value.

### Why this architecture is nevertheless strictly better

The 17-s K target required new quantitative geometry on a horizon shorter
than every existing physical-information theorem.  SW1 instead aligns the
contraction word with the proofs already available.  The remaining tasks are
now only:
1. convert the nominal signed-force/injection premises of G0 into
   radius-local same-history inequalities on W_100;
2. use the 32/64-s physical information plus recurring magnetic service and
   BA decay to exclude the unit-persistent exact-kernel equality on adjacent
   W_100 words;
3. evaluate the already-fixed canonical reader/action formulas to obtain an
   explicit K_100(c,r).

If K_100<=Kbar is obtained, the per-superword linear factor is

    rho_100 <= 1-1/Kbar.                                   (SW12)

For comparison only, a Kbar=250 certificate would give rho_100<=0.996 per
100 s; no such number is claimed here.

### Recurring finite-error composition

Let q_100=sqrt(rho_100).  At 100-s roots the retained-radius inequality is

    sqrt(V_(n+1)) <= q_100 sqrt(V_n)+E_100(r).              (SW13)

The scalar kernel ceiling propagates by SW6.  For every prefix
0<=tau<=100 s use the existing exact prefix composition

    sqrt(V(t_n+tau))
      <= G_tau(r) sqrt(V_n)+S_tau(r).                       (SW14)

Thus the longer proof word does not permit uncontrolled growth between
contraction epochs.  Regional practical stability follows once

    q_100 r+E_100(r) <= r,                                 (SW15)
    D_100(c,r) <= c,                                       (SW16)

and the prefix retained-domain inequalities close.  Startup/H18/release and
STILL/transition bridges remain the downstream obligations already present
in the single proof path.


## 64. Radius-local literal 100-s G0 attempt: exact reduction and obstruction

The 100-s superword removes the horizon mismatch, so obligation G100-1 can be
attacked on the SAME word.  The natural radius-local G0 substitution is

    m_perp(r) <= m0 + A1 r + A2 r^2,                       (RG1)
    delta_Q(r)<= q0 + Q1 r + Q2 r^2 + Q3 r^3,              (RG2)

where RG1 controls the two nominal 16-s force windows and RG2 controls the
literal injection-frame distortion relative to injection-free G0.

This calculation does NOT close with the current proved inequalities.

### Signed nominal AW mean

For normalized convex weights alpha_k on either 16-s window,

    mu_hat = sum_k alpha_k a_hat_k.                         (RG3)

The exact chronological AW mean can be written

    mu_hat = W0 a_hat_0 + sum_c W_c Delta_c,
    0<=W_c<=1,                                             (RG4)

with W_c generated by the literal OU/correction chronology.  Abel summation
moves Delta_c onto differences of W_c, but then requires

    TV(W)=sum_c |W_(c+1)-W_c|.                             (RG5)

The exact AW-loop identity instead gives

    sum_acc Gamma(e-eta)
      = e_0-e_N + sum xi
        -sum_pred[(1-phi)a_hat+Delta a],                   (RG6)

and sum Delta a telescopes physically.  However Gamma, xi and the effective
weights still depend on the adaptive covariance/sync chronology.  The
retained ball supplies pointwise

    |e_aw| <= sqrt(lambda_max(P_aw)) r
            <= 4 sqrt(1+eps_Q) r,                          (RG7)

but RG7 alone does not control the signed chronological mean because the
time-varying gains can rectify an oscillatory error.  The existing
sync-locked carried witness demonstrates that this mechanism is real.

Therefore no finite useful A1 in RG1 has yet been derived from the current
source assumptions.  A valid A1 requires a source-uniform gain/weight
variation inequality from the literal Riccati/sync recursion, or a direct
information proof that bypasses RG1.

### Literal injection frame

For each reset factor

    N_l=Exp(-X_l)(I+X_l/2)
       = I-X_l/2+R_l,
    ||R_l||<=|x_l|^3/6.                                    (RG8)

Lemma I* bounds the NET signed injection rotation by endpoint attitude errors
plus integrated gyro residual.  Hence the first-order signed term Q1 r is
radius-local.

But the exact ordered product remainder contains

    sum_l |x_l| |S_(l-1)|/4
      + exp(sum_l |x_l|^2/8)-1
      + sum_l |x_l|^3/6,                                   (RG9)

where S_n=sum_(l<=n)x_l.  The Loewner injection lemma controls each event,

    x_l x_l^T <= NIS_l P_theta,l,                          (RG10)

but the current theorem has no source-uniform complete-100-s bound on

    sum_l |x_l|^2                                          (RG11)

or the partial-sum weighted quadratic term in RG9 before contraction/action
is known.  Bounding RG11 from the desired contraction would be circular.
Thus Q2,Q3 in RG2 are not currently certified.

### Direct-information formulation avoids the scalar circularity

Consequently RG1--RG2 are NOT adopted as theorem premises.  The valid literal
100-s object is the nuisance-reduced information itself.

Construct the raw auxiliary record on W_100, apply the exact S-chain
annihilator to the LIN/AW/root/sync nuisance columns, retain every literal
reset factor N_l and nominal AW coefficient, and whiten with the FULL reduced
source covariance Sigma_red.  For slow coordinates x_s=(theta,b_g,b_a
compatibility quotient), define

    G_red(W)
      = O_s^T Sigma_red^-1/2
          (I-P_f)
        Sigma_red^-1/2 O_s,                                (RG12)

where P_f projects onto the whitened nuisance range.  Append the physical
kernel precision:

    G_red,mu(W)=G_red(W)+(1/c) nu nu^T.                    (RG13)

The exact radius-local 100-s G0 target is now

    G_red,mu(W) >= g_100(c,r) I >0                         (RG14)

for every retained same-history moving W_100.  This formulation has:
- no separate nominal AW mean assumption;
- no expansion of the reset product;
- no gain-TV bound;
- no event-count injection norm sum;
- exact S-chain nuisance cancellation;
- all source correlations retained before whitening.

### Remaining missing modulus

RG12--RG14 still do not close from the current assumptions.  After nuisance
projection, the missing quantitative implication is precisely

    physical MARINE attitude/gravity excitation
       => literal nominal accelerometer slow row
          stays a positive Sigma_red^-1 distance
          from the nuisance + magnetic-axis compatibility span.             (RG15)

MAGNETIC SERVICE supplies transverse magnetic information; Lemma T supplies
gyro chronology once the accelerometer-window geometry is positive; the
kernel row removes the final physical tilt/BA line.  But MARINE MOTION
constrains TRUE attitude/gravity while O_s contains the ESTIMATOR NOMINAL
force a_hat-g e_z.  The current assumptions and proved storage marginals do
not yet supply the positive physical-to-nominal separation in RG15.

Thus changing to 100 s successfully makes all horizons compatible, but the
literal G0 floor remains blocked by ONE physical-to-nominal accelerometer
modulus.  This is exactly the gap already identified by the direct-information
route in the corrected-word proof.

The next useful theorem is therefore RG15 itself, not another scalar
m_perp/u1 estimate.  It must use the SAME-HISTORY closed-loop accelerometer
identity and retained storage ball to show that a sequence with vanishing
reduced accelerometer distance would force either:
(a) the physical attitude/gravity span to vanish, contradicting MARINE
MOTION on the complete moving window; or
(b) nonzero correction/process action already counted in Sigma_red, yielding
a positive reduced information charge.

If neither implication can be proved, the limiting sequence supplies the
candidate admissible compatibility trajectory that blocks the 100-s
contraction theorem.


## 65. Physical-to-nominal separation lemma: exact compatibility and refutation

Assume the nuisance-projected literal accelerometer information tends to zero
on a retained same-history sequence.  The requested implication was that the
true MARINE attitude/gravity span must then tend to zero unless positive
process/correction action remains.  The exact accelerometer/BA compatibility
equations show that this implication is FALSE as a consequence of attitude
span alone.

### Two-epoch exact Schur compatibility

Take two applied accelerometer epochs inside one complete moving excitation
window and rotate both residual equations to common world coordinates.  After
the exact S-chain has removed the free LIN/AW root, write the slow
attitude/BA contribution as

    d(z)=D z,    z=(theta_0,b_a,0),                         (PN1)

    D=[ C_0                 B_0
        C_1 T_theta    B_1 phi_b ],                         (PN2)

where
- C_i=-[f_hat_i]_x is the LITERAL nominal specific-force attitude row;
- B_i is the transported invertible BA row;
- T_theta is literal attitude transport between the epochs;
- phi_b is the homogeneous BA decay.

With inherited AW/process nuisance retained, exact elimination gives

    Q_acc,red(z)=d(z)^T S_a^-1 d(z),                        (PN3)

where S_a is the positive two-epoch residual covariance containing actual
accelerometer noise plus the transported AW root/process action.  Thus zero
reduced accelerometer information is equivalent to D z=0 in the zero-action
limit.

Eliminate BA from PN2.  The first row gives

    b_a,0=-B_0^-1 C_0 theta_0.                              (PN4)

Substitution into the second gives the exact relative compatibility operator

    L_2 theta_0=0,                                         (PN5)

    L_2 :=
      C_1 T_theta
      -phi_b B_1 B_0^-1 C_0.                               (PN6)

(With the orthogonal BA convention B_0^-1=B_0^T.)  Therefore the nontrivial
two-epoch compatibility class is precisely ker L_2.

For finite action, completing the square gives, for 0<eta<1/phi_b^2,

    |D z|^2
      >= eta/(1+eta) |L_2 theta_0|^2
         +(1-eta phi_b^2)|e_0|^2,                          (PN7)

where e_0=C_0 theta_0+B_0 b_a,0.  Consequently a positive physical-to-nominal
modulus would require a source-uniform positive singular floor for L_2 on the
relevant quotient.

### Pullback to true gravity does not supply that floor

Let R_true,i be the true attitudes.  MARINE MOTION supplies a span between
some attitudes in every complete excitation window.  On a retained ball the
nominal gravity directions remain close to the true ones, schematically

    Delta_g,nom(r)
      >=2 sin(Delta_R/2)-4 sin(theta_err,max(r)/2).         (PN8)

But C_i is NOT the gravity cross-product map.  It uses the nominal specific
force

    f_hat_i = a_hat_i-g e_z                                (PN9)

in world convention.  Translational acceleration is an admitted physical
degree of freedom.  It can compensate the changed gravity direction so that

    C_1 T_theta
      =phi_b B_1 B_0^-1 C_0                               (PN10)

on the selected epochs, i.e. L_2=0, while R_true,1 differs from R_true,0.

The existing bounded collinear/same-cell constructions exhibit exactly this
mechanism: nonzero attitude span, bounded displacement/velocity/acceleration/
jerk and force directions chosen to be compatible at selected correction
times.  MAGNETIC SERVICE excludes a particular sparse magnetic cadence but
does not convert physical attitude span into a uniform two-epoch
specific-force separation.

Therefore

    MARINE attitude span
      -/-> sigma_min^+(L_2)>0,                              (PN11)

and hence

    MARINE attitude span
      -/-> positive nuisance-projected accelerometer
           information                                     (PN12)

under the current assumptions.

This refutes the proposed physical-to-nominal separation lemma in its
two-epoch/sensor-family form.  It does NOT produce a zero-action complete
shipping word: magnetic rows, all accelerometer epochs, S observations,
process penalties and terminal forgetting still act jointly.

### Correct replacement: complete-word joint compatibility

The proof must use the complete 100-s corrected word.  Stack exactly

    y=O_s x_s+O_f x_f+A s,
    x_N=T_s x_s+T_f x_f+B s.                               (PN13)

Eliminate nuisance root x_f with the full whitened Schur projector and retain
all fresh source factors once.  A sequence with vanishing TOTAL joint action
must simultaneously satisfy:
1. zero fresh LIN/AW/BA/AG process action;
2. four-S homogeneous LIN/AW compatibility;
3. every magnetic-service row;
4. EVERY accelerometer compatibility equation with one deterministically
   propagated BA root;
5. terminal persistence/forgetting.

Zero fresh action rigidifies the nuisance trajectory: translational/AW/BA
mimics cannot be retuned independently at the epochs used in PN10.  The S
rows kill the free homogeneous LIN/AW trajectory qualitatively.  Magnetic
rows restrict the attitude/gyro trajectory to at most one transported
field-compatible line.  Then all accelerometer equations determine at most
one common BA compatibility line.

Explicitly, after magnetic reduction let

    theta_k=F_k theta_0,
    b_a,k=phi_b(t_k) R_ba,k b_a,0.                          (PN14)

At every accelerometer epoch zero loss requires

    J_att,k F_k theta_0
      +R_ba,k phi_b(t_k)b_a,0=0.                            (PN15)

If the pulled-back magnetic lines intersect trivially, theta_0=0 and PN15
gives b_a,0=0.  Otherwise write theta_0=lambda theta_hat_0.  Then a nonzero
solution exists iff

    q_k :=
      phi_b(t_k)^-1 R_ba,k^T
      J_att,k F_k theta_hat_0                               (PN16)

is IDENTICAL at every accelerometer epoch.  If so the complete sensor
nullspace is the single word-dependent compatibility line

    nu_W=(theta_hat_0,0,...,0,-q_W).                        (PN17)

If the q_k are not all identical, the complete slow nullspace is trivial.

Thus the correct qualitative conclusion is

    Null(complete 100-s joint action)
       subset span(nu_W),                                  (PN18)

not that physical attitude span alone gives a positive accelerometer floor.

The rank-one Riccati-diameter theorem is basis-free and can use this
word-dependent nu_W.  The remaining quantitative problem is a RELATIVE
complete-word inequality between terminal persistence and total action on
the quotient of span(nu_W), plus the adjacent-superword scalar return.  This
is the finite-horizon detectability problem already identified in the
corrected-word proof.

### Consequence for the 100-s G0 plan

The longer superword remains useful because all physical/service lemmas now
live on one history, but G0 cannot be completed by deriving a standalone
physical-to-nominal accelerometer modulus from MARINE attitude span.  That
route is closed.

The next controlling calculation must instead attack the complete-word
relative action:

    H_eff^T Pi_100^-1 H_eff
      <= K_rel S_q,                                        (PN19)

where S_q is the nuisance-eliminated information Schur complement transverse
to the word-dependent kernel and H_eff is the corresponding terminal excess
map.  A vanishing S_q is harmless if H_eff vanishes at the same rate.  This
relative inequality is exactly what the Riccati diameter needs and is weaker
than a uniform Euclidean accelerometer-information floor.


## 66. Relative quotient inequality from the complete 100-s conditional Gaussian word

Fix one literal 100-s superword and perform the SAME nuisance elimination and
word-dependent kernel split as PN13--PN18.  In quotient coordinates v, write

    y = O_q v + A s,                                       (RQ1)
    x_N = T_q v + T s,                                     (RQ2)

with Sigma=A A^T>0 the complete auxiliary source covariance.  All prediction,
AW sync, accelerometer, S, magnetic, BA and reset chronology is already inside
O_q,A,T_q,T.

The known-root terminal covariance is the shorted source covariance

    Pi =
      T [I-A^T Sigma^-1 A] T^T.                            (RQ3)

Define

    J_q = O_q^T Sigma^-1 O_q,                              (RQ4)

    Ttilde_q =
      T_q - T A^T Sigma^-1 O_q.                            (RQ5)

RQ5 is the terminal root image after optimally reusing the SAME fresh-source
coordinates to explain the data.  The diffuse quotient adds exactly

    Ttilde_q J_q^dagger Ttilde_q^T                         (RQ6)

when the Moore-Penrose range condition is satisfied.

### Kernel-line Schur elimination

Before quotienting, split the augmented root information with normalized
word kernel n and quotient q:

    J_mu =
      [ j_nn+mu   j_nq^T
        j_nq      J_qq ],    mu=1/c.                       (RQ7)

The quotient Schur information is

    S_q =
      J_qq - j_nq j_nq^T/(j_nn+mu).                        (RQ8)

In Pi-whitened terminal coordinates split the excess map

    H=[h_n,H_q].                                            (RQ9)

Eliminating the same kernel coordinate gives

    H_eff =
      H_q - h_n j_nq^T/(j_nn+mu).                          (RQ10)

Therefore the exact quotient contribution to the kernel-bounded Riccati
diameter is controlled by

    H_eff^T H_eff <= (K_rel-1) S_q.                        (RQ11)

Restoring unwhitened coordinates,

    H_eff^T Pi^-1 H_eff
      <= (K_rel-1) S_q.                                    (RQ12)

The SHARP fixed-word constant is

    K_rel(W)-1 =
      lambda_max[
        S_q^dagger/2
        H_eff^T Pi^-1 H_eff
        S_q^dagger/2 ],                                    (RQ13)

with value infinity exactly when

    Null(S_q) not subset Null(H_eff).                       (RQ14)

Thus the requested relative inequality is not an additional relaxation: it
is exactly the quotient generalized eigenvalue in the complete conditional
Gaussian word.

### Exact range test

After kernel shorting, rewrite the quotient model again in the form RQ1--RQ2.
For any v in Null(J_q),

    0=v^T J_q v=||Sigma^-1/2 O_q v||^2
      => O_q v=0.                                          (RQ15)

Then RQ5 gives

    Ttilde_q v=T_q v.                                      (RQ16)

Hence the finite-extension condition is exactly

    Null(O_q) subset Null(T_q).                            (RQ17)

This is a finite-horizon detectability condition: every quotient root
direction invisible to the COMPLETE 100-s observation record must also have
zero deterministic terminal image.

The complete-word zero-action classification proves RQ17 for each fixed
nondegenerate word after quotienting by its exact compatibility line.
However pointwise injectivity is insufficient for a UNIFORM K_rel.  Along a
rank-changing sequence, singular values of O_q may tend to zero while
T_q v remains O(1).  Therefore compactness alone cannot bound RQ13.

### Regularized backward-reader representation

The relative ratio can be represented without forming small singular values.
Initialize the terminal quotient residual in Pi-whitened coordinates and run
the literal word backward.  At each operation:
- prediction: Y<-Y F and add the common process factor Y U;
- accepted correction: choose the observation reader block L_i and set
  Y<-Y-L_i H_i while adding L_i V_i;
- reset: Y<-Y G_i;
- PSD AW sync: identity mean pullback with its source factor retained.

Choose all L_i jointly by the regularized normal equations corresponding to
J_mu.  The minimum total fresh-source action is exactly the quadratic
numerator in RQ13.  The denominator is the quotient observation action S_q.
Thus a source-uniform reader estimate

    Action_terminal(v)
      <= C_det Action_observation(v)                       (RQ18)

for every quotient root v is equivalent to

    K_rel <= 1+C_det.                                      (RQ19)

No Euclidean information floor is required.

### What the current structural lemmas imply

The 100-s word gives the following zero-action chain:
1. zero fresh LIN/AW/BA/AG source action rigidifies nuisance trajectories;
2. four S rows kill the zero-source homogeneous LIN/AW root qualitatively;
3. magnetic service plus gyro transport restricts AG to at most one
   transported field-compatible line;
4. all accelerometer rows with one propagated BA root leave at most the
   word-dependent compatibility line nu_W;
5. quotienting by nu_W leaves no fixed-word zero-action root direction.

Therefore the only possible failure of a UNIFORM RQ18 is a NEAR-null sequence
whose observation action tends to zero faster than its terminal persistence.

### Remaining quantitative target

The relative quotient inequality has now been DERIVED exactly, but no finite
source-uniform numerical C_det is yet proved.  The next analytical obligation
is the near-null persistence lemma:

    for every retained sequence (W_n,v_n) with |v_n|=1,
    Action_observation(W_n,v_n)->0
      => Action_terminal(W_n,v_n)->0
         at a uniform linear rate.                          (RQ20)

Equivalently, construct C_det from the literal prediction/process structure.
This is weaker than proving sigma_min(O_q)>0 and is the shortest remaining
route to K_100<infinity.

A useful decomposition of RQ20 is by the last operation at which a normalized
near-null direction has appreciable amplitude.  If it persists to the
terminal state, backward propagation across the final 17-s nuisance-regular
suffix must either:
(a) generate S/process action through LIN/AW;
(b) generate magnetic/gyro action through AG;
(c) generate accelerometer/BA action outside nu_W; or
(d) lie asymptotically in the next word's compatibility kernel.
Cases (a)--(c) charge S_q; case (d) is exactly the adjacent-superword scalar
kernel-return problem.  Thus relative quotient detectability and kernel
return are the two complementary pieces of the recurring contraction proof.


## 67. Near-null persistence on the final 17-s suffix: dichotomy and reachability blocker

Let W_n be retained 100-s moving superwords and v_n normalized quotient root
directions, after shorting the exact word kernel nu_(W_n), such that

    A_obs,n :=
      v_n^T S_q(W_n) v_n ->0.                              (NP1)

Assume their terminal persistence does NOT vanish:

    A_term,n :=
      ||H_eff(W_n)v_n||_(Pi^-1)^2 >= eps_0>0.              (NP2)

We derive the consequences on the final regular 17-s suffix Z_n.

Because the complete joint action is a sum of squared COMMON source-factor
coefficients before any familywise norm, NP1 implies every individual
nonnegative block contribution on Z_n tends to zero after the optimal
nuisance/kernel shorting.  Pass to a convergent subsequence of the retained
coefficient histories and normalized directions.

### Case 1: LIN/AW amplitude persists on the suffix

If a nonzero homogeneous LIN/AW component survives on Z_n, zero fresh
LIN/AW/sync action rigidifies it to one deterministic homogeneous chain.
Zero S corrected loss at four distinct applied S epochs then forces the
root coefficients of

    S(t)=S0+p0 t+v0 t^2/2+a0 psi_tau(t)                    (NP3)

to vanish qualitatively because {1,t,t^2,psi_tau} is a strict Chebyshev
system for every finite tau>0.  Therefore a terminal-persistent near-null
sequence cannot retain an O(1) LIN/AW quotient component without paying
S/process action.  Quantitatively, a uniform C_det would require a lower
singular modulus for this four-S map in the ACTUAL normalized source metric;
the earlier explicit determinant floor is retracted and cannot be reused.

### Case 2: AG component persists away from magnetic compatibility

If the transported AG component has O(1) distance from the pulled-back
magnetic-compatible line at some qualified service interval, MAGNETIC SERVICE
and positive magnetic R force O(1) observation action.  Hence NP1 implies the
AG trajectory approaches the deterministic magnetic-compatible class on the
final suffix.  Zero fresh AG action simultaneously removes independent
gyro-bias retuning; Lemma-T chronology then leaves at most the transported
field-compatible attitude/gyro line.

A uniform rate again needs the residualized magnetic/gyro modulus; the
previous unprojected leverage inequality is retracted.

### Case 3: accelerometer/BA incompatibility persists outside the word kernel

With LIN/AW gone asymptotically and AG restricted to the magnetic-compatible
line, accelerometer loss reduces to one propagated BA root:

    J_att,k F_k theta_0
      +R_ba,k phi_b(t_k)b_a,0 ->0.                          (NP4)

If the vectors

    q_k =
      phi_b(t_k)^-1 R_ba,k^T
      J_att,k F_k theta_hat_0                               (NP5)

do not converge to one common q, at least one accelerometer row has positive
residual distance and NP1 fails.  Thus any terminal-persistent near-null
sequence must approach the complete-word compatibility line

    nu_W=(theta_hat_0,0,...,0,-q).                          (NP6)

This proves the qualitative suffix dichotomy:

    A_obs,n->0 and A_term,n not->0
      => terminal root approaches the NEXT compatibility
         family, unless one of S/process, magnetic/gyro or
         accelerometer/BA actions stays positive.           (NP7)

### Connection to adjacent-superword kernel return

Let N(r) be the compact family of normalized word-dependent compatibility
vectors.  Define the distance of the terminal quotient image to the next
kernel family,

    d_next(x)=inf_(nu in N_next(r),alpha)
                 ||x-alpha nu||_(Pi^-1).                   (NP8)

The suffix argument gives qualitatively

    A_obs,n->0 => d_next(H_eff v_n)->0.                    (NP9)

Therefore the only possible O(1) terminal escape from relative detectability
is INTO the next word's scalar compatibility kernel.  This connects C_det
and the scalar return exactly as desired.

A quantitative combined inequality would be

    d_next(H_eff v)^2
      <= C_det,perp A_obs(v),                               (NP10)

and

    sup_(nu in N_W,nu+ in N_next)
      nu_+^T P_(nu,W) nu_+
      <= D_100(c,r).                                       (NP11)

Then NP10 controls the terminal component transverse to the next kernel and
NP11 controls the scalar component along it.  The recurring contraction no
longer needs C_det and kernel return as unrelated assumptions.

### Why NP10 is NOT yet proved source-uniformly

The qualitative implications above do not provide a finite uniform linear
rate.  Three quantitative moduli remain:
- the normalized four-S homogeneous-chain singular modulus;
- the residualized magnetic/gyro modulus;
- the multi-epoch accelerometer/BA distance to the moving compatibility
  family.

More importantly, the third modulus has a reachability boundary case that
cannot currently be excluded.  The existing exact-compatible MOVING analysis
shows that a nominal field-compatible force history can satisfy all
accelerometer kernel equations if the closed-loop AW/BA/S mean recursion can
return to that manifold.  Prediction leaves the zero-BA manifold by

    P_b fhat_next^-=(phi-1)P_b g,                           (NP12)

so the intervening corrections must supply

    P_b Delta a_hat=(1-phi)P_b g + S/BA terms.             (NP13)

At an accelerometer correction

    Delta a_hat=K_aw r_acc.                                (NP14)

Local exact reachability therefore depends on

    rank(P_b K_aw)=2                                       (NP15)

plus longitudinal/S/BA cycle closure and physical measurement admissibility.
The current proof has neither a source-uniform lower rank/singular-value
certificate for P_b K_aw nor an invariant forcing it singular.  No exact
shipping-compatible MOVING trajectory satisfying NP12--NP15 and strict
MAGNETIC SERVICE has been constructed either.

Hence a sequence may in principle approach an exactly compatible reachable
manifold while carrying O(1) terminal kernel amplitude.  That behavior is
precisely the next-kernel branch of NP7, so it does not invalidate the
dichotomy, but it prevents promotion of NP10 to a numerical source-uniform
constant until the adjacent-kernel return/reachability problem is solved.

### Result

The near-null persistence lemma is CLOSED qualitatively in the combined form:

    near-zero complete observation action
      => terminal persistence either vanishes transversely
         or converges into the next compatibility kernel.  (NP16)

It is NOT closed quantitatively.  The correct next theorem is a COMBINED
quotient+kernel inequality, not separate C_det and D bounds:

    ||Proj_(N_next)^perp H_eff v||_(Pi^-1)^2
      + (1/c_next)|Proj_(N_next) H_eff v|^2
      <= C_joint(c,r) [ v^T S_q v + (1/c)|Proj_(N_W)v|^2 ].
                                                               (NP17)

If C_joint is finite uniformly on retained 100-s pairs, it directly gives
the two-superword Riccati diameter and scalar invariance in one estimate.
The exact-compatible reachability manifold is then allowed as the kernel
component rather than needing to be excluded.


## 68. Two-word kernel-augmented backward reader: exact Gram inequality

Let N_W and N_+ be the normalized compatibility-kernel subspaces at the
current and next 100-s roots (dimension zero or one).  Let P_+, P_+^perp be
their Euclidean projectors.  Define the terminal metric

    Q_+ =
      P_+^perp Pi^-1 P_+^perp
      +(1/c_+) P_+.                                       (JR1)

Factor Q_+=C_+^T C_+ (Moore-Penrose square root on its support).  The rows of
C_+ are terminal readouts: Pi-whitened transverse rows plus one next-kernel
row scaled by c_+^-1/2.

Run ALL rows backward through the literal current 100-s word.  For operation
j use the exact common-source recursion:
- prediction x+=F x+U w:
      Y_j=Y_(j+1)F,     Z_j=Y_(j+1)U;
- correction with row H and noise factor V:
      choose reader L_j,
      Y_j=Y_(j+1)-L_j H,   Z_j=L_j V;
- reset:
      Y_j=Y_(j+1)G;
- PSD sync:
      identity mean pullback plus its actual source factor.

Stack every Z_j in ONE chronological source row Z(L).  Then for root e_0 and
fresh source vector s,

    C_+ e_N = Y_0(L)e_0 + Z(L)s.                           (JR2)

Minimize the common source action ||Z(L)||_F^2 over all correction-reader
blocks L.  This is one finite quadratic least-squares problem.  Denote its
minimum by B_+ and the associated root row by Y_+.

The exact joint terminal quadratic after optimal source reuse is

    e_0^T Y_+^T Y_+ e_0.                                  (JR3)

At the current root define the action metric

    Q_0 =
      S_q +(1/c)P_W,                                       (JR4)

where P_W projects onto the current compatibility kernel and S_q is the
complete quotient Schur information from the same word.

The desired combined inequality is exactly

    Y_+^T Y_+ <= C_joint Q_0.                              (JR5)

Its SHARP fixed-pair constant is

    C_joint(W,W+)
      =lambda_max[
        Q_0^dagger/2 Y_+^T Y_+ Q_0^dagger/2 ],             (JR6)

with infinity iff

    Null(Q_0) not subset Null(Y_+).                         (JR7)

Thus the proposed backward-reader construction succeeds algebraically: it
combines transverse terminal persistence and next-kernel precision in one
generalized Gram pair, with every Joseph/process/sync factor shared exactly.

### Joseph identities do not by themselves prove a uniform constant

JR5 is NOT an automatic consequence of covariance monotonicity.  Joseph
identities evaluate the action of a chosen reader and make JR6 exact, but
they do not imply the range inclusion JR7 uniformly over changing kernel
pairs.  The remaining condition is two-word finite-horizon detectability.

The qualitative suffix result NP16 gives precisely the limiting range
classification needed:

    Q_0,n-action ->0
       => dist(C_+ e_N,N_+)->0.                            (JR8)

Because Q_+ itself assigns finite precision along N_+, zero RHS action also
requires the current-kernel component to be controlled by the propagated
kernel precision.  Therefore the only possible violation of JR7 is a
current zero-cost compatibility direction that maps to a nonzero next-kernel
component while carrying no current scalar precision.  The +(1/c)P_W term
removes that possibility for finite c.  Hence for each fixed admissible pair
with finite c,c_+,

    Null(Q_0) subset Null(Y_+).                             (JR9)

So C_joint(W,W+)<infinity pointwise.

### Uniformity and rank-changing kernel pairs

Pointwise JR9 still does not supply

    sup_(W,W+) C_joint(W,W+) <infinity.                    (JR10)

As kernels rotate, appear or disappear, the smallest positive eigenvalue of
Q_0 can approach zero.  The numerator must vanish at the same rate.  NP16
gives convergence to the next kernel but no linear modulus.  Therefore
compactness plus pointwise range inclusion is insufficient.

The exact source-uniform target is the regularized two-word range estimate

    ||Y_+ v||^2
      <= C_joint(c,r)
         [v^T S_q v +(1/c)||P_W v||^2]                    (JR11)

for all retained adjacent superword pairs.  JR11 is simultaneously:
- quotient finite-horizon detectability;
- scalar current-to-next kernel return;
- continuity through kernel rank changes.

It is strictly weaker than any uniform information eigenvalue floor.

### Useful decomposition of JR11

Decompose v=v_perp+alpha nu_W.  The exact minimum-action construction gives

    Y_+ v =
      Y_+ v_perp + alpha Y_+ nu_W.                         (JR12)

For the quotient piece NP16 says the transverse next-kernel component is
charged by complete observation/process action.  For the kernel piece,

    ||Y_+ nu_W||^2

is exactly the next-kernel readout action of the current compatibility mode;
its bound relative to 1/c is the old scalar return, now measured in Q_+
rather than separately.

The cross term is not bounded independently.  Complete the 2x2 block square:

    [v_perp;alpha]^T
      [A  b; b^T d]
    [v_perp;alpha]
      <= C
    [v_perp;alpha]^T
      [S_q 0;0 1/c]
    [v_perp;alpha].                                       (JR13)

The sharp C is again JR6.  Thus no factor two or separate C_det+D budget is
needed.

### Remaining analytical obligation

The augmented backward reader therefore DOES combine the two obligations in
one shot algebraically, but it does not manufacture a source-uniform constant
from Joseph alone.  To close JR11 one must prove a LINEAR RATE version of the
near-null suffix dichotomy over the compact retained adjacent-word family.

Equivalently, exclude sequences with

    v_n^T Q_0,n v_n ->0,
    ||Y_+,n v_n||^2 / (v_n^T Q_0,n v_n) ->infinity.        (JR14)

The structural limit of any such sequence is an exact-compatible
current-to-next kernel trajectory.  Unlike earlier approaches, this need not
be excluded.  Linearize the literal two-word mean/covariance chronology
TRANSVERSE to that compatibility manifold.  If the transverse derivative has
full rank in the Q_+/Q_0 quotient, the implicit-function/closed-range theorem
gives a local finite JR11 modulus.  A finite cover of the compact
compatibility family then gives C_joint(c,r)<infinity.

Thus the next concrete calculation is the TRANSVERSE JACOBIAN of the
exact-compatible manifold, including P_b K_aw control rank, S/BA cycle
closure and magnetic service.  Full transverse rank, not existence/nonexistence
of the compatible orbit, is the final qualitative condition needed for a
uniform two-word constant.


## 69. Literal transverse Jacobian: Schur elimination and explicit persistence direction

Linearize the exact-compatible two-word boundary equations about a regular
strict-margin compatible execution.  At accepted accelerometer epoch k use
the innovation u_k in R^3 as the physical/measurement control and let

    E_k(m_k,u_<k)=0 in R^2                                (TJ1)

be the two transverse compatibility equations after projecting the nominal
force/BA relation onto b^perp.  The literal mean recursion, including
prediction and endogenous S feedback, is

    m_(k+1)=F_k(m_k,u_k).                                  (TJ2)

### Transverse block

At fixed pre-correction covariance/gain, the AW mean increment is

    delta a_hat = K_aw,k u_k.                              (TJ3)

Therefore the derivative of TJ1 with respect to the two transverse innovation
components is, up to the invertible frame/projector factors already present
in E_k,

    D_perp,k = P_b K_aw,k |_bperp.                         (TJ4)

The full two-word transverse Jacobian is block lower triangular in
chronological innovation controls because u_k first enters the current
correction and affects later equations only through TJ2.  Its diagonal blocks
are D_perp,k.  Hence

    rank D_perp,k=2 at every required epoch                (TJ5)

implies full row rank of the complete transverse compatibility Jacobian.

This is the exact implicit-function condition previously identified by the
reachability analysis.

### Longitudinal/S/BA Schur elimination

Choose two transverse components of each u_k as dependent controls.  Under
TJ5 the IFT gives

    u_(k,perp)=psi_k(m_k,u_(k,parallel)),                  (TJ6)

and substitution produces the reduced compatible recursion

    m_(k+1)=Fhat_k(m_k,ell_k),
    ell_k=u_(k,parallel).                                  (TJ7)

There is NO independent longitudinal S endpoint equation.  The S residual is
endogenous to m and is already applied inside Fhat.  There is likewise NO
separate BA endpoint closure: the estimator BA mean is part of m, whereas
the compatibility ratio q is the homogeneous error-kernel coordinate.  The
two-word construction is not required to be periodic in S or BA mean.

Thus after transverse Schur elimination the longitudinal variables ell_k are
free controls subject only to physical/gate realization.  The alleged
longitudinal/S/BA overdetermination vanishes.

### Magnetic transport

MAGNETIC SERVICE acts transversely to the allowed field-axis compatibility
line.  Along an exact compatible field-axis kernel its measurement loss is
zero while its information Gram may remain strictly above the service floor.
Perturbing the compatible innovations changes transported magnetic rows
continuously.  Therefore at a base execution with strict service surplus

    lambda_min(G_M) >= mu_M+delta_M,   delta_M>0,           (TJ8)

all sufficiently small IFT controls preserve MAGNETIC SERVICE.  The magnetic
block does not add an independent row that destroys TJ5.

### Consequence: full transverse rank supports, rather than excludes, persistence

The remaining two-word boundary conditions for a unit-persistent kernel are

    theta_1=F_0 theta_0,                                   (TJ9)
    q_1=Phi_b,0 R_ba q_0,                                  (TJ10)

together with the compatibility equations TJ1 throughout W0 and W1.
After TJ6 these are propagated kernel identities, not extra estimator
S/BA-periodicity equations.

Hence, IF there exists one reachable recurring strict-margin A21 base that
simultaneously:
1. lies on an exact compatibility line;
2. satisfies TJ5 on the needed correction epochs;
3. has strict MARINE/gate and magnetic-service margins,

then the implicit-function theorem constructs a local exact-compatible
two-word family.  Physical accelerometer samples are realized by

    f_phys,k=f_hat,k+b_hat_a,k+u_k,                         (TJ11)

and sufficiently small controls preserve finite-horizon acceleration, jerk,
velocity, displacement and strict service margins by continuity/Hermite
interpolation.  Such a base would yield an admissible persistent
current-to-next kernel pair and therefore the equality branch of the
two-word return.

This means the hoped-for conclusion

    "full transverse rank => finite transverse loss away from persistence"

has the WRONG sign near a compatible base.  Full rank makes compatibility
locally CONTROLLABLE.

### Explicit null/tangent direction of the augmented two-word Jacobian

On the exact-compatible manifold, differentiate the IFT family with respect
to any free longitudinal control parameter ell.  Let

    delta m_k = d m_k/d ell,
    delta u_(k,perp)=D psi_k delta(m_k,ell).                (TJ12)

By construction,

    D E_k [delta m_k,delta u_k]=0                          (TJ13)

at every accelerometer compatibility row.  Differentiate the propagated
kernel identities TJ9--TJ10 and magnetic-compatible field-axis relation.
The resulting nonzero tangent vector

    z_tan =
      (delta m_0, delta u_0,...,delta u_N,
       delta theta_0,delta q_0)                            (TJ14)

lies in the nullspace of the TRANSVERSE compatibility Jacobian while moving
along the exact-compatible two-word manifold.  This is the explicit null
direction requested.  It is not an instability mode by itself; it is a
tangent/reachability direction in the augmented physical+estimator control
space.

After quotienting the compatibility manifold, TJ14 is removed.  The normal
Jacobian is full rank exactly under TJ5 and the ordinary nonsingularity of
the recursive state chart.  Therefore the local closed-range estimate needed
for JR11 holds NEAR ANY SUCH REGULAR COMPATIBLE BASE.  The problematic places
for uniformity are instead:
- rank loss of P_b K_aw;
- gate/service/event-stratum boundaries;
- appearance/disappearance of the compatibility line;
- lack of a reachable compatible base.

### Current theorem consequence

There is NO explicit transverse null direction outside the compatibility
manifold produced by the literal algebra when TJ5 holds.  Conversely the
current proof does not certify TJ5 source-uniformly on recurring A21 roots.
At constructor/diagonal covariance the AW correction block is full rank, but
that is not a certified recurring compatible root.

Therefore the remaining global question is NOT symbolic Jacobian rank.  It is
REACHABLE-BASE / STRATIFIED-RANK coverage:

    every exact-compatible retained two-word orbit
      either has rank-two P_b K_aw and hence a regular
      compatibility manifold with finite normal C_joint,
      or lies in a rank-deficient stratum for which the
      normal two-word Gram JR6 must be bounded separately.  (TJ15)

If a rank-deficient compatible stratum has an extra normal null vector with
nonzero terminal Q_+ image, it is the explicit obstruction to uniform
contraction.  If every such extra null vector is also terminal-null or moves
into the next kernel, C_joint remains finite.

This is the precise remaining linear proof obligation.


## 70. Classification of rank-deficient transverse AW-gain strata

At an accepted accelerometer correction the literal gain is

    K = P H_a^T Omega^-1,                                  (RD1)

with Omega=H_a P H_a^T+R_eff positive definite.  The AW block is therefore

    K_aw = N_aw Omega^-1,                                  (RD2)

    N_aw =
      P_aw,theta J_att^T
      +P_aw,aw R_wb^T
      +P_aw,ba I
      +P_aw,bg J_bg^T,                                     (RD3)

where the BA term is present when its mean is active in the correction and
the gyro-bias lever-arm term is present when enabled; frozen BA uncertainty
may remain in Omega but not in the mean-gain numerator.

Let E_b be any 3x2 orthonormal basis for b^perp.  The transverse control block
used in TJ4 is

    D_b = E_b^T K_aw E_u,                                  (RD4)

where E_u selects the two innovation coordinates used as transverse controls.
Because Omega^-1 is invertible, it is cleaner to classify the coordinate-free
map

    G_b = E_b^T N_aw Omega^-1 : R^3 -> R^2.                (RD5)

Its rank equals rank(E_b^T N_aw).  Choice of two control coordinates E_u can
be made after this classification.  Hence innovation covariance conditioning
cannot create or remove a rank-deficient transverse OUTPUT direction.

### Rank-two regular stratum

    rank(E_b^T N_aw)=2.                                    (RD6)

Then there exists a 2-column innovation subspace E_u with det(D_b)!=0.
This is the regular IFT stratum of TJ: the exact-compatible manifold is
locally controllable and, after quotienting its tangents, the normal
two-word Jacobian is full rank.

### Rank-one stratum

There exists a unique unit transverse AW-output direction d in b^perp such
that

    d^T N_aw=0.                                             (RD7)

Equivalently,

    P_aw,theta^T d acted through J_att
    + P_aw,aw d acted through R_wb
    + P_aw,ba^T d
    + P_aw,bg^T d acted through J_bg

cancel exactly in the accelerometer covariance numerator.  The missing
direction is a LEFT null direction of the AW correction map: no
accelerometer innovation can instantaneously change d^T a_hat_w.

The corresponding compatibility defect after prediction is

    delta_d =
      d^T P_b[(phi-1)g + S/BA/transport terms].             (RD8)

If delta_d!=0 at a required return epoch, exact compatibility cannot be
restored there by the accelerometer correction.  Therefore the rank-one
stratum splits:

R1-X (incompatible rank loss):
    d^T required transverse return !=0.                     (RD9)

This stratum cannot contain an exact-compatible orbit and is irrelevant to
the compatibility-manifold uniformity problem; its nonzero defect supplies
normal observation/action.

R1-K (kernel-aligned rank loss):
    d^T required transverse return =0 at every deficient
    epoch and the later S/BA/magnetic chronology preserves
    that equality.                                         (RD10)

Then the missing control direction is tangent to an enlarged exact-compatible
stratum.  Its terminal classification is obtained by propagating the
homogeneous AW-output covector d through the literal later mean maps.  Let

    h_d = Q_+^(1/2) M_(+<-k) E_aw d.                        (RD11)

If h_d=0, the extra direction is TERMINAL-NULL.
If h_d lies entirely in the next-kernel row of Q_+, it is NEXT-KERNEL.
If P_+^perp h_d!=0, it is TERMINAL-PERSISTENT and is an explicit obstruction
to JR11 unless some later non-accelerometer action charges it.

But exact S-chain cancellation gives the needed classification for a PURE
zero-source AW direction: a nonzero homogeneous AW perturbation generates

    S(t)=a_0 psi_tau(t)+...                                (RD12)

and four distinct S=0 rows force a_0=0 when all accompanying v,p,S root
coefficients and fresh process action are zero.  Therefore a rank-one missing
AW correction direction cannot remain a zero-action terminal-persistent PURE
AW mode through a suffix containing four distinct S observations.

Consequently any R1-K terminal-persistent obstruction must be a COUPLED slow
mode whose AW component d is accompanied by attitude/BA (and possibly AG)
components that cancel the S/accelerometer action.  Such a coupled mode is
exactly part of the complete word-dependent compatibility kernel already
represented by P_+, provided its magnetic rows are also zero.  If magnetic
rows are nonzero, MAGNETIC SERVICE charges it.

Thus, under the already-proved complete-word nullspace classification,

    rank-one deficient + zero complete action
      => terminal-null OR next-kernel.                     (RD13)

There is no extra fixed-word terminal-persistent normal null direction on the
rank-one stratum.

### Rank-zero stratum

    E_b^T N_aw=0.                                          (RD14)

Both transverse AW-output directions are uncontrollable by the instantaneous
accelerometer correction.

Again split by the required compatibility return vector r_perp.

R0-X:
    r_perp !=0.                                             (RD15)

Then exact compatibility is impossible at that epoch unless an intervening
S/other correction supplies the full return before the compatibility row.
The literal chronology decides this before declaring the stratum compatible.
If no such prior supply exists, the stratum is excluded and carries positive
normal action.

R0-K:
    r_perp=0 after ALL prior operations at every deficient
    epoch.                                                  (RD16)

Both missing AW directions are tangent candidates.  However the complete
zero-action word still has at most a ONE-dimensional slow compatibility
kernel by the established S+magnetic+accelerometer classification.  Hence two
independent terminal-persistent normal directions cannot survive.  At most
one combination can join the next compatibility kernel; every independent
complement is either S/process charged, magnetic charged, or terminal-null.

Therefore

    rank-zero deficient + zero complete action
      => at most one next-kernel direction;
         all other missing transverse directions terminal-null
         or positive-action.                               (RD17)

### Terminal image theorem for all rank-deficient compatible strata

Combine RD13 and RD17.  For any FIXED admissible literal word satisfying the
existing complete-word nullspace theorem,

    Null(normal compatibility/action Jacobian)
       subset Null(Q_+^(1/2) M_terminal)
              + span(next compatibility kernel).           (RD18)

Thus no rank-deficient P_b K_aw stratum creates an additional fixed-word
terminal-persistent normal null direction.  Rank deficiency changes the local
parameterization of the compatibility manifold but not the complete-word
zero-action terminal classification.

This answers the fixed-word classification requested:
- rank 2: regular compatible manifold; tangent nulls only;
- rank 1: missing direction is either incompatible/positive-action, or
  terminal-null/next-kernel after complete chronology;
- rank 0: same, with at most one next-kernel combination because the complete
  slow kernel dimension is <=1.

### What remains source-uniformly open

RD18 is qualitative.  Near a rank-one/rank-zero boundary, the normal action
can vanish faster than the terminal transverse image even though the exact
boundary null is terminal-null/next-kernel.  Therefore source-uniform
C_joint still requires a quantitative rate across these strata.

The useful algebraic stratification variable is the smallest nonzero singular
value of

    E_b^T N_aw.                                             (RD19)

On regions where it is >=delta, IFT/closed-range gives a local modulus.
Near delta=0, use the exact deficient-direction split RD9--RD17 and derive a
second-order/next-operation modulus from S/process or magnetic action.  A
finite semialgebraic/interval cover of these strata would then certify
C_joint.

No assumption that P_b K_aw is uniformly full rank is required.


## 71. Near-rank-loss quantitative modulus: four-S + residualized magnetic block

Let d in b^perp be a unit left near-null direction of the literal AW
accelerometer numerator at epoch k:

    ||d^T N_aw,k|| <= delta.                               (NM1)

Because Omega_k^-1 is bounded on the retained class, the instantaneous
accelerometer authority in d is O(delta).  A useful lower bound cannot be
taken from the NEXT S event alone: the carried (v,p,S) root can cancel one
S residual exactly.  The first source-valid object is the four-S
Schur-completed block.

### Four-S residualized modulus

Let t_1<...<t_4 be the next four distinct accepted S epochs in the regular
suffix and let h_j=t_j-t_k.  With zero fresh LIN/AW process action, the
homogeneous scalar LIN response in direction d is

    s_j =
      S_0+p_0 h_j+v_0 h_j^2/2+a_0 psi_tau(h_j),            (NM2)

where a_0=d^T delta a_w,k is the missing AW component and

    psi_tau(h)=tau^3(h^2/(2 tau^2)-h/tau+1-exp(-h/tau)).   (NM3)

Stack

    V_S =
      [1 h_1 h_1^2/2 psi_tau(h_1)
       ...
       1 h_4 h_4^2/2 psi_tau(h_4)].                        (NM4)

Let Rbar_S be the FULL four-event residual covariance after transporting all
common fresh source factors and shorting every nuisance source except the
homogeneous root (S0,p0,v0,a0).  Rbar_S>0 on the retained class.

Partition V_S=[V_0,v_a], V_0 in R^(4x3).  Eliminating S0,p0,v0 gives the exact
scalar Schur information for the missing AW direction

    gamma_S =
      v_a^T Rbar_S^-1/2
        (I-P_(Rbar_S^-1/2 V_0))
      Rbar_S^-1/2 v_a.                                    (NM5)

Strict Chebyshev independence implies gamma_S>0 for every fixed finite tau
and four distinct epochs.  On a compact cadence/tau class with a positive
minimum separation between the selected four epochs,

    gamma_S >= gamma_S,* >0.                              (NM6)

NM6 is the correct quantitative S modulus.  It does not use the retracted
raw Vandermonde determinant floor; gamma_S,* must be enclosed from the
literal normalized covariance/cadence ranges.

If fresh LIN/AW process action is allowed, the joint Schur complement simply
adds its normalized action.  Therefore for the missing direction amplitude
a_0,

    A_S+proc >= gamma_S,* |a_0|^2.                         (NM7)

### Residualized magnetic/gyro modulus

Let z_ag be the AG component induced by the same near-null compatibility
direction after the accelerometer correction and transport it to a qualified
magnetic-service block.  Stack the actual magnetic rows and gyro/process
sources over one service interval.  After whitening by their FULL common
source covariance and projecting out the transported field-compatible
attitude/gyro line, define

    G_M,res =
      O_M^T Sigma_M^-1/2
        (I-P_M,nuis)
      Sigma_M^-1/2 O_M.                                   (NM8)

For the component z_perp transverse to the magnetic-compatible line,

    A_M+gyro >= z_perp^T G_M,res z_perp.                   (NM9)

The raw MAGNETIC SERVICE premise lower-bounds an unshorted service Gram.  It
does NOT automatically lower-bound G_M,res.  Hence the valid source-uniform
constant is

    gamma_M,* =
      inf_(retained qualified service blocks)
      lambda_min^+(G_M,res).                               (NM10)

A useful proof requires gamma_M,*>0.  This is exactly the residualized
magnetic/gyro modulus previously identified as open; no unprojected leverage
constant is substituted.

### Joint near-rank-loss block

Propagate the unit missing AW-output direction d from epoch k to the selected
four-S block and magnetic block using the literal homogeneous mean maps.
After eliminating current/next compatibility-kernel coordinates, write the
resulting terminal-normal image as

    h_+(d)=T_S d + T_M d + r_delta,                         (NM11)

where r_delta is the contribution of the small but nonzero accelerometer
authority.  On the retained gain/noise class

    ||r_delta||_(Q_+) <= C_K delta.                        (NM12)

Define the two normalized residual maps

    B_S d = sqrt(gamma_S,*) A_S d,
    B_M d = G_M,res^(1/2) A_M d.                           (NM13)

and stack

    B_N=[B_S;B_M].                                         (NM14)

The exact source-valid near-rank-loss modulus is the generalized singular
value

    beta_* =
      inf_(retained deficient strata, |d|=1)
      ||B_N d||^2 /
      ||P_+^perp T_+ d||_(Pi^-1)^2.                        (NM15)

with the convention beta=infinity when the terminal transverse image is zero.
If beta_*>0, then

    ||P_+^perp T_+ d||_(Pi^-1)^2
      <= beta_*^-1 (A_S+proc+A_M+gyro)
         + C_delta delta^2.                                (NM16)

This is precisely the requested combined normalized-action versus terminal
Q_+ bound near rank loss.

### Can beta_*>0 be proved from the CURRENT certificates?

Not yet.  NM6 requires a literal normalized four-S interval enclosure and
NM10 requires a residualized magnetic/gyro floor.  Neither numerical
source-uniform constant is presently certified.  More importantly, even
gamma_S,*>0 and gamma_M,*>0 separately do not guarantee beta_*>0 if the
terminal-normal direction can approach the joint null of the propagated S
and magnetic blocks while moving into the next compatibility kernel.  After
P_+^perp projection that kernel motion is harmless; the fixed-word
classification RD13/RD17 shows the remaining exact joint null has zero
terminal-normal image.  Therefore pointwise beta>0 holds on every fixed
nondegenerate stratum.  Uniform beta_*>0 still requires continuity through
rank/cadence strata.

### Finite-cover formulation

The retained parameter set splits into finitely many combinatorial event
types once accepted-correction/gate boundaries are treated as separate
closed strata.  On each interior stratum:
- tau and S cadence lie in compact intervals;
- selected S epochs have positive separation;
- covariance/source factors are bounded and positive;
- the maps in NM5--NM15 are continuous.

If interval arithmetic certifies on each stratum

    gamma_S >= gs_j>0,
    G_M,res >= gm_j P_M,perp,
    P_+^perp T_+^T T_+ P_+^perp <= t_j I,                 (NM17)

and excludes an extra joint null by a lower singular enclosure of B_N on the
terminal-active subspace, then

    beta_j>0,   C_joint,j <= max(1/beta_j,C_regular,j).    (NM18)

A finite maximum gives the desired source-uniform C_joint.

Thus the near-rank-loss calculation is now reduced to TWO concrete interval
certificates from literal chronology:
1. four-S Schur information gamma_S,*;
2. residualized magnetic/gyro information gamma_M,*,
followed by one 2--3 dimensional joint generalized singular-value enclosure.
No eventwise TV, sqrt(N), or uniform full-rank AW-gain assumption is needed.


## 72. Interval-certificate audit for gamma_S, gamma_M and beta_*: current source ranges are insufficient for a rigorous positive number

The requested interval proof arithmetic was attempted from the literal
shipping/source ranges.

Certified literal ranges currently available include

    tau in [0.02,12] s,                                    (IC1)
    dt  in [0.004,0.006] s,                                (IC2)
    T_S <=0.15 s in OU-III Live,                           (IC3)
    sigma_aw>=0.05 m/s^2,                                  (IC4)
    sigma_acc>=0.05 m/s^2,                                 (IC5)
    sigma_S>=0.075 m s,                                    (IC6)
    Sigma_aw<=16 I,  R_S<=10000 I,                         (IC7)
    MAGNETIC SERVICE T_M=1 s, mu_M=1.                      (IC8)

These are sufficient to make every fixed carried four-S/magnetic block
finite.  They are NOT sufficient, by themselves, to certify positive
source-uniform interval lower bounds for NM5/NM10.

### Four-S interval obstruction

The exact modulus is

    gamma_S =
      v_a^T Rbar_S^-1/2
        (I-P_(Rbar_S^-1/2 V_0))
      Rbar_S^-1/2 v_a.                                    (IC9)

A direct interval box over four event times satisfying only
0<t_(j+1)-t_j<=0.15 has inf gamma_S=0: distinct epochs may coalesce.  The
shipping lower cadence clamp at small tau is not a useful fixed physical
separation for a source-uniform unscaled determinant.

This is NOT a failure of the S-chain theorem.  Select four proof epochs from
the many accepted S events using separated target cells.  Since the maximum
gap is 0.15 s, every interval of length 0.15 s contains an accepted S event.
For example choose one event in each cell

    I1=[0,0.15],
    I2=[0.30,0.45],
    I3=[0.60,0.75],
    I4=[0.90,1.05],                                       (IC10)

relative to a regular suffix start after allowing the first service gap.
Then selected epochs obey

    t_(j+1)-t_j >=0.15 s                                  (IC11)

and lie inside a 1.05-s block.  This converts the open event-time set to a
compact separated box.  On IC1+IC10 strict Chebyshev independence implies

    gamma_S,geom :=
      min_(tau,t_j) dist(v_a,span(V_0))^2 >0.              (IC12)

However converting IC12 to the NORMALIZED gamma_S in IC9 also needs a
source-uniform UPPER bound on the full transported residual covariance
Rbar_S.  IC7 bounds the applied local R_S but does not bound the complete
four-event residual covariance after common process/root/source transport.
That upper covariance is exactly part of the still-open complete detectability
comparison.  Therefore a numerical positive lower interval for IC9 cannot be
certified from IC1--IC8 without circularity.

### Magnetic interval obstruction

MAGNETIC SERVICE gives an unshorted normalized information floor mu_M=1 over
every 1-s service block.  The required quantity is instead

    gamma_M =
      lambda_min^+[
        O_M^T Sigma_M^-1/2
        (I-P_M,nuis)
        Sigma_M^-1/2 O_M ].                               (IC13)

Projection can remove an arbitrarily large fraction of an unshorted Gram.
No theorem currently supplies a source-uniform angle between the magnetic
attitude rows and the transported gyro/nuisance range.  Thus IC8 does NOT
imply gamma_M>0 numerically.  A naive interval enclosure using only IC8 has
lower endpoint zero.

### Consequence for beta_*

Since both normalized component certificates currently have rigorous lower
endpoint zero,

    gamma_S in [0,+infinity),
    gamma_M in [0,+infinity),                              (IC14)

the joint generalized singular-value enclosure from source ranges alone is

    beta_* in [0,+infinity).                               (IC15)

Therefore no positive numerical beta_* can honestly be exported yet.  Any
positive value obtained from carried words would be a diagnostic promotion.

### Non-circular certificate design

The interval task can still be completed, but the variables must be the
LITERAL finite-word factor matrices, not coarse scalar source boxes.

For each closed combinatorial event stratum:
1. choose four separated S epochs by IC10;
2. propagate interval enclosures of the exact common source-factor matrix A_S
   and observation/root matrix O_S through that <=1.05-s block;
3. compute the Schur complement gamma_S directly by verified QR/LDL, without
   separately bounding Rbar_S;
4. over one qualified magnetic service block propagate the JOINT AG/gyro
   source matrix [O_M,A_M] and compute the residualized Gram by verified
   QR/Schur elimination;
5. propagate the terminal-active map and next-kernel projector on the same
   stratum;
6. solve the resulting 2--3 dimensional verified generalized eigenproblem for
   beta_j;
7. bisect any interval box whose lower beta bound contains zero, splitting on
   tau, event times, covariance-factor entries and kernel angle;
8. treat exact rank-changing faces with the analytic RD13/RD17 terminal-null/
   next-kernel classification rather than forcing a positive Euclidean
   singular value there.

This is rigorous computational proof arithmetic.  It does not promote
carried minima: the carried word is used only to choose a subdivision/order,
while every accepted box must be enclosed from literal source recurrences.

### What must be added to the proof infrastructure

The current repository does not yet export interval enclosures of the
complete local source-factor matrices needed in steps 2--5.  Existing
certificates export scalar covariance/noise bounds and carried matrices, not
source-uniform interval matrices for arbitrary event/covariance histories.
Thus the requested gamma_S/gamma_M/beta_* numerical certificates cannot be
completed honestly in this turn by algebra alone.

The next implementation-proof task is precise:
- add a literal interval factor propagator for a <=1.05-s four-S block and a
  1-s magnetic-service block;
- use outward-rounded interval arithmetic;
- verify every accepted box against shipping tau/dt/cadence/noise/gate
  ranges;
- emit gamma_S_lower, gamma_M_lower and beta_lower only after all boxes close.

Until that tool exists, theorem status must remain

    gamma_S_source_uniform_numeric = false,
    gamma_M_residualized_numeric   = false,
    beta_rank_loss_numeric         = false.                (IC16)
