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

    P^+=P-P h^T(hPh^T+R)^-1 hP                             (TP2)

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
separate posterior covariances: C_H(P1+P2) != C_H(P1)+C_H(P2).

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
