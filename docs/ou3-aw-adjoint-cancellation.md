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
