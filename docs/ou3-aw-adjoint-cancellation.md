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
