# Held-BA LIN BIBO reduction and the displacement boundary row

The displacement completion introduces no new reader architecture.  Its only
root term is q^T v0, one coordinate of the existing held-H18 LIN state
z=(v,p,S,a_w).

For the literal homogeneous comparison use the covariance storage
V_L=e_L^T P_LL^{-1}e_L only as a shorthand for the corresponding full-state
variational restriction; do not discard cross covariance in the parent
certificate. Prediction with positive process covariance is nonexpansive in
the carried covariance metric. Every covariance-matched S or accelerometer
correction is also nonexpansive by the exact Joseph measurement-energy
identity. PSD synchronization cannot increase fixed-error inverse-covariance
storage. Thus accelerometer updates need NOT be assigned an independently
bounded Euclidean gain or treated as damping.

The zero-action kernel is structural.  For the neutral chain
v'=0, p'=v, S'=p, three distinct S observations have rows
[t^2/2,t,1]. Their determinant is

    (t1-t0)(t2-t0)(t2-t1)/2 != 0.

Hence v0=p0=S0=0.  The remaining homogeneous integrated-OU a_w root contributes
the nonpolynomial exponential extension to S.  A fourth distinct S row kills
that root by the same extended-Chebyshev/four-S argument already used in the
corrected-word proof, uniformly for tau in the compact interval [.02,12].
Therefore no nonzero 12-state homogeneous LIN root can have zero S action over
a four-event separated word. Accelerometer action can only add nonnegative
measurement loss in the covariance metric.

This proves STRUCTURAL DETECTABILITY, not yet a source-uniform contraction
number. To pass from pointwise strict loss to a common rho_LIN<1 by compactness,
the family of literal held-H18 words must itself be compact independently of
how long magnetic refinement delays release. Required premises are:

* a source-uniform upper bound on P_LL and the relevant full covariance/cross
  covariance during held H18;
* compact reachable tau,sigma_aw,R_S,T_S and scheduler phase on that same
  all-time held history;
* bounded actual correction coefficients/gains with innovation covariance
  uniformly positive.

The current repository has a fresh-process covariance LOWER bound and finite
carried 17-s diagnostics. Neither supplies the first item. Using P>=L as an
upper/action bound would reverse the covariance order. Thus the requested
uniform BIBO theorem cannot honestly be promoted yet.

Conditional on those compactness premises, continuity plus the zero-kernel
result gives a finite separated four-S word with

    V_L,end <= rho_L V_L,root,   rho_L<1.

For affine bounded input d (AG/BG, held BA, physical/sensor forcing), variation
of constants over that fixed word gives

    sqrt(V_L,end) <= sqrt(rho_L) sqrt(V_L,root) + c_L ||d||_word,

and iteration gives the standard BIBO radius c_L dmax/(1-sqrt(rho_L)).
The displacement boundary row then follows immediately from
|q^T v0| <= ||E_v^T P_LL^(1/2)|| sqrt(V_L,root), with the same all-time
covariance upper bound. It is therefore one consequence of H18 BIBO, not a
separate proof obstruction.

This reduction closes the structural/nullspace question and identifies the
single remaining H18 premise: ALL-TIME HELD-H18 COVARIANCE/COEFFICIENT
COMPACTNESS. The next calculation should attack that covariance upper bound
from recurring S corrections plus the stable OU process, retaining
accelerometer corrections as covariance-decreasing Joseph operations and
retaining full cross covariance. If such an upper bound fails, a uniform
release-time theorem is required instead.
