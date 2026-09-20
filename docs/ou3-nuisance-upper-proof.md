# Recurring covariance upper comparison for the 15 nuisance coordinates

This comparison supplies the LIN/BA part of the covariance upper bound needed
for coercivity of `V=e'P^-1 e`. It does not supply the remaining attitude/gyro
bound, a full-state contraction factor, or a nonlinear retained region.
All covariances below are the literal real-arithmetic covariance recursion in
the default bounded SpectralMSE profile used by the existing LIN certificate.
The regular A21 segment excludes frame/relock reconfiguration, uses steps in
[.004,.006] s, and carries the actual construction/hold/release state.

## Proposition

After 17 s of regular A21 operation, at every operation boundary,

`P_nn <= 5 diag(b_v^2 I_3,b_p^2 I_3,b_S^2 I_3,156^2 I_3,(1/1600) I_3)`.

Here `n=(v,p,S,a_w,b_a)` and the rational constants below are reproduced by
`nuisance_upper_certificate.py`. This is a **principal upper comparison** with
all nuisance cross covariance allowed. It is not a full-state lower bound or
an assertion that nuisance coordinates are independent.

## 1. Bounded stable forcing, including the actual covariance sync

In the default profile, `Sigma_aw <= 16 I`, `.02 <= tau <= 12`, and applied
`R_S <= 10000 I`. The clamps, convex parameter smoothing, S_factor=1, unit
SpectralMSE cadence normalization, and horizontal S factors .72/.50 give
these inequalities. The band noise floor cannot defeat the sigma upper
clamp: with a unit-variance auxiliary white input, the time-varying lowpass
has L2 norm at most 1, the highpass at most 2, and its following lowpass at
most 2, by convexity at each recursion from zero state. Thus its tracked
variance gain is at most 4; the default .12 noise floor is at most .24.
This argument freezes realized coefficients; it makes no independence claim
about physical input and adaptive coefficients.

Let epsilon be the rational relative process-covariance defect from
`small_x_source_defect()`. In natural step coordinates that proof gives
`Q_source <= (1+epsilon) Q_ideal`; the closed-form branch has zero defect.
The source small-x transition polynomials obey
`0 <= Psi_j(h) <= h^j/j!`, j=1,2,3, while the AW multiplier is `phi=exp(-h/tau)`.
The neutral 3x3 transition is the exact polynomial integrator transition.

The pending sync adds only `Delta=Pi_+(Sigma_target-P_aw)` to AW. Since
`Sigma_target <= 16 I` and `P_aw >= 0`, spectral calculus gives
`0 <= Delta <= 16 I`. There is at most one pending addition per prediction.
This is not a claim that the matrix positive-part map is operator monotone.
Put

`B=(1+epsilon)16+16(1+12/(2*.004)) < 156^2`.

The scalar recursion

`phi^2 B+(1+epsilon)16(1-phi^2)+16 <= B`

follows from `1/(1-exp(-x)) <= 1+1/x` and `x=2h/tau >= 2*.004/12`.
Initial AW covariance is 2.2^2 I; pre-live sync sets it to Sigma, and regular
optimal corrections only decrease its marginal. Attitude resets leave it
unchanged. Hence `P_aw <= B I` is inherited from construction.

For BA, the initial variance is .004^2. Holding leaves that decoupled block
unchanged; the single release floor is the same initial variance. In A21,
`phi_b=exp(-h/5000)` and
`Q_b=(1/1600)(1-phi_b^2) I`, because the driving density is 2.5e-7.
Prediction preserves `P_ba <= (1/1600) I`; corrections decrease the marginal,
attitude resets leave it unchanged, and mean projection does not change P.
Repeated external hold/release commands are not covered by this tail.

## 2. Exact cancellation of the unknown neutral root

Fix the realized coefficients, measurement schedule and **actual** sync
increments. Form a comparison covariance that keeps S corrections and omits
acc/mag corrections. Initialize it with the actual LIN marginal at the first
chosen observation. Its covariance dominates the actual LIN covariance:
LIN prediction is autonomous, S correction is the monotone Riccati map of
that marginal, other optimal corrections decrease the marginal, and attitude
reset is identity on LIN. Add the same realized Delta to both recursions;
do not recompute Delta from the comparison covariance.

This comparison admits an auxiliary zero-mean Gaussian linear model with
independent fresh process increments and measurement noise, used solely to
prove the deterministic matrix inequality. Preserve all within-step process
cross covariance and the whole initial LIN covariance. Its unconditioned AW
variance stays below B even if all measurements are ignored.

At target time T choose the latest applied S observation t2<=T, then the
latest t1<=t2-8, then the latest t0<=t1-8. The progress-preserving source
scheduler gives every applied S gap <=.156 s in this regular exact-arithmetic
profile. Thus `a=t1-t0,b=t2-t1 in [8,8.156]`, `d=T-t2 in [0,.156]`, and
`T-t0<=16.468<17`. The 17 s entry delay ensures all three observations exist
in the regular segment. Include t0's observation after initializing at its
pre-correction covariance.

For a noise-free neutral trajectory, S(t) is quadratic. Its interpolation
rows at t2, in observation order (t0,t1,t2), are

```
v:  2/[a(a+b)]       -2/(ab)       2/[b(a+b)]
p:  b/[a(a+b)]       -(a+b)/(ab)   (a+2b)/[b(a+b)]
S:  0                0             1
```

Propagate these rows to T with the neutral transition:
`L_v=v`, `L_p=p+d v`, `L_S=S+d p+d^2 v/2`.
Then **exactly** `L O=F_N(a+b+d)`, where O contains the three rows
`[t_i^2/2,t_i,1]` relative to t0. Thus this trial estimate cancels the entire
unknown initial (v,p,S), irrespective of its magnitude and cross covariance.
It does not restart either the physical history or the actual estimator.
Uniform row-l1 bounds, with g=.156, are

`W_v=4/8^2`, `W_p0=4(8+g)/8^2`,
`W_p=W_p0+g W_v`, `W_S=1+g W_p0+g^2 W_v/2`.

## 3. Trial error bound

Write the auxiliary neutral state as `z(t)=F_N(t-t0)z(t0)+w(t)`.
All neutral transition entries and AW injection coefficients are nonnegative.
Replacing each injection by `h^j/j!` and summing with the exact neutral
semigroup gives coefficient sum at most `17^j/j!`. Minkowski therefore bounds
the AW-driven standard deviation in any fixed spatial unit direction by
`sqrt(B)17^j/j!`. This permits arbitrary temporal correlations in AW.

For the separate fresh neutral process increments, the ideal driving density
is at most `2*16/.02=1600`. After propagation to a target at lag at most 17,
the j-th neutral impulse due to noise within one step of age s<=h is at most
`s 17^(j-1)/(j-1)!`. Indeed the transported exact OU convolution is bounded
by the integral of `(lag+u)^(j-1)/(j-1)!` over 0<=u<=s.
Independence across fresh increments, the source relative Q upper comparison,
and `sum h^3 <= .006^2 sum h` give the standard-deviation bound

`C 17^(j-1)/(j-1)!`, `C^2=(1+epsilon)1600*.006^2*17/3 < 1`.

Within-step AW/neutral correlation is not discarded: add the AW and fresh
neutral contributions by Minkowski, not by adding their variances.
Consequently one may take

`E_v=156*17+1`, `E_p=156*17^2/2+17`,
`E_S=156*17^3/6+17^2/2`.

The trial neutral error is `w_j(T)-sum_i L_ji(w_S(t_i)+noise_i)`.
The S noise standard deviation is at most 100, giving

`b_j=E_j+W_j(E_S+100)`, j=v,p,S.

The optimal S-only estimator has no larger error covariance than this causal
linear trial estimator. Take AW trial estimate zero. Each LIN marginal is
therefore bounded by its stated scalar squared standard deviation times I_3.
Together with the BA marginal, block Cauchy--Schwarz gives, for arbitrary
five-block vector x,

`x'P_nn x <= (sum_j b_j ||x_j||)^2 <= 5 sum_j b_j^2 ||x_j||^2`,

where `b_aw=156,b_ba=1/40`. This proves the proposition, including every
cross block, interleaved optimal corrections and every operation prefix
following the 17 s entry delay. It establishes a finite nuisance upper bound;
its conservative size is not promoted to a useful full-state decay rate.
