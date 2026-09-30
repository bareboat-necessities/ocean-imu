# OU-III AW mean: analytical variation and marginal passivity

Scope: the default bounded, real-arithmetic, regular A21 profile, with the
actual carried construction, covariance, tuner, clocks and magnetic history.
No estimator behavior is changed. This note supplies analytical identities
and an explicit but unusably large uniform variation bound. It does NOT
prove the small nominal-mean bound, capture, O1/O2, or the stability theorem.
Its only role in the existing proof path is to supply the nominal-mean
premise of Corollary A*, and ultimately the existing finite-error inequality.

The earlier four-scalar lift, the omission of direct magnetic mean forcing,
and the claim that pairing events automatically removes raw second variation
are not valid and are replaced below. The beta diagnostic on this branch is
not a certificate for this corrected lift; its output must not be promoted.

## 1. Literal mean lift: twelve coordinates, not four independent scalars

Write z=(a_w,v,p,S), each entry a world-frame three-vector. E_a and E_S
select its acceleration and integral-displacement blocks. Let L select these
coordinates, in this order, from the full shipping state. At each accepted
correction use the actual full-state gain and set K_L=L K. Neither the gain
nor the covariance is recomputed using a reduced filter.

For a prediction with applied h and tau, put phi=exp(-h/tau) and

    c1 = tau (1-phi),
    c2 = tau h - tau^2 (1-phi),
    c3 = tau h^2/2 - tau^2 h + tau^3 (1-phi).

The ideal analytic transition is the following block matrix, tensored with I3:

    F = [[phi, 0, 0, 0],
         [c1,  1, 0, 0],
         [c2,  h, 1, 0],
         [c3, h^2/2, h, 1]].

For literal source arithmetic use IntegratedOUChain's coefficients, including
its small-x polynomial branch, rather than silently substituting ideal c2/c3.
The source bounds 0<=cj<=h^j/j! used below cover that branch.

The S correction occurs inside time_update, after prediction and application
of pending AW covariance inflation. Its mean map is exactly

    z+ = A_S z-,       A_S = I-K_L,S E_S.

For a lever-arm-disabled accelerometer row let R and Rhat be true and nominal
world-to-body rotations, and let g0 be the filter gravity reference. Then

    r_acc = Rhat(a_phys-E_a z) + eta,
    eta = (R-Rhat)(a_phys-g0) - R(g_phys-g0)
          + b_phys-bhat_temp + sensor_error.

Consequently

    z+ = A_a z- + B_a a_phys + K_L,a eta,
    A_a = I-K_L,a Rhat E_a,
    B_a = K_L,a Rhat.

Equivalently use B_a=K_L,a R and move the corresponding attitude/acceleration
term out of eta. One convention must be retained throughout the calculation.
The expression above uses the Rhat convention.

Magnetic measurements DO change the LIN mean through cross covariance:

    z+ = z- + K_L,m r_m.

They are not merely modifications of subsequent gains. Quaternion reset and
pending covariance synchronization do not themselves rotate or add to these
world-frame LIN means in the qualified regime, but they change subsequent
coefficients. Magnetic residuals, the carried bias estimate, and eta are not
independent physical controls.

Freeze coefficients only on one realized history. Each operation has the form

    z_i = A_i z_(i-1) + B_i a_i + d_i.

Here d_i retains direct magnetic forcing, the accelerometer eta contribution,
and any explicitly charged implementation defect. Predictions and S updates
have B_i=d_i=0 in exact arithmetic. A rejected correction is the identity.
Select the output at the literal pre-accelerometer epochs, after any S update.

For a fixed unit world vector u perpendicular to the reference field, let
c_i=alpha_i E_a' u at output epochs, and c_i=0 otherwise, with alpha_i>=0
and sum alpha_i=1. The backward recursion

    lambda_N = c_N,
    lambda_(i-1) = c_(i-1) + A_i' lambda_i

gives the exact identity

    u' mu_W = lambda_0' z_0
               + sum_i beta_i a_i + sum_i lambda_i' d_i,
    beta_i = lambda_i' B_i.                                  (L)

beta_i is a row with THREE physical-acceleration input coordinates. Projecting
onto u does not make the three-axis recursion a closed four-scalar filter.
The identity is exact on the same history; it is not the derivative of an
independently rerun adaptive filter. A covariance bound does not bound z_0,
and the root term in (L) cannot be discarded at an arbitrary word boundary.

## 2. Exact first Abel bound, including sampling

Index the accelerometer sample grid by j=1,...,n, with h_j=t_j-t_(j-1).
Keep zero beta rows for rejected/missing corrections on this grid. Define

    w_j = beta_j/h_j,
    D1 = ||w_1|| + ||w_n|| + sum_(j=1)^(n-1) ||w_(j+1)-w_j||.

All row norms are Euclidean. With physical v'=a, ||v||<=V and ||a'||<=J,

    v_j-v_(j-1) = h_j a(t_j) + q_j,
    ||q_j|| <= J h_j^2/2.

Thus summation by parts gives exactly

    sum_j beta_j a(t_j)
      = w_n v_n-w_1 v_0
        - sum_(j=1)^(n-1) (w_(j+1)-w_j) v_j - sum_j w_j q_j.

In particular

    |sum_j beta_j a(t_j)|
       <= V D1 + (J/2) sum_j h_j ||beta_j||.                  (A1)

The last term cannot be dropped merely because acceleration has zero mean.
Different sampling/trapezoidal conventions require their own quadrature
identity rather than reusing the right-endpoint constant without checking.

## 3. Exact second Abel bound, including endpoints

Let d_j=(w_(j+1)-w_j)/h_j, j=1,...,n-1, and define

    D2 = ||d_1|| + ||d_(n-1)||
         + sum_(j=1)^(n-2) ||d_(j+1)-d_j||.

The physical relation p'=v gives

    v_j = (p_j-p_(j-1))/h_j + e_j,
    ||e_j|| <= A h_j/2,          A = sup ||a||.

Applying summation by parts to the interior sum in (A1), while keeping its
velocity endpoints, gives

    |sum_j beta_j a(t_j)|
      <= V (||w_1||+||w_n||) + P D2
         + (A/2) sum_(j=1)^(n-1) h_j ||w_(j+1)-w_j||
         + (J/2) sum_j h_j ||beta_j||,                       (A2)

where P=sup ||p||. The previous raw P*D2 statistic omitted both velocity
endpoints and this acceleration quadrature term; it was not the complete
second-Abel bound.

For the cited envelopes V=5.5, P=8.1, A=8.8 and J=100, even ignoring every
root, defect and sampling contribution, the 1.96133 target requires

    D1 < 0.3566054546 1/s,
    D2 < 0.2421395062 1/s^2.

These are necessary gates for these upper-bound certificates, not necessary
conditions for actual filter stability.

## 4. A fully explicit uniform bound exists, but does not close the target

This section proves finiteness rather than claiming a small constant. It uses
an outer bound of the reachable set and therefore is not evidence that a bad
history is reachable. It is included to distinguish an actual uniform bound
from a conditional formula containing an unevaluated supremum.

After 17 s of regular A21, the individual nuisance marginal comparisons in
ou3-nuisance-upper-proof.md give

    b_v = 340861/32,
    b_p = 71236727/800,
    b_S = 212911701659/800000.

Lemma B of ou3-world-frame-rows.md sharpens P_aw to (1+epsilon)16 I, with
0<=epsilon<1. Four-block Cauchy--Schwarz therefore gives P_LL<=D^2 for

    D = diag(12 I3, 21304 I3, 178092 I3, 532280 I3)

in (a_w,v,p,S) order. The innovation floors in the qualified profile are
R_acc>=0.05^2 I and R_S>=0.075^2 I. The optimal gain on the LIN rows obeys

    K_L S_innov K_L' = P_LL^- - P_LL^+ <= P_LL^- <= D^2,

so, if R>=r^2 I,

    ||D^-1 K_L|| <= 1/r.

Hence the scaled literal mean maps satisfy

    ||D^-1 A_a D|| <= 1+12/.05 = 241,
    ||D^-1 A_S D|| <= 1+532280/.075 < 7100000,
    ||D^-1 B_a|| <= 20.

For a prediction, 0<=cj<=h^j/j! and h<=.006 show

    ||D^-1 F D|| <= exp(ell h),
    ell = 12/21304 + 21304/178092 + 178092/532280
          + .006/2 (12/178092 + 21304/532280)
          + .006^2/6 (12/532280) < 1.

This follows by bounding the off-diagonal operator norm by the sum of its
nonnegative block coefficients; the diagonal has norm at most one.

Let N_a,N_S count applied acc and S corrections on the word, and H be the
sum of its prediction durations. Equation (L), output normalization, and
submultiplicativity give the explicit uniform coefficient-mass bound

    sum_j ||beta_j|| <= B0,
    B0 = 240 N_a exp(H) 241^(N_a) 7100000^(N_S).              (U0)

The overcount of the correction at the input event is conservative. The
regular scheduler bounds both counts by a fixed IMU-cell count on a fixed
word, with at most one acc correction and one S correction per cell. For a
16 s word rooted and ending at the pre-acc epochs, h>=.004 gives at most
4001 sample epochs; using 4002 and H<=16.006 is a conservative endpoint pad.
This bound is inherited from the existing recurring covariance theorem, not
from a restarted or freely chosen root covariance.

The triangle inequality now proves

    D1 <= 2 B0/h_min,
    D2 <= 4 B0/h_min^2,             h_min=.004.               (U1/U2)

These are uniform analytical bounds on the corrected literal LIN beta rows.
They are vastly too large to establish (A1) or (A2) below 1.96133. As already
recognized for the coarse nuisance approach, improving arithmetic precision
cannot fix this loss of correlation. No covariance-box closure is promoted.

## 5. An S update has exact restoring action without any gain-sign premise

Let X=P_LL^->0 and H_S=E_S. The actual LIN S gain and marginal covariance are

    K = X H_S' (H_S X H_S'+R_S)^-1,
    X+ = X-K (H_S X H_S'+R_S) K',
    A = I-K H_S.

Direct substitution, or the matrix inversion lemma, yields

    A' (X+)^-1 A = X^-1 - H_S' (H_S X H_S'+R_S)^-1 H_S.       (S)

Thus the homogeneous LIN storage decreases by exactly the measured S action.
Prediction also satisfies F'(F X F'+Q+Delta)^-1 F<=X^-1 when Q+Delta>=0.
There is no need to prove that each entry of P_aS is positive. Identity (S)
already holds for every positive definite X, including the reachable ones.
It does not imply componentwise monotonicity of a_w or small raw beta variation.

For error relative to the physical vessel, the S pseudo-measurement has
measurement mismatch nu=-S_phys. In a linear Kalman correction with residual
H e+nu the exact information balance is

    V+ - V- = nu' R^-1 nu - (H e+nu)' S_innov^-1 (H e+nu).

The physical S contribution is therefore a supply term, not automatically
zero or negative.

## 6. The missing correlation term can be written exactly

Partition the full pre-correction covariance and row as

    P = [[X,C],[C',N]],       H = [H_L,H_n],
    D_c = H_n C' X^-1,
    H_tilde = H_L+D_c,
    R_eff = R + H_n (N-C'X^-1 C) H_n'.

Schur complementation gives R_eff>0 and the exact marginal Kalman gain

    K_L = X H_tilde' (H_tilde X H_tilde'+R_eff)^-1.

But the literal nominal LIN map uses

    A_L = I-K_L H_L = A_tilde+K_L D_c,
    A_tilde = I-K_L H_tilde.

For any adjoint vector lambda, put t=K_L'lambda and
lambda_tilde=A_tilde'lambda. The marginal Joseph identity implies exactly

    ||A_L'lambda||_X^2 + ||t||_(R_eff)^2
      = ||lambda||_(X+)^2
        + 2 t' D_c X lambda_tilde + ||D_c' t||_X^2.          (C)

At an S update, H_n=0 and D_c=0: (C) is loss only. At an accelerometer
correction, H_n includes attitude and active BA and D_c need not vanish.
At a magnetic update H_L=0, but the LIN mean still receives K_L,m r_m and
its marginal covariance changes. Treating that update as only a future-gain
change misses an actual term.

Equation (C) identifies the exact correlation supply that a sharp complete-word
passivity argument must bound or cancel. A proof for the full 21-state
homogeneous error cannot be applied to the 12-state mean lift by deleting
H_n: that would silently discard D_c. No sign or reachability counterexample
for D_c is claimed here.

## 7. What the coupled tuning law supplies, and what it does not

For the unclamped TARGET law in SeaStateFusionFilter_OU_III.h,

    r_S = C_J (2 r_a)^(1/14) (sigma/c_sigma)^(6/7)
                   tau^(24/7)/sqrt(T_S).

A useful dimensionless identity retains all the parameters jointly:

    r_S^2 T_S/(sigma^2 tau^7)
       = C_J^2 c_sigma^(-12/7)
         (2 r_a/(sigma^2 tau))^(1/7).                       (T)

When T_S is proportional to tau, its logarithmic differential is

    d log r_S = (6/7) d log sigma + (41/14) d log tau.

Here sigma/c_sigma is the tuner's wave-amplitude estimate, not a proved
identity with true physical wave RMS. The source applies different smoothing
coefficients to tau/sigma and r_S, then commits candidates at the next sample
boundary. Floors, clamps and actual S cadence also matter. Thus (T) is NOT
an instantaneous equality of all four applied values. A joint proof must
carry the tuner states and this target-to-applied discrepancy explicitly.

Even exact target equality (T) does not bound D_c in (C), which also depends
on actual attitude/BA cross covariance and the measurement row. Therefore
substituting the tuning law into a scalar gain formula is not yet a uniform
bound on the literal beta variation.

## 8. Genuine jumps cannot be cancelled inside raw total variation

On a uniform grid suppose w has an isolated interior jump J_e and is constant
on both neighboring sides. Then Delta w/h has one nonzero entry J_e/h, and
its augmented variation is exactly 2||J_e||/h. In beta units the same term is
2||Delta beta||/h^2. This is a genuine contribution to the defined norm.
Pairing an event's two sides cannot remove it while claiming to bound that
same raw norm. An alternative signed supply can retain cancellation, but is
a different certificate.

A small observed signed AW mean, such as the older 0.371 m/s^2 tracking-error
statistic, neither estimates D1/D2 nor bounds them. That statistic also is
not the isolated physical-input term in (L): roots and other sources remain.

## 9. Exact cycle/jerk split: a different, cancellation-preserving certificate

The following algebra avoids differentiating rapid cadence modulation. It is
not a claim that its remaining uniform constants are already known.
Partition the acceleration grid into contiguous blocks C, for example using
actual scheduler events, and define

    wbar_C = sum_(j in C) beta_j / sum_(j in C) h_j,
    beta_low,j = h_j wbar_C,
    gamma_j = beta_j-beta_low,j.

Each block has sum gamma_j=0. For C=[a,b], set
G_j=sum_(i=a)^j gamma_i. Summation by parts gives exactly

    sum_(j=a)^b gamma_j a_j
      = -sum_(j=a)^(b-1) G_j (a_(j+1)-a_j).

Consequently the actual jerk bound yields

    |sum_j beta_j a_j|
      <= V D1(beta_low/h)
         + (J/2) sum_j h_j ||beta_low,j||
         + J sum_C sum_(j=a_C)^(b_C-1) h_(j+1) ||G_j||.       (CJ)

This uses bounded velocity on the slow component and bounded jerk on the
zero-sum cadence component. Unlike raw second variation, the cadence
component is charged through a primitive, not a derivative. It is valid for
nonuniform steps, nonperiodic schedules and history-dependent coefficients.
A uniform theorem would need to bound these two structured coefficient
functionals on the same reachable Riccati/tuner history, together with the
root and all sources in (L). Finite trace maxima cannot replace that proof.

## 10. Threshold qualification and present result

Corollary A* needs

    ||mu_W x b|| < ||g0 x b|| = g sigma_w.

The number 1.96133 uses g=9.80665 and sigma_w>=1/5. The newer physical
inclination allowance |I|<=80 degrees only gives sigma_w>=cos(80 degrees),
and hence 1.7029069 m/s^2 at that gravity. A bound below 1.96133 alone would
not exclude a 1.70-scale compatibility state over that larger domain. Field
defects must also be charged without double-counting contributions already
included in mu. Neither physical assumption is silently changed here.

Proved analytically in this note: the full-vector literal lift (L), complete
first/second-Abel inequalities (A1/A2), finite source-uniform outer bounds
(U0/U1/U2), sign-free marginal S restoring identity (S), exact correlation
supply (C), joint target normalization (T), and cycle/jerk identity (CJ).

Not proved: a useful source-uniform nominal-mean bound below g sigma_w, a
uniform small bound on D1/D2, a construction-reachable pathological execution,
or any end-to-end stability conclusion.

Failure classification: quantitative relaxation failure and correction of
an invalid reduced-model argument, not a physical counterexample. The prior
raw-TV tactic is not repeated with more samples, narrower arithmetic or a
silently different norm. The controlling unclosed quantity is the joint
word-level correlation/physical-supply bound, with actual tuner lag, magnetic
mean forcing and the carried root retained.
