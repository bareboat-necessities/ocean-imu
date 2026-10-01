# Linked finite-error supply: exact identity, feasibility, and scope correction

Controlling continuation of PR #639, audited at `6c96893640abbbc49caeb68a9fe27a27680850b3`.
The shipping sources are also unchanged on main `17df0e8d64ca07dd6d35ee9a8fea72a85e3eb0d7`.
The original algebra is retained. The current physical IMU domain is the SLOW + FAST
model below; estimator, tuner, scheduler and capture gates are unchanged.
The path remains construction -> capture -> magnetic H18 -> refinement/release ->
A21 -> regional practical stability. This note supplies the linked finite-error
inequality in that path; it does not claim its source-uniform premises.

## 1. Correct the prerequisites before using them

The earlier PR continuation incorrectly made three implications.

* **Homogeneous action is not base innovation energy.** For a covariance-matched
  correction, zero loss implies H e_h=0 for the homogeneous comparison. It does
  not imply r_base=0. See `ou3-corrected-word-proof.md`, ZF-5--ZF-8, and
  `ou3-field-axis-regional-exclusion.md`, section 1. Consequently it does not
  imply zero nominal correction increments, zero nominal S, or zero nominal AW.
  The four-zero Rolle argument is valid for an unforced exact-integral chain;
  it is not applicable to the nominal mean of an arbitrary zero-homogeneous-action
  base history. In addition the literal small-x phi_pa/phi_Sa polynomials are
  not exactly the continuous exponential chain. The claimed outer invariant-set
  closure based on this argument is withdrawn, not silently used below.
* **Live/H18 handoff is not A21 release.** The pre-Live LIN mean can be zero.
  H18 thereafter predicts and corrects it before BA activation. The new carried
  replay has Live at 31.84 s, A21 activation at 120.08 s, and a release LIN mean
  norm 0.711963... . Thus the zero-mean construction seed is not an A21 release
  box. A finite H18 duration on each execution is not a uniform duration bound.
* **A lower covariance bound does not give compactness of a large V sublevel.**
  Already in one dimension P=n^2, e=n has P>=1 and V=e^2/P=1 while both grow
  without bound. A compact joint physical/filter history class needs additional
  coercivity/upper bounds, closed event strata and retention. It cannot be
  obtained from P>=P_min alone. The finite-horizon image of a compact class,
  when that input class and continuity are actually proved, remains a useful
  conditional construction, not an invariant outer set.

The source-audited AW marginal ceiling, conditional radius-local exclusion at
r_FA=.15, correct BA marginal elimination, exact operation energy identities
and same-history requirement remain. The new outer exclusion and homogeneous
entry statements from the invalid implications above are not retained.

## 2. Current physical domain: unchanged MARINE, two-timescale IMU

Use the existing MARINE MOTION complete-episode windows and actual MAGNETIC
SERVICE. Do not add the former EXCITED_MOVING research proposal as a premise.
For both sensors use `e_i=b_i,s+b_i,f` with slow amplitude/rate and the
all-placed-window fast contract in `ou3-imu-two-timescale.md`, SF1--SF5.
Both fast horizon/cap pairs and assembled slow-budget qualification are OPEN.

At every operation, the BA/BG physical coordinates in e are SLOW bias errors.
The same slow history generates the bias-prediction mismatch, and the fast
history enters the actual calibrated sensor rows. Compose those rows before
bounding them. Let E_i(W|past) be the restriction of the ONE carried admissible
slow/fast history, including windows crossing W's boundary. The needed bound is

    chi_* = sup { chi_gamma(J0,JN,M,b) :
        same literal continuation, (e_a,e_g) in E_a(W|past) x E_g(W|past),
        the SAME physical MARINE trajectory and applied MAGNETIC SERVICE }.

The product notation does not authorize independent choice of physical
kinematics, M, b, loss, tuner, covariance or scheduler. SF5 gives signed matrix
support bounds for sensor terms, with primitive endpoint and rotating-weight
variation retained. `linked_supply.imu_supply_outer` supplies this conditional
algebra for both sensors and refuses an unknown fast profile; it is NOT a
source-uniform estimate of chi_*. SF6 gives the joint ambiguity equations.
Raw two-epoch bounds use reachable cell differences, not a fictitious small
noise mean; averaged or transported rows must use their own functionals.

The all-slow sin-cubed family survives for an existing symbolic MARINE pair.
Therefore strict gauge breaking does not follow from a renamed residual class.
No quiet-packet bias-only threshold is a universal physical entry test.

## 3. Exact signed word and optimal linked additive constant

At every operation of ONE realized history define

    e_i = A_i e_(i-1) + d_i,
    M = A_N ... A_1,
    b = sum_i A_N ... A_(i+1) d_i.

This is an identity, with d_i the exact discrepancy from the actual operation.
It includes nonlinear coordinate/reset, physical OU/BA mismatch, projections,
reference changes and arithmetic. No second independently adapted filter is
substituted. Do not replace b by the sum of the norms of its summands.
Let J0=P0^-1, JN=PN^-1. Then

    e_N = M e_0+b,
    D = J0-M' JN M.

For a fixed 0<gamma<1 suppose G=D-gamma J0 is positive definite. Define

    z = M' JN b,
    chi_gamma = b' JN b + z' G^-1 z.                         (LS1)

Completion of the square gives EXACTLY

    V_N = (1-gamma)V_0 + chi_gamma
          -(e_0-G^-1 z)' G (e_0-G^-1 z).                   (LS2)

Therefore

    V_N <= (1-gamma)V_0+chi_gamma,                          (LS3)
    sqrt(V_N) <= sqrt((1-gamma)V_0+chi_gamma)
              <= sqrt(1-gamma)sqrt(V_0)+sqrt(chi_gamma).   (LS4)

For fixed (J0,JN,M,b,gamma), chi_gamma is the SMALLEST additive constant in
LS3 valid for unrestricted e_0: equality is attained at e_0=G^-1 z. On a
physical execution b generally depends on e_0 and the carried history; LS2
is still pointwise exact, but applying it uniformly requires bounding
chi_gamma over the actual linked reachable pairs, not independently selecting b.

For P0=L0 L0' and PN=LN LN', define H=LN^-1 M L0 and beta=LN^-1 b.
Then Delta=I-H'H and the equivalent normalized formula is

    chi_gamma = |beta|^2
      +(H'beta)'(Delta-gamma I)^-1(H'beta).                 (LS5)

This retains the direction of the forcing relative to the SAME loss matrix.
For example H=diag(.99,.5), gamma=.01 and unit forcing in the first/second
output direction give chi=100 and chi=99/74 respectively. Equal forcing norms
do not imply equal robustness. This is precisely what independent worst-case
q and E bounds discard.

If b=B_W w is an exact source parametrization, keep the full matrix

    Xi_gamma = B_W'[JN+JN M G^-1 M' JN]B_W.                 (LS6)

The corresponding block matrix with upper-left G, upper-right -M'JN B_W,
and lower-right Xi_gamma-B_W'JN B_W is positive semidefinite, with zero Schur
complement. This provides a linked matrix supply certificate without a product
of independently optimized gain and loss constants. The physical admissibility
of w, its linkage to M, and any nonlinear remainder remain obligations.

## 4. Boundary retention, strict entry, and prefixes are different tests

If gamma>=gamma_*>0 and chi_gamma<=chi_* on a proved retained class, put
r_*^2=chi_*/gamma_*. Root retention in V<=C follows when chi_*<=gamma_* C.
Finite entry into V<r_in^2 requires the STRICT reserve r_*<r_in, not merely
chi_*<=gamma_* r_in^2. Indeed

    V_n <= r_*^2+(1-gamma_*)^n (V_0-r_*^2).

With a positive gap r_in^2-r_*^2 this crosses r_in^2 in finitely many words.
At equality it may approach the target forever without reaching it.

Endpoint retention alone says nothing about intermediate operations. At each
prefix compute its own M_l, b_l, J_l from the same history. For a fixed b_l,
a sufficient directly linked test for V_0<=C implying V_l<=C_prefix is
existence of mu_l>=0 with

    [[mu_l J0-M_l'J_l M_l,          -M_l'J_l b_l],
     [-b_l'J_l M_l, C_prefix-b_l'J_l b_l-mu_l C]] >= 0.    (LS7)

This is the quadratic implication obtained by subtracting
mu_l(C-V_0). On the actual reachable class the coefficients/defects must again
be treated jointly. LS7 is not asserted verified. The local .15-ball supply
conditions still need every-prefix slack with r_in<.15.

## 5. Executed 80-digit non-promoting feasibility

`linked_supply_source_diagnostic.py` observes the inherited `ag_readout_source.cpp`
fixture through construction, Live, refinement and BA release, then exports
225.00--225.32 s without reseeding. It uses a temporary instrumented header and
an untapped control; all terminal means/covariances/quaternions/stage values
agree exactly. The physical wave is roll=.02 sin(t/2), p_z=.4 sin(.6t), with
its continuous derivatives and one S origin at Live carried into the word.
Source headers are hashed. The exported covariance's exact symmetric part is
used as the metric; this is not a float32 positivity/roundoff enclosure.

For gamma=one half of the smallest normalized loss eigenvalue:

| Quantity | Quiet | Wave |
|---|---:|---:|
| Recorded operations | 281 | 283 |
| Homogeneous rho | .9995812776622473 | .9993070523483296 |
| gamma | .0002093611688763 | .0003464738258352 |
| sqrt(V_root) | 0 | 42.0590575901 |
| sqrt(V_end) | 0 | 45.0966333790 |
| Signed endpoint forcing norm | 0 | 9.58984071777 |
| chi_gamma | 0 | 306.1085110793 |
| chi_gamma/(gamma*.15^2) | 0 | 39266523.7608 |
| Fixed-forcing sufficient radius sqrt(chi/gamma) | 0 | 939.945096597 |

The separated fixed-forcing expression |beta|/(1-sqrt(rho)) is 27673.6030 on
the wave word: linking improves it greatly, but does not make this a .15
retention certificate. LS2's residual is below 4.4e-78. The endpoint defect is
computed retrospectively as e_N-M e_0; it is not an a priori source supply
bound valid when e_0 is changed. The finite wave replay has no all-time magnetic
service certificate. Thus its failed numerical budget is NOT a counterexample
to all admissible entry. No interval refinement of this failed budget is made.

Reproduce:

    python -m tools.stability.ou3_theorem.linked_supply_source_diagnostic \
      --export --eigen /usr/include/eigen3 \
      --directory /tmp/ou3-linked --output /tmp/ou3-linked/report.json

The exact-rational `linked_supply.py` checker separately verifies LS1--LS3,
sharpness, correlation retention, signed source cancellation and fail-closed
relative positivity. These are algebra checks, not a source-uniform theorem.

## 6. Historical V3 norm-only finite-residual obstruction (not new-model admission)

This obstruction does not use the failed wave budget. It uses the literal
norm-only deterministic residual contract: |n_a|<=.3 m/s^2, |n_g|<=.02 rad/s.
The PRE-MIGRATION V3 constants stated no temporal cancellation condition.
`finite_residual_obstruction.py` freezes that historical class explicitly.
The current constants have a different two-timescale contract with OPEN temporal
parameters; the construction below is not thereby admitted to the new class.

Set alpha=.01 rad, nu=.5 /s, c=.01 m/s^2 and, for all time,

    phi(t)=alpha sin(nu t),  R_bw(t)=Rx(phi(t)), B=75 e_x,
    p=v=a=S_true=0,  b_g=0, b_a=c e_x,
    n_g=-phi'(t)e_x,
    n_a=g R_bw(t)'e_z-g_model e_z-c e_x.                   (LS8)

Here g=9.80665 is physical gravity and g_model=9.8066501617431640625 is
separately the literal float parameter value. The exact measured packets are

    gyro=0, accel=-g_model e_z, mag=75 e_x.                (LS9)

The biases are constant, translation/jerk/primitive are zero, and physical
angular speed is at most .005 rad/s. Residuals obey the all-time rational bounds

    |n_a| <= g alpha+c+|g-g_model| < .108067 < .3,
    |n_g| <= alpha nu = .005 < .02.                       (LS10)

The period is 4*pi<60 s and each 60-s window has true gravity span .02 rad,
strictly greater than 1 degree. This satisfies the suggested 1-degree/60-second
EXCITED_MOVING example and also exceeds the earlier bias-only cutoff. Those
round numbers are tested examples, not newly imposed universal qualifications.

The packet history is the existing quiet construction. Consequently Mahony,
tuner, covariance, nominal magnetic reference, accepted updates and release
chronology are identical to that quiet execution; they are not independently
chosen. In this witness actual innovations vanish by the explicit packet and
mean induction, NOT because homogeneous D=0. The nominal LIN and bias means
remain zero. The inherited real-arithmetic quiet-service certificate has
normalized actual-information lower bound 3.720373823... . At each physical
root the true-heading/axial-bias z projection is cos(phi); independent quiet
axis groups add positive semidefinite action. Retaining the z group gives
at least (1-alpha^2/2)^2 times that floor, still greater than one. Thus the
magnetic premise is checked by information, not callback count.

The source BA covariance ceiling is P_ba,ba<=(1/40)^2 I. As b_hat_a=0,

    V >= e_ba' P_ba,ba^-1 e_ba >= (40c)^2=.16,
    sqrt(V)>=.4>.15                                         (LS11)

on the regular A21 continuation. Cross covariance cannot remove this principal
lower bound: it already minimizes storage over all other error coordinates.
The active BA prediction demonstrates the supply mechanism exactly:
phi_b c+(1-phi_b)c=c. Its physical/model mismatch replenishes the homogeneous
OU decay. There is no possible complete-word inequality with a strict uniform
budget chi_gamma<gamma*.15^2 on a class containing LS8--LS11 and guaranteeing
retained full-physical-error entry. This does not refute the conditional local
homogeneous theorem, claim filter divergence, or prove an all-time float32 fact.

`finite_residual_obstruction.py` checks the historical rational bounds and
reuses the existing actual-service certificate. SF7 re-tests LS8 under the
new model: all-slow compensation violates the candidate rate bounds; ANY
mixed split must meet the opposing-window accumulation requirements. Admission
and exclusion remain OPEN until both temporal profiles are independently
qualified. LS11 is about the old chosen slow BA coordinate and is not invariant
under reassignment of an error component to a different slow/fast split.

## 7. Current conclusion

The linked identity and its root/prefix conditions are proved. The diagnostic
improves the separated bound but fails the small-radius feasibility gate.
The HISTORICAL norm-only class contains a quiet-packet history outside the
old inner ball forever. Its new-model mixed allocation is OPEN; the smaller
all-slow sin-cubed ambiguity still survives the symbolic MARINE condition.
Therefore robust physical entry below .15 is not closed by this continuation.
Outer compactness/release and zero-action transfer also remain open following
the prerequisite corrections in section 1. Preserve the existing local result
with its original conditions; do not revive the invalid global extension or
promote any end-to-end flag.
