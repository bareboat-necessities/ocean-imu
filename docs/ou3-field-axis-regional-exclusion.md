# OU-III: source-audited regional field-axis exclusion

Source snapshot: PR #637, `1bbd3b3bc58879d7f95863cef5726ee10c80616f`.
This note supersedes QR's 25-Hz accelerometer substitution and its use of a
congruent-sync covariance ceiling for the default additive-floor implementation.
It establishes a conditional regional geometric bound, not capture, prefix
retention, a complete-word contraction coefficient, or end-to-end stability.
No runtime, physical assumption, tuner law, cadence, or quality gate is changed.

## 1. Exact object and exact epoch

In `measurement_update_acc_only`, the implemented attitude Jacobian is

    J_theta,i = -[R_wb,i (a_hat_i - g_model)]_x.

Here `a_hat_i` is the nominal AW state AFTER prediction and any due S correction,
but BEFORE this accelerometer correction. Bias and lever-arm terms are not in
this attitude Jacobian. The magnetic attitude row uses `R_wb,i B_ref,i`, where
`B_ref,i` is the actually stored nominal reference, not automatically the true
physical field. Write `b_i = B_ref,i / ||B_ref,i||` and `P_b = I - b b'`.
The field-axis unit attitude vector has accelerometer-row magnitude exactly

    ||J_theta,i R_wb,i b_i|| = ||P_bi (a_hat_i - g_model)|| =: c_i.       (FA1)

Thus the particular pathological base execution has `c_i = 0` at every relevant
accelerometer epoch. This is a statement about the BASE nominal geometry. A
homogeneous test vector, the actual base estimation error, and a base innovation
are three different objects; zero homogeneous action does NOT zero base
innovations. We do not identify all possible zero-action modes with FA1 here.

The wrapper calls time_update and measurement_update_acc_only at EVERY valid
MEKF-driven IMU sample. On the qualified regular branch, applied accelerometer
gaps are at most 0.006 s. The 25-Hz example belongs to magnetometer service.
An invalid input or failed solve must be charged through the ACTUAL applied
accelerometer gaps; invocation alone is not an accepted-update certificate.

## 2. A marginal bound for the literal default covariance path

The default wrapper has `S_factor_ = 1`, so the stationary AW covariance and
every queued sync target are isotropic. Its default sync is NOT congruent:

    A <- A + (s^2 I - A)_+,       A = P_aw,aw.                         (FA2)

If `A = U diag(lambda_j) U'`, FA2 gives
`U diag(max(lambda_j,s^2)) U'`. Therefore it preserves `A <= C I` whenever
`A <= C I` and `s^2 <= C`. This argument does not require an upper bound on the
full 21-state covariance or observability of the field-axis attitude mode.
Do not extend this commuting-target argument to arbitrary anisotropic targets.

The AW row of prediction is independent: `A^- = phi^2 A^+ + Q_aa`.
Accepted Joseph corrections decrease this marginal, because its decrement is
`PCt_aw S_innov^-1 PCt_aw' >= 0`. The quaternion reset and body-prime retargeting
leave the world-frame AW marginal unchanged. Live initialization explicitly
seats it at the stationary covariance and clears its cross terms.

For an exact OU covariance, `Q_aa = sigma^2(1-phi^2)`. The literal small-x source
uses polynomials, so an exact ceiling of 16 should not be imported without
checking that branch and its two regularizations. A simple conservative audit
of the actual polynomials gives, in regular REAL arithmetic,

    Q_aa <= 1.03 sigma^2 (1-phi^2) I.                                (FA3)

Here is a reproducible uniform bound, rather than a parameter grid. Put
`x=h/tau`, `h<=.006`; the polynomial branch has `0<x<.01`. For each scalar
entry factor out `sigma^2 h^p x`, replace coefficients by absolute values,
and divide by `2(1-.01)`, using `1-exp(-2x) >= 2x(1-x)`.
All remaining powers of h and x are nonnegative, so evaluation at the two
upper endpoints bounds the entire branch. For the initial [v,p,a] marginal,
use its largest absolute row sum. Spectral clipping cannot increase operator
norm. Adding the S row increases the bound by at most
`|q_vS|+|q_pS|+|q_Sa|+|q_SS|`; the second regularization again cannot increase
operator norm. The resulting EXACT RATIONAL bound is less than
`1.023336613 < 1.03`. The closed-form branch is the exact covariance in real
arithmetic. This audits the literal source formula, including its polynomial
branch; it does not certify IEEE-754 sanitization, eigensolver failures, or
roundoff over an infinite execution.

Consequently, for the already qualified default profile `sigma_aw<=4`, including
its applied noise floor, and a Live seed `A<=16I`, every regular prefix obeys

    P_aw,aw <= 16.48 I.                                              (FA4)

A different configured amplitude ceiling, nonisotropic profile, or numerical
fallback needs its own bound. The full coupled tau/sigma_aw/r_S/T_S history
remains the shipping one. FA4 is a necessary bound valid on that history; it
is NOT a license to choose the four channels independently in later arguments.

For the actual base error e, with SPD covariance P at the SAME prefix, define
`V_base=e'P^-1 e`. Matrix Cauchy--Schwarz, retaining every cross covariance,
gives

    ||a_hat-a_phys|| <= sqrt(16.48) sqrt(V_base) < 4.06 sqrt(V_base).  (FA5)

This is an algebraic component bound, not a theorem that V_base decreases.
The assumed/proved retained ball must cover all pre-accelerometer prefixes,
not only word endpoints or post-correction epochs.

## 3. Sharp same-history geometric excitation inequality

First take a fixed nominal reference b and fixed model gravity g. Define

    G=||P_b g|| > 0,       u=P_b g/G.

Let `t_0<...<t_N`, `d_i=t_(i+1)-t_i`, `L=sum d_i`. These are the ACTUAL applied
accelerometer epochs from one physical history. Its acceleration is the
derivative of velocity, `||v||<=Vmax`, `||a_dot||<=Jmax` almost everywhere.
Let w be normalized trapezoidal endpoint weights:

    w_0=d_0/(2L), w_N=d_(N-1)/(2L),
    w_i=(d_(i-1)+d_i)/(2L),               sum w_i=1.

Cellwise integration by parts proves

    |sum w_i u'a_phys(t_i)|
       <= 2 Vmax/L + Jmax/(4L) sum d_i^2 =: C_phys.                  (FA6)

Indeed the cell error is `int_0^d (s-d/2) u'a_dot(s) ds`, with bound `J d^2/4`.
This is sharp for a triangular tent. These are proof weights, not a new
integrator, resampler, or pseudo-measurement cadence.

If `V_base,i <= r^2` at all these prefixes, then FA1 and FA5 give

    u'a_phys(t_i) >= G - c_i - kappa_a r,   kappa_a=4.06.

Averaging and using FA6 yields the quantitative bound

    sum w_i c_i >= G - kappa_a r - C_phys =: M(r),
    sum w_i c_i^2 >= max(M(r),0)^2.                                 (FA7)

The second line is weighted Cauchy--Schwarz. This proves more than the absence
of exact collinearity: it gives a strictly positive BASE geometric row-energy
floor whenever M(r)>0. It is NOT yet the complete innovation-weighted Schur
loss, which must also account for time transport and all nuisance states.

With maximum applied gap h_acc, `sum d_i^2/L <= h_acc`. Exact field alignment
cannot persist over a complete-cell interval with

    L > 2 Vmax / (G - kappa_a r - Jmax h_acc/4),                      (FA8)

provided the denominator is positive. Neither translation excitation nor a
minimum wave amplitude is used: this particular exclusion also covers still
water within the same error ball.

A sharper, still same-history version replaces kappa_a in the fixed-reference
argument by `sum w_i sqrt(u'P_aw,aw,i u)`. The scalar global ceiling merely
makes FA7 source-uniform without solving the full Riccati-return problem.

## 4. Variable/biased nominal reference: do not substitute the physical field

For a fixed comparison axis b0 and gravity g0, let

    G0=||P_b0 g0||,  u=P_b0 g0/G0,
    ||g_model-g0|| <= delta_g,
    ||P_bi-P_b0|| <= eta_ref <= 1,
    ||a_phys|| <= Amax,  ||g_model|| <= g_model_max.

Decompose `a_hat_i-g_model = P_bi(a_hat_i-g_model)+lambda_i b_i`.
Since `|u'b_i|<=eta_ref` and
`|lambda_i|<=Amax+kappa_a r+g_model_max`, the same proof gives

    sum w_i c_i^2 >= max(M_ref(r),0)^2,
    M_ref(r)=G0-delta_g-eta_ref(Amax+g_model_max)
             -kappa_a(1+eta_ref)r-C_phys.                            (FA9)

Thus a fully explicit sufficient radius is

    0 < r < [G0-delta_g-eta_ref(Amax+g_model_max)-C_phys]
             / [kappa_a(1+eta_ref)],                                (FA10)

when the numerator is strictly positive. A scalar additional geometric defect
is subtracted in the same numerator. Sensor bias and lever arms do not appear
as independent FA1 terms; their effects on reference construction, retention,
and physical-model consistency cannot be discarded from those other proofs.
The physical field's epsilon_B alone is NOT a bound on the learned reference's
eta_ref. They must be connected through the literal acquisition/refinement and
continuous hard-iron chronology, or the bound must be stated directly for the
actually stored reference.

## 5. Certified substitutions and scope

With the declared `Vmax=5.5`, `Jmax=100`, `h_acc<=.006`, and `L>=100`,

    C_phys <= .11+.15 = .26 m/s^2.                                 (FA11)

Use the actual complete-cell duration, not a rounded clock label. Inside an
arbitrarily aligned 100-s proof word, trimming to accepted endpoints leaves
L>=99.988 s when both boundary gaps are at most .006 s. The same formula then
uses 11/99.988+.15 instead of .26. This does not change the 100-s proof clock.

For the fixed-reference specialization with `g=9.80665` and a certified nominal
inclination at most 80 degrees, `G>=g cos(80 deg)=1.7029069015...`.
The pre-radius margin is therefore `1.4429069015...`, not QR's `0.5929069015...`.
At the explicit Mahalanobis radius `r=1/4`, with kappa_a=4.06,

    M(1/4) > .4279 m/s^2.                                          (FA12)

With the just-described boundary trimming, the certified bound is >.4278
m/s^2 instead; the loss is less than .000014 m/s^2.

The regression encloses sin(10 degrees) with exact-rational Machin/arctangent
and sine series; it does not label an ordinary floating-point evaluation an
interval proof. The illustrative radius is a FULL-STATE COVARIANCE-METRIC
radius, not an angle or a replacement for the six-degree attitude domain.
For variable references FA12 remains positive only after the FA9 deductions
have been bounded; unspecified terms are not assigned zero.

The source also offers a weaker direct startup fact: MagAutoTuner requires
`horizontal_fraction >= .05`, hence a fixed accepted startup reference has
`G>=.05 g`. Under the same 100-s bounds, r=.04 gives M>.0679 m/s^2 without
assuming an 80-degree nominal inclination. But the default continuous
hard-iron tracker is enabled, rewrites the reference, and does NOT reapply this
fraction gate. Its setter only checks finite values and a positive horizontal
magnitude against the configured norm threshold. Therefore neither the 5%
startup cone nor the physical 80-degree cone has been proved an all-time
nominal-reference cone by this calculation. Disabling that tracker or adding a
new gate would change runtime and is NOT proposed here.

## 6. Status and exact next obligation

Established: FA2--FA5 for the stated regular default covariance profile; the
same-history sharp inequality FA7, its variable-reference form FA9, and the
conditional substitutions FA11--FA12. These exclude the exact field-aligned
candidate on any execution satisfying their prefix and reference premises.

Not established: universal exclusion on every shipping-reachable execution,
entry and retention in this ball, exclusion of every other zero-action mode,
source-uniform complete-word contraction, nonlinear disturbance closure, or
IEEE-754 stability. Do not set any theorem-success flag from these checks.

The next calculation is to bound eta_ref (or an equally strong direct nominal
cone plus its variation) from the literal magnetic-reference chronology, then
combine FA9 with the ACTUAL prefix-retention/storage inequality and the coupled
adaptive schedule. Merely reasserting a true-field inclination bound, or citing
visual AW tracking, does not discharge those two premises. No shipping
counterexample has been constructed by pointing out this missing implication.

Validation: 13 new exact scalar/rational regressions pass. The full existing
native/CI suites were not run in this continuation. Earlier O1/O2 results and
failures remain in the archived ledger; none is promoted by this note.
