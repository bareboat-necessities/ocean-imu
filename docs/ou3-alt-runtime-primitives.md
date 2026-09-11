# ALT shared OU runtime primitive lemma

## Scope

The finite prediction path carries one physical segment, one full 21-state
covariance and one joint24 error/true-bias state.  This lemma removes three
previously conditional prediction operands from the higher-level runtime entry:
per-axis OU mean coefficients, active accelerometer-bias decay `phi_hat`, and
active accelerometer-bias process covariance `Q_BB`.

It is an algebraic runtime relation, not COMPLETE-BRMM/BIAS admission, finite
precision, contraction or storage feasibility.

## Integrated OU mean coefficients

For one shipping step `h`, post-clamp time constant `tau`, and the scalar runtime
decay

`alpha = exp(-h/tau)`, `x=h/tau`,

the current shipping implementation uses

`phi_va = tau (1-alpha)`.

For `|x| < 0.01`, `safe_phi_A_coeffs` uses the literal polynomial branch

`phi_pa = tau^2 (x^2/2 - x^3/6 + x^4/24)`,

`phi_Sa = tau^3 (x^3/6 - x^4/24 + x^5/120)`.

Otherwise it uses, with `expm1(-x)=alpha-1`,

`phi_pa = tau^2 (x + alpha - 1)`,

`phi_Sa = tau^3 (x^2/2 - x - alpha + 1)`.

The same tuple `(phi_va,phi_pa,phi_Sa,alpha,h)` is used on all three axes in
the current shared-tau shipping configuration.  `finite_ou_runtime_primitives`
constructs these coefficients and feeds exactly those values into both the
finite mean predictor and `F_LL`; callers cannot supply a second coefficient
set at that entry point.

The remaining source obligation is to bind `alpha` itself to the deployed
`exp(-h/tau)` evaluation and to the same tuner/runtime predecessor.  The finite
real-arithmetic lemma does not claim a binary32 exponential enclosure.

## Accelerometer-bias mean/covariance root

When BA updates are active, shipping uses

`phi_b = exp(-h/tau_b)`

and

`qd_scale = -tau_b/2 expm1(-2h/tau_b)`.

Therefore, in exact real arithmetic,

`qd_scale = tau_b/2 (1-phi_b^2)`.

The runtime primitive constructor uses this identity directly and builds

`Q_BB = Q_bacc * tau_b/2 * (1-phi_b^2)`.

Consequently the finite BA mean/error recurrence, BA covariance block, and both
BA cross-covariance blocks consume one and the same `phi_b`.  H18 is the literal
held branch `phi_b=1`, `Q_BB=0`.  The post-clamp requirement `tau_b>=1e-3` is
retained.

This removes an independent `Q_BB` matrix and an independent `phi_hat` from the
higher-level paired runtime predictor.  It does not yet prove that the deployed
exponential, tuner state, `Q_bacc` and branch guard are reached from every
admitted source history.

## Remaining prediction attachment

The ordinary covariance block identity and this runtime-root lemma now leave the
following prediction obligations open:

- actual `F_AA,Q_AA`, including the constant-rate rotation/B matrix, fast versus
  structured process-noise branch, Simpson evaluations and PSD hygiene;
- analytic `Qaxis` construction and its regularization branch;
- same-history tuner/frontend ancestry for `tau`, `Sigma_aw`, `tau_b`, `Q_bacc`
  and the scalar decays;
- pending `a_w` covariance positive-part inflation, final covariance hygiene and
  periodic S service;
- deployment floating-point residuals.

Until these and the remaining measurement/frontend/hybrid branches are composed
into the complete finite runtime word, storage search remains forbidden and all
ALT final gates remain false.
