# OU-III P4 parallel route: path-dependent whole-word contraction

This route is developed in parallel with the operation-matched Joseph ledger.
It is intentionally not a minor variant of that route.

The object being certified here is the **entire fixed-mode source word**.  For a
reachable source-word edge from source node `g` to source node `h`, let

- `F_w(x)` be the exact nonlinear H- or A-mode one-word return map in Cayley and
  physical linear coordinates;
- `J_w(x)` be its Jacobian, or generalized Jacobian where the deployed A-mode
  accelerometer-bias ball projection is active;
- `M_g` and `M_h` be the source covariance-information metrics at the actual
  start and endpoint source nodes.

The target certificate is

`J_w(x)^T M_h J_w(x) <= gamma^2 M_g`, with `gamma < 1`,

for every admitted state `x` in the useful finite-angle/full-state domain and
every required reachable source-word edge.

Equivalently, if `L_g^T L_g=M_g` and `L_h^T L_h=M_h`, the numerical backend may
certify

`||L_h J_w(x) L_g^-1||_2 < 1`.

The factors must be validated and must correspond to the actual source nodes.
A node-specific arbitrary rescaling is forbidden.  The existing P4 metric uses
one global positive normalization per H/A mode, so source-node switching cannot
manufacture contraction.

## Why this is genuinely different

This route does **not** ask for a useful scalar decrease from each measurement.
It therefore does not care that a same-sample vector packet can be rank
deficient.  Rank may be gained by dynamics and later measurements before the
word endpoint.

It also does not require an accelerometer correction to remain inside the
0.80-rad attitude sector by itself.  An intermediate coordinate may grow as long
as the exact complete return map contracts in the endpoint information metric.

The route avoids four scalarization losses at once:

1. no `N * global Lipschitz defect`;
2. no packet-by-packet minimum information eigenvalue;
3. no translation/nontranslation Schur split in the final acceptance test; and
4. no conversion of the P3 source-uniform Riccati `delta` into a nonlinear state
   radius.

All 18 H coordinates or all 21 A coordinates are tested in one matrix
inequality, including attitude-linear and translation/nontranslation cross
terms.

## Why P3 delta is not the nonlinear budget

The current P3 word certificate gives a strict source-uniform linear endpoint
margin.  That proves the linearized origin is on the stable side of the
boundary, but its worst source-uniform scalar `delta` is extremely small.  If
one converts it to a norm gap,

`1 - sqrt(1-delta)`, 

the result is again microscopic.  This route records that value only as a
**diagnostic fallback**.  It does not require the nonlinear Jacobian perturbation
to fit inside that scalar gap.

Instead, the future backend computes the actual complete-word Jacobian on each
reachable source cell and evaluates the full generalized matrix inequality
directly.  A source-uniform scalar is taken only after all within-word
cancellations and cross-coordinate effects have been retained.

## Metric choice

The route reuses the exact Cayley-lifted source information metric already bound
to P3:

`M_g = s_mode * Sigma_KF(g)^-1`.

Important properties are preserved:

- one fixed positive scale for every H node and one for every A node;
- full attitude-linear cross terms;
- endpoint metric tied to the actual endpoint source node;
- no block-diagonal surrogate; and
- no arbitrary path/node scale factors.

Thus a path-dependent metric is allowed only because the source state itself is
part of the theorem state.  The metric switching is source-correlated, not
chosen after observing which metric gives the smallest number.

## Exact map that must be differentiated

The one-word derivative must follow the shipping operation order, including:

1. previous tune commit;
2. dormant transparent vibration-guard branch;
3. prediction;
4. pending `a_w` covariance increment;
5. periodic `S=0` correction or identity, followed immediately by quaternion
   injection/reset when accepted;
6. accelerometer accepted/rejected branch and immediate reset when accepted;
7. tuner evolution/staging;
8. `a_w` covariance-sync staging; and
9. asynchronous magnetometer accepted/rejected/not-due branch and immediate
   reset when accepted.

The fixed H/A word does not contain dimension-changing hybrid events.

For A mode, the deployed accelerometer-bias projection is nonsmooth at its ball
boundary.  The backend must therefore include a **generalized Jacobian** for
that projection.  Euclidean nonexpansiveness alone is not silently promoted to
nonexpansiveness in the full covariance-information metric.

## Numerical backend design

The intended backend should use interval automatic differentiation or a
validated Taylor model over source/state cells.  For each complete-word edge:

1. propagate the exact state and derivative together through every operation;
2. retain full 18x18 or 21x21 derivative matrices;
3. keep sequential reset derivatives instead of merging resets;
4. retain all cross terms;
5. outward-enclose the endpoint derivative over the full cell;
6. form the generalized contraction matrix using the actual `M_g` and `M_h`;
7. validate its largest generalized eigenvalue, or an outward spectral-norm
   equivalent, below one; and
8. refine only cells that fail the interval test, without shrinking the theorem
   domain.

The expensive partition should be source/reachability aware.  Unreachable
Cartesian combinations of tuner states must not be generated just to make the
interval problem larger.

## Current scope preserved

The branch keeps the same deployment-facing scope as main:

- P5 gauged entrance: **45 deg** SO(3) geodesic;
- useful outer attitude geometry: exactly **0.80 rad / 45.84 deg**;
- no `a_w/sigma_applied` consistency assumption;
- no filter/gain/schedule/adaptation-law change;
- no trajectory fit or replay promotion;
- raw tuner sigma states below the filter-side 0.05 m/s^2 floor remain in P2;
- no unqualified `powf/sqrtf` R_S path tightening; and
- the active/transitioning vibration guard remains a separate hybrid/source
  obligation.

## Promotion rule

`ou3_p4_path_dependent_word_contraction.py` is a route contract and diagnostic.
It can validate the metadata contract of a future outward-Jacobian producer but
**cannot promote P4 itself**.

A future numerical producer may promote this route only after it supplies all
of the following:

- source-only, trajectory-independent cells;
- all required reachable/recurrent source-word edges;
- exact shipping complete-word operation order;
- full 18/21-state derivative matrices;
- sequential reset derivatives;
- accepted/rejected/not-due branches;
- A-mode bias-projection generalized Jacobians;
- actual endpoint source-node metrics;
- outward validation; and
- `max ||L_h J_w L_g^-1||_2 < 1` for H and A.

If this route closes, it yields a finite-angle P4 certificate without relying on
per-operation information scalarization.  P5 can then use the certified
complete-word contraction factor to derive finite capture from the 45-degree
outer entrance toward the inner stochastic localization set.
