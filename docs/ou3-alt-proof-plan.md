# OU-III ALT proof execution plan — theorem work only

This plan is normative for the parallel ALT route.  It exists to prevent
research time being spent on attractive diagnostics that cannot discharge the
shipping stability theorem.  The original P2/P3/P4/P5 route remains separate
and independently continuable.

## Target theorem chain

ALT work proceeds only through this dependency chain:

1. **Physical finite error map.**  Materialize the true-minus-estimate shipping
   map on the corrected COMPLETE-BRMM source, including model mismatch,
   physical-reference forcing, BIAS0/1/2 ancestry, actual covariance/frontend
   coefficients and every shipping branch.
2. **Literal source-uniform Live word.**  Compose every prediction, floor,
   due/rejected S update, accelerometer, magnetometer, reset, projection,
   clamp/acceptance branch and H18/A21 mode edge on one same-history source.
3. **Common joint24 storage feasibility.**  Only after (2), search one common
   coercive joint24 metric with bounded neutral/source supply.  Use rank-three
   arithmetic as an exact implementation optimization, never as a state
   reduction.
4. **Rigorous source-uniform certificate.**  Outward-check the common storage,
   all source sectors, endpoint and every-prefix inequalities, metric
   compatibility, first-exit/chart retention and finite precision.
5. **Live basin and hybrid continuation.**  Certify H18, A21, H18->A21 and
   every allowed subsequent edge indefinitely.
6. **Startup capture.**  Separately prove Mahony/proxy/learning reaches the
   certified Live basin for both measured-period takeover and prior-frequency
   timeout.  Only then compose an end-to-end theorem.

No later phase may be started because an earlier one “looks numerically good.”

## Mandatory anti-dead-end guards

Every ALT task, commit and PR must name the theorem obligation it changes from
OPEN toward CLOSED.  If it cannot do that, it is not ALT proof work.

### G1 — Physical admissibility

Proof evidence must quantify the corrected COMPLETE-BRMM family and one of the
admitted BIAS0/BIAS1/BIAS2 histories.  Replay, finite seeds, synthetic
unreachable state/covariance perturbations and captured trajectories may only
falsify an algebraic implementation; they can never qualify a theorem set,
metric, source bound or basin.

### G2 — Same-history ancestry

Physical source, true bias, estimator state, covariance, frontend/tuner,
scheduler/guards, geometry, actual R_S and every event coefficient must descend
from one source lineage.  Independent P/H/R/K, f/sigma/tau/T_S/R_S, source
moment or guard boxes are forbidden unless a proved outer-relation theorem
explicitly permits the relaxation.

### G3 — No homogeneous-truth shortcut

A filter model is not a physical truth model.  Prediction must retain the exact
same-history defect.  For translation this includes

`u=[J0-phi_va a0, J1-phi_pa a0, J2-phi_Sa a0, a1-alpha a0]`.

The existing joint 15D `(a0,a1,J0,J1,J2)` sector is the required source object;
its moments may not be Cartesianized.  A homogeneous error prediction is
allowed only where a source identity proves its defect is zero.

### G4 — Physical S=0 residual

For `e_S=S_phys-S_hat`, the shipping pseudo-measurement residual is

`r_S = -S_hat = e_S - S_phys`,

not `e_S`.  `S_phys` must retain the one-time Live origin and physical primitive
ancestry.  Wordwise re-zeroing, the legacy 300 m*s error ball, or silently
setting `S_phys=0` is forbidden.

### G5 — Bias is one physical history

The true bias and bias-error recurrences share one driver/root.  Equivalently,
for an estimator factor `phi_hat`,

`e_ba+ = phi_hat e_ba + beta+ - phi_hat beta`.

H18 hold is `phi_hat=1`.  H18 accelerometer residuals retain the held
`e_ba`; A21 retains the projection and true-bias cross information.  Two
independent bias drivers or fresh per-sample truth slots are forbidden.

### G6 — Complete word before storage

No rho search, common metric optimization, parameter-dependent metric search,
or high-precision contraction calculation is allowed until the literal
physical word graph contains every shipping event/branch and its source
forcing.  High precision cannot repair a missing term.

### G7 — Common storage first

Once the master is complete, search a single common joint24 storage with
bounded neutral/source supply first.  Metrics may not be fitted to replay.
If two source-uniform attempts fail by the same limiting mechanism, stop and
perform the AGENTS architecture review: record the failure, compare at least
three qualitatively different alternatives, and only then consider
source-dependent or piecewise storage.

### G8 — Rank-three means arithmetic, not theorem reduction

Use `Psi-K(H Psi)`, thin Joseph factors, 3D innovation variables and rank-at-most
six storage corrections where exactly equivalent.  Never drop state
coordinates, covariance/source dependencies, motion-bias cross terms,
`delta N*q`, `delta S*q`, or finite reset/projection coupling merely because a
measurement has dimension three.

### G9 — Hybrid and prefix closure

Endpoint decay alone cannot promote ALT.  Every literal prefix must remain in
its nonlinear/chart domain, every accepted/rejected branch must be covered,
and H18/A21 plus H18->A21 storage/guard compatibility must close.  No hidden
metric restart, bias root restart or S-origin restart is permitted.

### G10 — Finite precision is a promotion prerequisite

Real-arithmetic identities are intermediate lemmas only.  ALT_LIVE_PASS cannot
become true until actual binary32 Eigen/LDLT, Joseph, reset, projection, floor
and relevant frontend numerical defects have an outward additive enclosure.

### G11 — Startup remains downstream of a real Live basin

Live-word proof work must not import Mahony startup or the old route's PE/vector
mismatch.  Startup is addressed only after a quantitative ALT Live basin exists,
and must cover both measured-period takeover and prior-frequency timeout.

### G12 — Stop rule

Before executing a computation, answer both questions:

1. Which theorem obligation can this make strictly stronger or close?
2. Is every input an admitted analytic/source-uniform object, or is the result
   only diagnostic?

If (1) has no concrete answer, do not run it.  If (2) is diagnostic-only, run it
only when it can falsify the *current complete theorem candidate*; never open a
proof PR whose main deliverable is such a diagnostic.

### G13 — A derivative cocycle is not a finite physical master

The six Phase-1 assembly flags are necessary bookkeeping, not an equality for
finite errors. Before a common-M, rho, high-precision or endpoint-refinement
attempt, `proof_plan.assert_finite_storage_master` must validate the declared
representation. A Jacobian product, a source token, or a collection of TRUE
metadata fields cannot replace an exact finite graph or a proved anchored
mean-value representation. Nonzero physical S and other reference offsets must
remain present, as must all coefficient product and branch graphs.

The common-metric entry point and its candidate helper invoke this guard before
building a covariance outer or loading fresh-entry scales. They reject the
current `physical_word` Jacobian cocycle. Substituting a thinner gain enclosure
or a different metric does not discharge this missing identity.

## Current physical-map obligation

The accepted, finite, real-arithmetic accelerometer, magnetometer and S graphs
are proved in [the finite measurement note](ou3-alt-finite-measurement-proof.md).
They include the exact Cayley injection, both quaternion branches, same-beta
radial projection and a thin-factor full joint24 storage identity. In held mode,
BA uncertainty remains in the innovation: the reduced accelerometer covariance
uses R_acc+B0, not bare R_acc. This is not an error-state reduction.

The physical predictor in `finite_physical_prediction.py` retains the actual
physical quaternion increment and its nonzero defect, correlated translation
moments and the shared bias recurrence. All-coefficient polynomial checks prove
the local finite algebra. The older sampled-shadow predictor is not a physical
truth model. Passive shipping traces check implementation correspondence only;
they do not satisfy G2, G6 or G13 or authorize a storage search. The branch-correct
rank-three Joseph identity also remains conditional on the actual inverse-free
solve; it cannot cancel a deployment solve defect.

The complete finite word remains OPEN. Connect these graphs to the actual
covariance/frontend successors, preserve physical continuous-rotation versus
sampled-gyro prediction forcing, and cover every configured accepted/rejected
and asynchronous branch with the same source/primitive/bias ancestry. Either
compose finite descriptors directly or prove the complete anchored mean-value
bridge. Do not present pointwise Jacobian products as the finite endpoint map.
Only a closed finite master may enter the first common joint24 storage search.

## Promotion state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Existing `P4_PASS` / `P5_MAY_START` are not controlled by ALT and remain
untouched by this plan.
