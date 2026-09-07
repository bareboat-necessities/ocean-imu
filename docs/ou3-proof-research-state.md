# OU-III proof research state

## Current status and scope

The canonical source is `COMPLETE_SEA3_NORMAL_LIVE_WORD`. The target is
conditional nonlinear finite-window dissipation (P4), then finite capture (P5),
for the declared complete SEA3 family. Production filter code and its quality
thresholds are not to be changed for proof convenience.

| Stage | Current meaning |
| --- | --- |
| Physical SEA0 -> SEA3 inclusion | OPEN: hard phase-continuous shaping/excitation and the joint translational/rotational response from that history remain separate deployment obligations. Compactness itself is not open. |
| P3 conditional SEA3 | CLOSED by the retained complete-word H18/A21 chain at `delta=1e-18`; keep this mathematical chain frozen. |
| P3 deployment | OPEN. `P3_DEPLOYMENT_PASS=false`; conditional P3 does not establish physical left inclusion. |
| P4 | OPEN. Exact covariance identities and pointwise coordinate identities exist; their full nonlinear shipping attachment/transport and strict complete-word dissipation are not closed. |
| P5 | BLOCKED. No finite-capture claim is authorized before strict P4 closure. |

`P3_CANONICAL_PASS` is a deprecated alias of `P3_CONDITIONAL_SEA3_PASS`, not a
claim of deployment closure. The canonical P4 result must keep
`P4_CANONICAL_PASS=false`, `P4_FINITE_WINDOW_CLOSED=false`, and
`P5_MAY_START=false` until the actual word inequality closes. A green
architecture/contract CI job does not mean P4 has passed.

## Source and frozen P3 contract

The same phase-continuous `zeta=(x^s,lambda,z^t,q)` history must generate the
coupled three-partition directional JONSWAP/PM response, translation and
rotation, private front end, WavePeriodEstimator, tuner candidates and staged
commits, committed `tau/sigma/R_S`, pseudo cadence, and shipping `F,Q`.

Retain every valid accelerometer update, every due S=0 update with its actual
applied per-axis SpectralMSE R_S, asynchronous magnetic/vector PE, full process
Q, covariance-floor events, and immediate Joseph/injection/reset ordering.
The deployed R_S standard-deviation factors are `[0.72,0.72,1]`. Target R_S is
not a substitute for the independently smoothed *applied* R_S. Four selected
S records can witness information inside the complete word; they cannot
replace that word or its scheduler.

The H18 prior-free 3-second completion is followed by margin-preserving shipping
events. The H-to-A dimension change is separate; the configured magnetic
refinement hold gives H18 time to close before the A21 continuation. Consume
the complete full-matrix comparison `Omega_W-delta P_W >= 0`, `delta>=1e-18`.
Do not substitute the retired scalar-tube ratio or rework P3 to make P4 pass.

Preserve the separate 0.8-rad outer Cayley geometry and the declared inner
candidates `[30,25,20,15]` degrees. No candidate is currently a certified P4
basin. Keep lever arm zero/disabled and the vibration guard dormant-transparent.

## Retained facts and exact limitation

For a shipping operation `S=H P H^T+R`, `K=P H^T S^-1`, `d=K y`, and
`J=y^T S^-1 y`, the dimension-independent full-matrix inequality is

`d d^T <= J K S K^T <= J P`.

Thus `||E_i d||^2 <= lambda_max(E_i P E_i^T) J`. For `y=H phi`, also
`J<=phi^T P^-1 phi`. These are per-operation identities, not a source-uniform
correction radius or a proof of physical nonlinear storage transport.

With `Q_aw=R_hat^T E R_hat` and `u_aw=Q_aw delta_a_w`, the pointwise residual
rewrite has an a_w-free rotation remainder relative to an **auxiliary H0**
whose a_w column is `R_hat`. Separately, the exact covariance congruence
`T_E=diag(I,...,Q_aw)` requires

`P_u=T_E P T_E^T`, `H_u=H T_E^T`, `K_u=T_E K`.

Its a_w column is `R_hat Q_aw^T`, not `R_hat`. Hence `H0 != H_u` at finite
angle. The residual rewrite and metric congruence are individually valid;
combining them as an unchanged shipping Joseph operation is not established.
The corresponding helper status `shipping_Joseph_binding_closed` is false.

The corrected measurement-linearizing coordinate uses the ORIGINAL shipping
H,P,K,S and the full shift, not a mixture of auxiliary H0 and congruent H_u:

    epsilon_aw=(Q_aw-I)delta_a_w+e_eta,
    e_eta=R_hat^T((E-I)-[c]_x)f_hat,
    Phi(z)=z+E_aw epsilon_aw,
    y=H Phi(z).

This pointwise identity does not make Phi-storage isometric to the original
physical-error storage. For d=K y, t=z-d, z_plus=G t+rho, and A=I-KH, the
correct complete shift transport is

    xi=rho+E_aw(epsilon_plus-epsilon_minus),
    epsilon_mixed_plus-epsilon_mixed_minus
      =(Q_plus-Q_minus)delta_a_w-(Q_plus-I)d_aw.

With u=A Phi and b=G^-1 xi, the exact signed Joseph/reset ledger is

    V_plus-V_minus=-J+2 u^T P_J^-1 b+b^T P_J^-1 b.

Retain the FULL posterior precision P_J^-1, including all cross terms. The
pure e_eta bound alone does not bound epsilon_aw. The nominal specific-force
bound used by the invariant-coordinate helper is derived from the existing
SEA3 true-force ceiling plus the declared `delta_a_w` error ceiling; the old
physical-force-only reuse is forbidden.

The true-minus-estimated rotation error after left quaternion injection obeys
E_plus=E E_correction^-1. Its Cayley numerator is c-a+0.5 a cross c and its
denominator is 1+0.25 a dot c, against the reset target G(d)(c-d).
Exact-rational tests in both 18 and 21 dimensions with correlated covariance
check the full-shift identity, signed energy ledger, and the failure of the
old pure-shift transport. These fixtures are operation tests, NOT legal SEA3
word generators, universal source coverage, or nonlinear word certificates.

Exact covariance congruence is a metric isometry. It is not a Euclidean
covariance upper: `G=I+[d]_x/2` can enlarge transverse covariance. Likewise a
lower on `P_theta,theta` is not an upper on `(P^-1)_theta,theta`; use the full
Loewner bound or the appropriate conditional-covariance Schur complement.
The correction helper therefore emits no purported global energy balls,
source-uniform ceilings, or posterior inverse-metric floors.

## Failure analysis, dead ends, and architecture review

Classification: **proof-method/attachment gap plus point-diagnostic
implementation defects**; neither diagnostic failure is a counterexample to
the complete SEA3 theorem or frozen P3.

The first exact operation-ledger run halted at an accelerometer covariance
reconstruction after only 43 H18 predictions. The cause was a pre-commit
snapshot: the ledger read the pending a_w floor, applied R_S and pseudo cadence
before the wrapper's staged tuner commit, whereas shipping commits the previous
sample's pending tune before the current Riccati prediction. After moving the
host-only snapshot to that exact shipping boundary, the ledger executes the
full word in both modes: 600 predictions, 600 accelerometer updates, the exact
selected S/vector counts, 28 covariance-floor events, and zero telescoping
error.

The repaired point ledger gives

- H18: `rho=0.999863862991`, with 137 due S updates and 75 vector events;
- A21: `rho=0.9958766098`, with 108 due S updates and 75 vector events.

The old point-map scan gives H18 `0.9998645113825368`, agreeing within about
`6.5e-7`, but A21 `0.9958291821145342`, differing by about `4.74e-5`.
Inspection shows the supposedly independent exact-map observer has the **same
pre-commit snapshot defect**: it reads R_S, pseudo cadence, pending floor and
A21 bias-process parameters before `fusion_.update()`, whose first shipping
action is the staged commit. Its covariance endpoints are shipping-real, but
its reconstructed sample map can use stale current-sample schedule data.
Therefore the old `OU3MAP3` point trace is not an admissible independent
reference across tuner-commit samples.

**DEAD END:** duplicated host observers that independently reconstruct the
current-sample schedule before calling the wrapper. The same mechanism failed
twice; do not repair it by maintaining two parallel copies of the schedule
logic.

Architecture alternatives considered after the second strike:

1. Apply the same staged-commit patch to the old exact-map observer. This is the
   smallest diff but preserves two duplicated current-sample schedule models
   and therefore shares the failed mechanism.
2. Make the corrected **single shipping operation observer** the sole point-map
   source. It propagates the full 21x21 deterministic map and records block
   boundary shipping covariances while using the same exact current-sample
   schedule boundary as its signed event ledger. This removes the duplicate
   schedule model and directly serves the master complete-word ratio.
3. Reconstruct each current sample only after shipping from scratch matrices
   and covariance deltas. This avoids private staged-commit access, but the
   pending covariance-floor request has already been consumed and would need
   additional production-adjacent instrumentation.

Choose alternative 2. The old point-map observer may remain for unrelated
legacy diagnostics, but it must not be the reference for this P4 feasibility
experiment.

Invalidated proof claims remain: pure e_eta transport dropping
`(Q_aw-I)delta_a_w`; zero-cost combination of H0 with the H_u covariance/gain;
prediction-only covariance ceilings through nonorthogonal resets; reusing an
H18 nonlinear ceiling in A21 without transport; selecting the minimum posterior
floor from the cell minimizing a different scalar ratio; and reading inverse
metric bounds from a marginal covariance lower. Unsupported numerical
reset-radius closure has been withdrawn, not repaired by a smaller domain.

Retained: the exact residual algebra, individual congruence identities, pure
vector Cayley identities, magnetometer radial cancellation, the full
per-operation correction inequality, the repaired single-estimator event
ledger, actual shipping covariance endpoints, and the frozen complete-word P3
result. The 30-degree information-headroom calculation is not signed-word
dissipation.

Do not revive replay fitting, arbitrary bounded-input/source boxes, independent
`tau x sigma x R_S x T_S` or sea x RAO products, endpoint/history graphs,
selected-four-S replacement words, scalar information-beta, blockwise minimum
ratios, or packet-count-times-worst-remainder bounds. No new micro-certificate
should be added without its place and quantitative effect in the master word
inequality being demonstrated.

## Next falsifiable experiment

The controlling object is `rho_W=V_after(F_W(x))/V_before(x)` on the actual
complete H18/A21 word, not positivity of an isolated lemma.

Extend the corrected single shipping operation observer with a scan mode that
propagates the full 21x21 map and records its own P0/P1 covariance at every
600-sample block boundary. Preserve the existing `OU3MAP3`/`OU3COV1` binary
contract if practical so the current generalized-eigenvalue analyzer can be
reused. Run the genuine same-history source, select the worst legal H18/A21
blocks from that one observer, then replay each maximizing direction through
the same observer's event ledger. Require, for each mode:

- 600 predictions and 600 accelerometer updates;
- exact due-S and vector event counts;
- exact telescoping of event `delta_V`;
- identical actual-R_S anisotropy/history;
- event-ledger rho and block-map generalized rho agreeing to numerical
  roundoff.

If that self-check fails, stop and analyze the mismatch. Do not build a
nonlinear enclosure on top of inconsistent point-map machinery.

Only after the single-observer point experiment passes should the three P4
treatments be compared: (1) unchanged shipping coordinates with the exact
signed full residual/Joseph/reset word; (2) corrected full-Phi nonlinear
storage with all mixed transport and uniform metric comparison; (3) full-word
Schur elimination of the coupled a_w direction while retaining every actual-R_S
S event. The existing interval-AD primitive is appropriate only if this point
experiment justifies a complete-word Jacobian/transport enclosure; do not fall
back to packetwise scalar remainder sums.

Keep the complete-word feasibility test non-promoting and same-history. Report
both H18/A21 worst ratios, source/phase, maximizing error direction,
operation-by-operation margin consumption, and distance to rho=1. A point
sample is a falsification/feasibility diagnostic, never universal source
coverage. Do not silently replace a missing source realization with independent
bounded sequences. Follow `EXECUTE -> FAILURE ANALYSIS -> REPLAN -> EXECUTE`.

Only a strict complete-word result may authorize enclosure of the nonlinear
transport, selection of the widest certified candidate, and then P5 finite
capture including the separate H-to-A hybrid event and forcing terms.

## Code map and reproduction

Retained proof scripts and committed proof JSON inputs live in
`tools/stability/`. The domain files are
`tools/stability/ou3_proof_operating_domain.json` and
`tools/stability/ou3_sea3_directional_response_domain.json`; the committed bridge is
`tools/stability/ou3_sea3_spectral_moment_bridge.json`. Relocation changes no JSON bytes.

Entry points in that directory: `ou3_sea3_riccati_metric_p3.py` is the frozen
conditional gate; `ou3_sea3_riccati_metric_p4.py` is the non-promoted canonical
P4 status. The `ou3_p4_complete_sea3_*coordinate.py` helpers contain the trial
coordinate algebra. `ou3_p4_complete_sea3_correction_information_bound.py`
contains the exact operation identity, full-posterior precision evaluator,
and explicit open obligations; candidate metric-energy radii remain null.
`ou3_p4_exact_reset_transport.py` owns the parameterized reset remainder.

Run `make -C tests/validation build`, the source-foundation unittest command in
`.github/workflows/ou3-proof.yml`, and all `test_ou3_p4_*.py` tests. Run each
canonical gate with `--output /tmp/<name>.json`; inspect its pass flags, not
just its exit code. `python3 tools/ou_evidence_contract.py --check` is separate
from the stability theorem. Fresh validation/robustness replay evidence is
still required for the changed source inherited from main; never hand-edit
provenance hashes. The opt-in `ou-full-evidence-branch.yml` workflow owns a full
branch refresh. No filter or pitch-gate edit is part of this handoff.

## Next conversation / PR boundary

PR #496 is the complete-word P4 feasibility/closure continuation from the
merged #494 handoff. Keep P3 frozen and P4/P5 open until the actual nonlinear
word inequality closes. `tests/validation/test_ou3_p4_coordinate_transport_algebra.py`
locks the exact full-shift, reset, and correlated-precision corrections. The
full source-correlated H18/A21 transport and storage comparison remain the
decisive obligations; green algebra/architecture tests do not close them.
