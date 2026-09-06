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
correct measurement/reset shift transport is

    xi=rho+E_aw(epsilon_plus-epsilon_minus),
    epsilon_mixed_plus-epsilon_mixed_minus
      =(Q_plus-Q_minus)delta_a_w-(Q_plus-I)d_aw.

Prediction is different. For the physical map `z_plus=F z+rho_p`, retain

    Phi_plus=F Phi+xi_p,
    xi_p=rho_p+E_aw epsilon_plus-F E_aw epsilon_minus.

`F E_aw` has nonzero v,p,S as well as a_w rows. Substituting the reset-only
shift difference loses the path into the S=0 regularizer. The supplied full
physical defect must include source/model and chart terms; it is not declared
zero by the transport evaluator. Source-only events use their actual map, and
H18->A21 uses the actual rectangular lift and separate covariance/seed operation.
No H18 covariance ceiling is transferred to A21 by this identity.

With u=A Phi and b=G^-1 xi, the exact signed Joseph/reset ledger is

    V_plus-V_minus=-J+2 u^T P_J^-1 b+b^T P_J^-1 b.

Retain the FULL posterior precision P_J^-1, including all cross terms. The
pure e_eta bound alone does not bound epsilon_aw. A physical specific-force
bound does not automatically bound nominal f_hat: the nominal-force/nominal-a_w
bounds remain explicitly conditional, not certified from physical SEA3 alone.

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

## Failure analysis and rejected uses

Classification: **proof-method/attachment gap**, plus implementation/contract
errors; this is not a counterexample to the complete SEA3 theorem or frozen P3.

Invalidated claims: pure e_eta transport dropping (Q_aw-I)delta_a_w;
zero-cost combination of H0 with the H_u covariance/gain;
prediction-only covariance ceilings through nonorthogonal resets; reusing an
H18 nonlinear ceiling in A21 without transport; selecting the minimum posterior
floor from the cell minimizing a different scalar ratio; and reading inverse
metric bounds from a marginal covariance lower. Unsupported numerical
reset-radius closure has been withdrawn, not repaired by a smaller domain.

Retained: the exact residual algebra, individual congruence identities, pure
vector Cayley identities, magnetometer radial cancellation, the full
per-operation correction inequality, and the frozen complete-word P3 result.
The 30-degree information-headroom calculation is not signed-word dissipation.

Do not revive replay fitting, arbitrary bounded-input/source boxes, independent
`tau x sigma x R_S x T_S` or sea x RAO products, endpoint/history graphs,
selected-four-S replacement words, scalar information-beta, blockwise minimum
ratios, or packet-count-times-worst-remainder bounds. No new micro-certificate
should be added without its place and quantitative effect in the master word
inequality being demonstrated.

## Complete-word feasibility admission: current limiter

The numerical experiment has not yet executed an admitted SEA3 word. The
tracked JSON inputs contain no 601-transition numerical witness with its
same-history frontend entry and Live covariance seed. The retained
`ou3_sea3_complete_window_executor.py` rejects materialization while the
code-owned hard-window provider is open. H18/A21 rho, maximizing direction,
source/phase and signed word margins therefore remain **unmeasured**, not PASS
and not evidence that rho exceeds one.

Classification: an executable source-admission gap, not a counterexample to
conditional P3 or to nonlinear stability. This does not make global physical
SEA0 left inclusion a prerequisite for testing one independently validated
legal point. It does prohibit self-asserting provider closure or treating an
algebra fixture, replay or unrelated bounded sequence as that point.

The operation regressions verify full prediction/reset/lift shift transport
and the information identity for each supplied S covariance. They are not the
required feasibility experiment and authorize no new proof-bound module.
The strongest critic objection remains: additional coordinate identities do
not demonstrate a complete-word margin. The next executable result must be a
legal same-history word and its signed margin, not another claimed radius.

## Next falsifiable experiment

The controlling object is `rho_W=V_after(F_W(x))/V_before(x)` on the actual
complete H18/A21 word, not positivity of an isolated lemma.

First choose between three genuinely different treatments: (1) retain shipping
coordinates and accumulate the exact signed full residual/Joseph/reset form;
(2) evaluate the corrected full-Phi nonlinear storage with all mixed
transport and uniform metric-comparison obligations; (3) eliminate the coupled a_w direction only in
the full-word Schur form, retaining all actual-R_S S events. Start with (1) as
the unchanged-shipping reference against which the others must be checked.

Before further interval refinement, run a non-promoting high-precision
complete-word feasibility test on legal same-history SEA3 realizations. Report
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

## Continuation boundary

Continue the complete-word feasibility experiment above from the reviewed
branch/checkpoint, preserving newer main changes. Keep P3 frozen and P4/P5 open.
Do not merge without explicit authorization.
`tests/validation/test_ou3_p4_coordinate_transport_algebra.py` locks the exact
full-shift, reset, and correlated-precision corrections. The full source-
correlated H18/A21 transport and storage comparison are still the decisive
obligations; green algebra/architecture tests do not close them.
