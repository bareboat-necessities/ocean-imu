# OU-III proof research state

## Current hypothesis: bounded-bias practical motion stability

The user-authorized primary theorem now asks for bounded accelerometer-bias
error and practical stability of the other 18 errors, not full A21 nonlinear
contraction or bias convergence. Both physical SEA3+ branches and both actual
filter modes remain required; A21 still executes its full 21-state P/K/reset
history with every actual R_S. The performance storage is 18-dimensional;
this is not a reduced estimator or a copied held-mode H18 certificate.

BIAS0/1 gives a conditional same-history true-bias bound B_true; projection
gives the estimate ball and hence ||e_b||<=B_true+0.4. The weaker target uses
the closed 0.4 ball, including projection boundaries, instead of assuming
the old 0.35 interior. The legacy domain and full-state flags remain intact;
their P3/source coverage is not automatically extended to the new domain.
BIAS2 is optional gain sharpening, valid only with a certified graph sector.
Bias corrections remain internal coupling, not renamed sensor noise.

The controlling new inequality is
`W_H,N <= rho_H W_H,0 + gamma_b D_b + gamma_s D_s + gamma_n D_n`, rho_H<1,
with finite every-event prefix gain, chart/source retention and hybrid entry.
The gain/floor master uses the same full nonlinear graph and actual forcing.
Zero sensor noise does not remove ocean/OU model mismatch. A nonzero residual
floor is allowed and expected; the observed 2--3% is not yet a certified
uniform constant or an identified output normalization.

Current evidence: the conditional compactness and geometric-series/retention
lemmas are in LaTeX; exact composition arithmetic and an outward full-graph
motion-master assembler are tested. Synthetic algebra witnesses are not
Ocean-IMU motion-gain certificates. New P4_MOTION_PASS/P5_MOTION_MAY_START
remain false independently of the stronger full-state flags.

Current limiter and next falsifiable experiment: reconstruct exact joint
error/source/forcing subevents, then evaluate the motion master with explicit
channel gains on the same complete word before any uniform enclosure search.
Report the resulting motion-accuracy floor and test that it lies within the
retained chart and requested physical performance tolerances. Preserve P3's
delta=1e-18 while resolving projected-domain source coverage separately.
No new metric grid, packetwise scalar budgets, longer-window search or
filter changes are authorized. P5-motion cannot start from compactness alone.

Local build limitation: `make all` stops at `KalmanQMEKF.h:30` because
Eigen/Dense is absent; the include path is unchanged. The theorem smoke
harness initially loaded an unavailable, unused siunitx package; it now loads
only packages used by its actual inputs. Neither issue is a theorem failure;
GitHub CI supplies the independent proof and LaTeX rebuild results.

## Retained stronger-target source-capture evidence

The source-connected **actual finite** endpoint capture passes CI run
`34157292548` on code head `456fce0f`. One shipping observer owns the complete
power-on history, nominal state, P, resets, tuner commits and actual R_S.
The new conditional Stokes capture uses zero true bias/driver and makes no
state/covariance intervention. It is not the original retained witness.

| Fixed window | Sample indices | Source parity defect | Actual forced VN/V0 | Actual sample-prefix maximum |
| --- | --- | --- | --- | --- |
| H18 | 6600--7199 | 1.56e-15 | 7.97176522677814 | 7.97176522677814 |
| A21 | 228600--229199 | 1.68e-14 | 0.496821403940055 | 1.46706257365206 |

Both have 600 IMU calls and 75 magnetic calls; H18/A21 retain 137/108 due
S calls. The checked Live/mode/magnetic-transition, accepted-accelerometer,
physical acceleration/rate, bias-interior and Cayley-chart predicates pass.
This is point evidence, not full Normal-Live/PE or source-uniform admission.
Endpoint storage is evaluated by exact rational elimination with 80-digit
decimal output from the float shipping coordinates; call margins telescope
to reported roundoff. The capture artifact is `10031432066`, SHA256
`ce787b17bcf700730f6905c87f2e5a0b303dd18e20da697de3b63a27afc7ff46`.

The failed identification is **zero sensor/bias noise => homogeneous P4
trajectory**. The same-root physical prediction defect
`u_l[k]=x_true[k+1]-F_LL[k]x_true[k]` is nonzero: the largest acceleration
coordinate increments are 0.018034304915645463 (H18) and
0.022571626343973056 (A21). Classification: proof-map attachment gap, not
a failed contraction theorem or interval enclosure. Neither raw ratio proves
or falsifies the homogeneous P4 storage; do not optimize either one.
BIAS1 does not remove ocean/OU model forcing. Assembled-sensor qualification
is not a gate for this conditional work.

Critic: abandon using an actual zero-noise simulation as the unforced master;
its nonzero source-model input makes the ratio answer a different question.
The current limiter is the exact joint error/source/forcing transition,
including subevent nominal/P/reset attachment, not an eigenvalue margin.
Next falsifiable experiment: construct that transition and require it to
reproduce every captured shipping subevent at the recorded inputs. Then
test whether its zero-input section preserves the claimed source attachment
before composing the fixed complete-history homogeneous endpoint. If it
does not, the qualitatively different alternative is the full forced
window ISS master with explicit disturbance coordinates/cross terms, not
deleting inputs and retaining the old covariance cells. BIAS2 must enter
the correctly attached nonlinear graph, not a post-hoc PSD or sampled ratio.
No interval/metric search is justified yet; P4/P5 remain open/blocked.

All 27 focused tests and the retained-byte audit pass in CI. The first capture
run hit a NumPy-bool JSON emission defect, repaired without changing the
history or mathematical gates and covered by a complete-report regression.
Locally `make all` remains infrastructure-blocked at `KalmanQMEKF.h:30`,
missing `Eigen/Dense`; the local lint executable is absent. CI installs those
dependencies; the source probe compiles with warnings as errors.

The retained-witness CI run `34148309564` rejected a different payload SHA256
(`c6b99250...` instead of `db38b812...` for H18). The regenerated A21 scan also
selected 1195.01025390625--1198.01318359375 s instead of the retained
1143.96044921875--1146.96337890625 s. This is an evidence-reproduction failure,
not a BIAS theorem failure. Regeneration is not content-addressed evidence.
The next check must consume the retained bytes, retaining the identity gate.

The finite attachment hypothesis that deleting covariance resets while keeping
F/Q/H unchanged is a coordinate change has failed. For a frame T,
the required maps are F'=T_next F T^-1, Q'=T_next Q T_next^T,
H'=H T^-1 and the exact finite map f'=T_next f(T^-1 z). Omitting these
transports invalidates that justification; it does not invalidate the P3
congruence identity or establish a canonical P4 counterexample. The next
falsifiable check compares these required transforms with the retained
diagnostic at every operation. No metric or interval search is authorized by
the diagnostic's expanding ratios.

Canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`, now with the explicitly authorized response union `LINEAR_VESSEL` union `STOKES_WAVE_FOLLOWING`. The legacy linear domain is unchanged. Conditional P3 uses the same response-independent Normal-Live implication at `delta=1e-18`, with explicit H18/A21 coverage for both branches. This does not certify physical word admission. P4 is **OPEN** and P5 is **BLOCKED**; H18->A21 remains a separate rectangular hybrid event, with zero lever arm and dormant-transparent vibration.

The only production/proof-domain change made on PR #496 is the user-authorized accelerometer-bias projection-radius tightening from `0.5` to `0.4 m/s^2`. The declared startup/handoff accelerometer-bias error envelope is also `0.4 m/s^2`, while the Normal-Live active-bias interior bound is `0.35 m/s^2`, preserving a `0.05 m/s^2` projection margin. No other filter tuning, quality gate, source-language parameter, or P3 mathematics was changed.

The retained stronger paper target uses the finite full-state source-indexed quadratic storage

`V(e,zeta)=e^T M(zeta)e`

on one complete same-history SEA3 word, with a strict finite endpoint inequality, a finite every-prefix gain, and every-prefix chart/source-domain retention. A replay, point experiment, fitted endpoint metric, reduced-state coercivity test, or alternate source language cannot promote P4.

## Exact structural work retained on this checkpoint

The branch retains the theorem-facing algebra and execution machinery needed by a continuation PR:

- exact finite accelerometer coordinate shift with original shipping H/P/K/S and `H0 != H_u` at finite attitude error;
- exact prediction transport with the literal full `F E_aw`, retaining v/p/S/a_w rows;
- exact Joseph/reset signed-information identity
  `Delta V = -I_y + E_eta + X_reset + E_reset`;
- exact deployed Cayley reset transport and A21 `0.4 m/s^2` bias projection with generalized Jacobian;
- every due S=0 event with actual applied anisotropic SpectralMSE `R_S` retained; S events have `eta=0` exactly and therefore contribute favorable information;
- every valid accelerometer update, applicable vector update, full Q, covariance floors and immediate resets retained;
- separate H18->A21 hybrid lift;
- outward interval/differential AD and generalized mean-value machinery for the exact nonlinear physical map;
- complete-word accelerometer covariance channel and exact homogeneous Cayley residual-sector factorization;
- the existing finite-`tau_b` A21 detectability module tied to complete SEA3, which closes the paper-level finite-bias detectability/UES hypothesis but explicitly does **not** close the canonical full 21x21 implementation-word P4 inequality.

Do not replace these with selected-S words, independent tuner/R_S boxes, independent per-sample source boxes, packet-count nonlinear budgets, scalar correction radii, inverse-metric-floor arguments, state elimination, or replay-derived source families.

## Authoritative single-observer linear evidence

The corrected single shipping observer owns one source/tuner/Riccati history and retains every operation. On the genuine PM+Stokes Hs=1.5 m history, the pre-0.4 legal 3-second point words were:

- H18: `rho_linear=0.9998658024147671`, 600 predictions, 600 accelerometer updates, 137 actual-R_S S updates, 75 vector updates;
- A21: `rho_linear=0.9958536807113242`, 600 predictions, 600 accelerometer updates, 108 actual-R_S S updates, 75 vector updates.

Same-observer event ledgers reproduce those ratios within a few `1e-6` and telescope to roundoff. The old duplicate observer that staged tuner-commit data at the wrong boundary is retired; its older A21 rho values are stale.

On the 0.4 head, the 3-second physical reset-normalized diagnostic keeps strict zero-state parity. H18 remains contracting over the retained tested scales. A21 first crosses one at scale `8.0`. The bias-admissibility audit below excludes the approximately `1.086` boundary case from the 0.35 Normal-Live estimated-bias interior, but retains expanding +/-8 and +/-16 cases with inactive projection. The 0.4 clamp does not remove that interior diagnostic obstruction.

## 6-second and 9-second falsification result: longer-window route stopped

The source-contiguous physical long-window diagnostic completed successfully in GitHub Actions run `34083473117` (`ou3-p4-physical-long-window-feasibility`), artifact `10004621877`, artifact SHA256 `8b3c1e3b83e6774895cf20770764609afd38e044131204c568d49d7f3544aa19`.

It uses the same single shipping observer, the canonical source label, actual applied R_S, all accelerometer/S/vector events, strict zero-state parity, and no source/domain/filter substitution. It is explicitly non-promoting.

6-second result:

- H18 linear rho `0.9917547078907433`; worst retained physical rho `0.993720034085038`;
- A21 linear rho `0.991109058400818`;
- A21 first finite-scale crossing `7.5`;
- A21 worst retained physical rho `1.1357228916109403`;
- A21 worst prefix ratio `1.1358285446129752`.

9-second result:

- H18 linear rho `0.9786791982272148`; worst retained physical rho `0.978951454184692`;
- A21 linear rho `0.9868275426222095`;
- A21 first finite-scale crossing `16.0`;
- A21 worst retained physical rho `1.143708327664902`;
- A21 worst prefix ratio `1.1653466626592186`.

Therefore the proposed 3->6->9 second extension does **not** repair the finite A21 mechanism. Per the research protocol, stop the longer-window route here. Do not spend another iteration optimizing window length or point storage around this replay.

## A21 mechanism and dead ends

The finite A21 obstruction is primarily finite accelerometer curvature in a coupled attitude / latent-acceleration / accelerometer-bias cancellation direction. The exact lever-arm-off residual is

`y = (E-I) f_hat + E R_hat delta_a_w + delta_b_a`,

with first-order row

`H e = [c]_x f_hat + R_hat delta_a_w + delta_b_a`.

These first-order terms can nearly cancel while second-order attitude curvature remains. The A21 bias projection is not active at the problematic scales, and reset terms are small; tightening reset bounds or eliminating a_w attacks the wrong mechanism.

Retired / forbidden rescue routes include:

- estimator-pair shadow promoted as theorem map;
- raw or full-Phi endpoint optimization after their A21 finite-scale failures;
- arbitrary single-map converse-metric fitting;
- extending the failed physical replay to still longer windows merely to seek a green point;
- reduced-state/Schur certificates or a_w elimination;
- selected-S words or independent R_S/tuner schedules;
- scalar Lipschitz, correction-radius, inverse-metric-floor, or packet-count-times-worst-remainder bounds;
- further proof-driven filter/domain tightening beyond the authorized 0.4 change.

## Retained finite-`tau_b` premise

`tools/stability/ou3_sea3_a21_detectability_completion.py` establishes the
paper-level finite-bias detectability/UES hypothesis from:

- complete-SEA3 H18 contraction;
- exact finite residual-bias Gauss-Markov decay;
- bounded full-state H18<->b_a coupling on the compact word;
- full A21 process UCC;
- no alternate estimator and no state elimination.

Its stronger implementation-word flags remain deliberately false. The joint
master consumes this result only through the validated canonical P3 chain; it
does not promote the comparison observer or replace the full A21 matrix test.

## Retained stronger full-state hypothesis

Use the exact complete-word endpoint identity and the full 21-state finite-`tau_b`
P3/detectability result to build one correlated nonlinear graph-sector master.
For `z=[x;w_W]`, the controlling matrix is
`L_W=[[-D_W,M_W^T J_N B_W],[B_W^T J_N M_W,B_W^T J_N B_W]]`.
Admissible graph sectors `z^T Pi_j z>=0` enter only through the full
S-procedure test `-(L_W+sum lambda_j Pi_j)>0`; the same construction is
required at every prefix with finite gain, not prefix contraction. Set
`D_k(Gamma)=Gamma J_0-M_k^T J_k M_k`, with finite `Gamma>=1`, in the
same augmented builder. Every-prefix domain/chart retention remains separate.

## Evidence and current limiter

P3 was recomputed with the deployed `0.4 m/s^2` accelerometer-bias projection
limit: H18 and A21 retain `delta=1e-18`, the first active A21 bias full-matrix
margin is `1.2499987189052501e-9`, and the H18 worst interval LDLT pivot is
`4.987499868870966e-14`. `P3_DEPLOYMENT_PASS` remains false for the
physical-language inclusion and assembled-sensor BIAS0 qualification obligations.

The canonical 6 s A21 payload has a strict small-error margin
(`rho=0.9911176` at scale `0.125`) but crosses one at scale `7.5` and
reaches `1.1357229` at the retained-domain boundary. The finite-`tau_b`
detectability rerun passes with bias energy gap `1.1992803e-3` and
asymptotic A21 gap `1e-18`. Thus linear bias decay is not the limiter;
the limiter is the correlated attitude/`a_w`/`b_a` accelerometer curvature.

## Failed approaches / DEAD_ENDS

On the same canonical payload, diagonal bias precision, fixed-frame and
source-framed `a_w`/`b_a` cross penalties, and separate latent diagonal
energies all failed after their allowed refinement. The closest result was
`rho_linear=0.996797290`, `rho_finite=1.00009035805` at
`(beta_aw,beta_b)=(500,310000)`. Do not resume metric-grid tuning.

## Retained facts, alternatives, and next experiment

Retain the complete same-history source, frozen P3 `delta=1e-18`, full H18/A21
states and cross terms, exact Cayley/reset residuals, all valid accelerometer
and vector events, and every due S event with actual applied `R_S` inside its
suffix. Packetwise radii, state elimination, replay fitting, and further
domain/filter changes remain forbidden.

The remaining alternatives are a dense source-structured storage LMI, a
path-dependent joint storage, or direct falsification of a source-uniform
inner funnel. First obtain CI evidence for the connected execution fixture and
resolve the complete-source witness. Before any new sector search, reconcile
the retained A21 counterexample with the exact proposed source/error domain,
physical map and storage. If it is admitted, valid graph sectors cannot change
its `rho>1`. Conversely, failure of an S-procedure relaxation alone does not
falsify the storage. Only a demonstrated feasible full-word formulation
justifies source-uniform enclosure work.

## Complete-source obligation still open

Direct physical-generator admission to the **linear branch** has a structural failure. The retained
input names the v1.1.3 PM--Stokes surface-particle generator, pinned at
`oceanography-waves-lib` commit `c5ddd8ddba6e062bb131d92efcd672dfa189455a`.
It uses 128 logarithmic fundamental frequencies in [0.02,0.8] Hz, order 3,
common seeded phases/directions, and attitude from the same advected surface
slopes. This is an actual correlated generator, not 600 independent boxes.
It is not itself the declared continuum linear vessel-response model.
Relative to its fundamental-only sea, nonzero higher harmonics cannot be
created by a linear response. Relative to its full Stokes elevation, the
surface-particle displacement response is h=(i cos(theta),i sin(theta),1)
(up to phase convention), so ||h||^2=2. At the highest third harmonic,
f=3*(4/5)=12/5 Hz, the largest allowed SEA3 squared gain is
[4*( (6/5)/(12/5) )^2]^2=1. The exact failed inequality is 2<=1.

This is a model/premise incompatibility, not an interval or A21 theorem
failure. It invalidates identifying Stokes with the linear RAO, not the
user-authorized union `G_SEA3+ = G_vessel_linear union G_wave_following_Stokes`.
That explicit new model branch makes the old rejection inapplicable to the
union; the exact failed linear inequality remains a regression test.
The new root retains all phases/directions and bound-harmonic coefficients,
not just seed 42 or one frequency grid. Fixed finite N is a Stokes source
coordinate; it does not replace or certify the legacy continuum branch.
Per-root hard bounds follow from the actual coefficients; linear rolloff
moments and surface-period identities remain scoped to the linear branch.

Current limiter: response-model admission is weaker than full retained-word
membership. The source clock/root, same-history frontend/tuner, BIAS0/1
physical bias, actual Normal-Live bounds and finite reset/storage attachment
are still required. The pinned generator's centered-difference gyro also
needs an explicit discretization-defect enclosure; it is not exact Cdot.
No bias feedback is silently moved into ISS. BIAS2 must be certified jointly
on each full corrected-error history, not inferred from a Stokes label.

Next falsifiable experiment: CI checks all 128 actual root amplitudes,
steepness and dependent harmonic coefficients against the Stokes branch,
retains rejection by the linear branch, and rebuilds conditional P3 plus
endpoint/prefix/retention P4 obligations for **both** branches/modes.
Then attach the complete source and finite storage to the retained word.
A Stokes-only certificate cannot promote the enlarged union.

The upstream complete-source contract still does not materialize the correlated finite-window SEA3 realization as an executable outward source family. Parameter compactness, RAO/moment envelopes, hard pathwise acceleration/body-rate caps, frontend parity, and adaptive-state rate/jump bounds do not replace a same-history transition for the correlated source state. P4 may not substitute independent per-sample boxes, a replay record, or a finite harmonic/grid surrogate.

The full-state master now identifies the needed quantities: same-history prefix
state selectors, stacked nonlinear residual selectors, reset/boundary terms,
and the A21 projection graph. A cell/lineage identifier is not a hard SEA3
membership or inclusion proof. The translational response envelope does not
supply the missing joint rotational realization or hard correlated driver.
No additional linear-vessel rotational envelope or coherence budget is
adopted. The added Stokes branch has its own joint slope-derived rotation.

The A21 projection also needs the same-history true temperature-centered
residual bias. At zero bias error its output is
`b_true-project_0.4(b_true)`. The homogeneous map must prove this is zero;
an OU covariance or the estimated active-state cap does not provide that
true-bias premise. Source forcing and the separate H->A event remain explicit.

## Bias-premise audit and current proof plan

The fixed-witness audit uses the archived directions/scales and SHA256-pinned
3 s H18/A21 payloads. It assumes the same zero true-bias root, zero unmatched
deterministic offset and zero homogeneous driving as the retained evaluator;
it does not wait for assembled-sensor qualification. BIAS1 admits this true
bias algebra. It checks the initial state and every event against the existing
0.35 Normal-Live estimated-bias interior as well as the old error-domain checks.

The point result is:

| A21 scale | Endpoint rho | Maximum estimated-bias norm | Sufficient BIAS2 point ratio, Xi=V0 | Checked bias/domain result |
| --- | --- | --- | --- | --- |
| -8 | 1.0022511802 | 0.1060604400 | 0.0067782449 | retained |
| -16 | 1.0213229467 | 0.2121209584 | 0.0068726407 | retained |
| -30.1714758181 | 1.0860153089 | 0.4000000000 | 0.0072294990 | initial estimate outside 0.35 interior |

The positive +8/+16 cases also expand within the checked bounds. All four
interior expanders have inactive projection and 600 predictions, 600
accelerometer events, 108 actual-R_S S events and 75 vector events. Corrected
bias-error recurrences agree with the retained finite evaluator; at scale -16
the corrected error differs from a free GM path by up to 0.0001764025 m/s^2.
The negative boundary case has two active projections and was never an
inactive-projection witness. H18's worst retained ratio is 0.9998370258.

BIAS2 is evaluated jointly: a and d are the full stacked whitened nonlinear
attitude/acceleration and corrected bias-error channels. At scale -16,
kappa_point=0.9999717099, so separation is positive but cancellation remains
strong. These point ratios (and their sampled minimum) are not uniform
constants. Positive separation coexists with positive Delta V in the retained
storage. For every sector valid at this point and lambda>=0,
`Delta V + sum lambda_j*z^T Pi_j*z >= Delta V > 0`. Thus an admitted expansion
cannot be repaired by stronger use of valid BIAS0/1/2 sectors in that storage.

The source membership decision is **undetermined**. OU3PHY1 stores derived
P/H/R/F/Q, resets and timing; it lacks the common hard-driver/phase/response
realization, frontend/tuner root, raw measurements, nominal-state history and
physical-bias/forcing decomposition. The original SEA3 hard-driver set and
joint output map also lack an executable membership characterization. This
is missing mathematical source data, independently of assembled-sensor
qualification. A replay label, norm cap or zero bias root cannot fill it.

The proposed finite gauge attachment is **rejected**. Set T0=I and propagate
T_next=T G^-1 at each captured reset. At event 2 the next H needs a nonzero
transport: Frobenius difference 0.00897303817 for H18 and 0.0000609319492 for
A21. Maximum required H changes over the retained words are 0.4790199261 and
0.0104868913 respectively. F and Q require transport too. The reset-deleted
diagnostic keeps the original matrices, so it is a rebuilt Riccati history,
not the captured word in new coordinates. This invalidates using its expanding
rho as a canonical counterexample; it does not invalidate the conditional P3
congruence theorem, finite Cayley algebra or BIAS0/1/2.

The finite storage identity is executable at every retained measurement:
q=e-Ky, r=e_plus-Gq, b=G^-1 r gives
V_plus=q^T Pj^-1 q+2 q^T Pj^-1 b+b^T Pj^-1 b.
It includes projection defects and retains the signed cross term. This
one-event attachment does not authorize reuse of subsequent reset-deleted
covariances. The audit also implements exact finite frame pullback and tests
dense H18/A21 Joseph, prediction/noise and nonlinear energy identities.
An initial numpy/list adapter error was an implementation defect (ambiguous
array truth value in the retained signed-ledger helper); converting its inputs
to the helper's sequence representation resolves it without changing algebra.
The prefix CSV retains each checked state and residual.
No filter, domain, direction, scale, metric or window-length search is added.

Critic decision: abandon the claimed reset-deletion gauge equivalence.
The next falsifiable experiment requires a common source witness and either
(1) an actual finite shipping covariance/nominal history with each correction's
reset, or (2) its fully transported finite representative with all transformed
F/Q/H and floors. These are attachments, not metric searches. A third option
is to prove a conservative outer source inclusion and certify the whole
nonlinear graph there; success would suffice without inverting a replay.
Source-only names and independent sample boxes cannot establish that inclusion.
Only after this attachment gate may a full-word BIAS2/storage feasibility
decision justify interval certification or abandoning the storage itself.

The BIAS0/1/2 package in `w3d-mems-bias-preconditions.tex-part` is conditional.
BIAS0 separates qualified sensor physics from the filter's 5000 s setting and
0.4 projection; assembled-unit tau/sigma, true-root and mismatch bounds remain
unassigned. BIAS1 retains one true-bias root and GM parameter through every
prefix. Main's nonlinear lineage and correlated window-cell data model are
retained; the nonlinear consumer now checks the GM graph, not just matching IDs.

The proposed substitution `delta_b_j=exp(-t_j/tau_b)*delta_b_0` through a
corrected word is invalid: every measurement changes the error by `-K_b*y`,
with the projection defect added afterward. This is a mathematical premise
failure, not a new filter instability. It invalidates treating corrected bias
error as a free slow source or charging its feedback to exogenous ISS forcing.
It does not invalidate true-bias GM transport or the existing finite-tau P3
detectability argument. BIAS2 is therefore a conditional full nonlinear graph
separation lemma; its uniform positive `mu_sep` is still missing.

The controlling contribution is `lambda_sep*(C_y^T W C_y-mu_sep X)` in the full
augmented master. There is no certified numerical contribution yet, so this
package authorizes no new interval or metric search. First check the retained
expanding A21 witness against the exact BIAS1 source/error graph. If admitted,
valid sectors cannot rescue that storage. The alternatives remain dense
source-structured storage, path-dependent storage, or direct falsification.
The stage-by-stage BIAS0-to-P5 dependencies are in `tools/stability/README.md`.
Rebuild conditional P3 at delta=1e-18 and canonical P4 with these explicit
premises; keep physical qualification, source cover and P4/P5 closure separate.

## Connected execution experiment

The rejected `P0_H=2 I_18`, `P0_A=2 I_21` two-prefix fixture loses innovation
pivot 1 at `[-31.51349023865584,33.33589453040376]` in the unchanged outward
Joseph inverse. This is an enclosure failure of that arbitrary fixture, not
a filter, conditional-P3, or P4-contraction result. Do not regularize the pivot,
drop a prefix, or select a favorable frontend successor.

The single motivated refinement uses Live-structured covariance: the same
committed stationary `a_w` variance, tilt/yaw handoff, and shipping `b_g/v/p/S`
seeds. A21 appends the shipping bias-floor variance for the fixture only.
It does **not** execute pre-release H18 history or certify A21 entry reachability.
The controlling diagnostic is completion of both prefixes with unchanged
Joseph arithmetic and uninterrupted event-local covariance ancestry.

CI produces `ou3-same-history-prefix-execution`: per-prefix H18/A21 covariance
diagonals and interval widths, event order, and captured actual `R_S` diagonals.
Its two-sample scope and unclosed source/entry/P4 flags are mandatory. The
canonical P3->P4 dependency chain is unchanged. The critic boundary is that
this repairs execution evidence, not the known A21 endpoint inequality.
If the Live-structured fixture also fails, stop this refinement; the distinct
alternatives are a covariance representation change, direct source-witness
construction, or abandonment of this execution representation after audit.

## Current CI / evidence status

The MEMS-bias integration rebuild exposed an inherited CI lint defect:
`ou3_sea3_correlated_window_cells.py` imported unused `dataclasses.replace`
(ruff F401). Removing the import repairs that implementation defect without
changing the source model or quality gate; the next check is the unchanged CI
python job. The focused five bias algebra/negative-contract tests pass locally.

The expanded same-history CI suite passed 46 of 49 tests, including the bias
and nonlinear lineage tests, but exposed the inherited source-cell splitter's
endpoint defect. Reapplying `outward_bounds` to stored endpoints changes a
zero lower hull endpoint to `-1e-323` after two roundings and creates a nonzero
overlap at the cut. This is an intersection-construction defect, not failure
of SEA3 or a contraction inequality. Copying existing parent endpoints and
one exact binary64 cut restores the exhaustive binary partition; exactly
representable smoke-domain endpoints need no arithmetic rounding. The cover
validator and outward arithmetic stay unchanged. This repairs structural
partitioning only; the hard joint SEA3 oracle remains the limiting quantity.
The next check is the same CI split tests, including repeated non-midpoint cuts.

On head `6ed27d24`, GitHub CI passes 73 source-foundation tests, 50 connected
execution tests, and the rebuilt conditional H18/A21 P3 at delta=1e-18 with
46 contract tests. All quality checks pass. The canonical P4 artifact builds,
but its expanded discovery suite passes 104 of 105 tests: the mocked P3 input
in the joint-sector test lacks the newly required `mems_bias_preconditions`
field and raises `KeyError`. This is an incomplete test-fixture integration,
not a nonlinear inequality failure. Add the actual bias contract to that
fixture and check its propagation; retain the production dependency check.
The next check is that focused test followed by the same canonical CI chain.

`make all` is locally blocked before compilation because
this environment has neither `/usr/include/eigen3/Eigen/Dense` nor a vendored
Eigen tree; the build command and include policy were not changed.

The SEA3+ local regression initially misplaced the cache-isolation assertion
inside the new branch test (`NameError: cached`). This was a test-edit defect,
not a failed matrix inequality: the source and P3/P4 builders completed.
Restore the assertion's original scope and rerun that same suite; no proof
constant or admission gate changes. GitHub CI remains the C++ build check
while local `make all` fails at `KalmanQMEKF.h:30`, missing `Eigen/Dense`.

SEA3+ CI run `34153666156` passes the C++ public-API audit and 21 tests:
all 128 atoms satisfy the Stokes model, maximum component steepness
0.02503536704829765 < 0.2. Response admission is ADMITTED on the new branch;
full word membership stays UNDETERMINED. The coupled proof run `34153666115`
stopped before P3: the package-integrity test requires the response-domain
JSON reference to include `tools/stability/`. This is a path-metadata
integration failure, not a matrix/source inequality failure. Correct the
reference, preserve the integrity gate and rerun the same certificate chain.

`ou-validation` is red for a known evidence-provenance reason, not because its numerical unit-test body found a new filter failure: `tools/ou_evidence_contract.py --auto` reports replay dependencies changed relative to committed validation/robustness provenance, including the OU-III filter and WavePeriodEstimator dependencies. Genuine validation/robustness evidence regeneration is therefore still required before a later proof PR is declared final/ready. Do not hand-edit provenance hashes.

No claim is made that P4 or P5 is complete.

## Promotion boundary

The master and graph-sector machinery are non-promoting. P4 remains open until
the source-uniform endpoint and every-prefix augmented LDLT plus domain
retention close on the same complete SEA3 history. P5 remains blocked until
strict canonical P4 contraction closes.
