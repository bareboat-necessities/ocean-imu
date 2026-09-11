# Parallel OU-III whole-word dissipativity proof

> **Continuation:** after PR #517, start a NEW PR from latest `main` and read
> [`ou3-alt-contraction-handover.md`](ou3-alt-contraction-handover.md) first.
> The existing P2/P3/P4/P5 route remains independently continuable through
> [`ou3-brmm-main-handover.md`](ou3-brmm-main-handover.md).

## Status and coexistence

This is the independent **ALT** proof route. Its current research question is:
can an inverse-free, full joint24, finite-increment/IQC word formulation support
a shorter Live stability argument for the unchanged shipping implementation?

The existing P2/P3/P4/P5 route remains authoritative and independently
continuable through [the canonical handover](ou3-brmm-main-handover.md).
Nothing in ALT changes its algorithms, theorem assumptions, proof tools,
workflow gates, frozen P3 threshold `1e-18`, or next experiments. Shared source
contracts are consumed read-only. The [research ledger](ou3-proof-research-state.md)
contains separate ALT failure analysis. Neither route's tests substitute for
the other route's certificates.

**Current result:** exact inverse-free algebra, shared-bias lifts, a projection
IQC, structural anti-shortcut regressions, and a captured frozen-map diagnostic
are executable. A source-uniform shipping joint24 master has NOT been solved.
`ALT_LIVE_PASS`, `ALT_STARTUP_PASS`, and `ALT_END_TO_END_PASS` are false.
Existing P4/P5 are not promoted. A green exploratory workflow means successful
execution of the experiments, not proof completion.

## 1. The right state and target

Use true-minus-estimate errors and retain

```
z = (e_theta, e_bg, e_v, e_p, e_S, e_aw, e_ba, beta_true) in R^24.
```

The physical/source, covariance, frontend, candidate/active tuner, scheduler,
clock and hybrid mode states form `xi`. This indexing is NOT permission to
freeze endogenous dependence: an incremental comparison must retain the
actual change of these quantities, or prove an alternative same-trajectory
storage argument that legitimately does not differentiate them. Add internal
IQC-filter storage when a valid hard dynamic IQC requires it.

Keep corrected COMPLETE-BRMM exactly as specified by the canonical handover,
including the 8.8 m/s^2 acceleration cap, 1100 m*s all-time centered-primitive
cap, one-time Live S origin, zero lever arm and actual anisotropic R_S. Do not
replace physical source bounds with convenient estimator-error radii. Maintain
the same-signal WPE/bandpass/statistic relation through f, sigma, tau, T_S and
R_S, and preserve every staged commit and scheduler decision.

### A necessary correction to strict full-state contraction

H18 holds the estimated accelerometer bias. BIAS2 admits constant physical
bias. For a common SPD metric, any linearized augmented map with `A v=v!=0`
violates `A^T M A <= rho M`, `rho<1`, since the quadratic gap on v equals
`(1-rho) v^T M v > 0`. A uniformly bounded source-dependent metric cannot
manufacture indefinite decay of a persistent nonzero coordinate either.
This rules out the over-strong homogeneous target, not motion ISS.

Do not fix this by discarding A21 motion/bias cross-information. Instead retain
coercive joint storage and charge independently justified bounded bias/source
ports through a supply. The desired complete-word inequality is

```
m ||z||^2 <= V(z,xi) <= Mbar ||z||^2,
V(z_next,xi_next) <= rho V(z,xi) + w^T Gamma w + c^T Beta c,
0 < rho < 1,  Gamma >= 0, Beta >= 0.
```

Here w contains declared disturbance/model/roundoff ports, and c contains
explicitly bounded neutral bias or physical-reference ports. Every bound must
come from the admitted history, not from the desired conclusion. In particular,
`||e_ba|| <= R_projection + B_true` needs the actual estimate projection/hold
invariant, including entry, temperature centering and finite precision.
Charging all unknown motion errors as bounded c would be circular and is
forbidden. Keep cross terms inside V; do not reset storage at word boundaries.

If the total word supply is uniformly at most C, iteration gives

```
V_k <= rho^k V_0 + C (1-rho^k)/(1-rho),
||z_k||^2 <= [rho^k V_0 + C/(1-rho)]/m.
```

This is a regional practical bound with a separately bounded bias, not a claim
that physical wave motion or bias tends to zero. To obtain motion ISS with a
useful small ultimate bound, report the actual supply coefficients and motion
coercivity; a merely finite but enormous C is not basin closure.

An incremental theorem alone does not establish error to physical truth:
truth need not be a solution of the estimator, notably at the S=0 pseudo
measurement. Retain the reference-map defect as same-history physical forcing.

## 2. An inverse-free measurement graph, including the H18 mask

Let N denote the **actual** `PCt` gain numerator after shipping masks, and S the
actual innovation matrix after its declared numerical handling. In real
arithmetic, the innovation correction is described by

```
S q = r,       correction = N q.
```

No independent interval enclosure of S^-1 or K is needed in the final graph.
The shipping code already uses LDLT solves. This is a change of proof
representation, not an algorithm replacement.

For a fixed-coefficient true-minus-estimate perturbation, r=H_residual e+w,
so e_next=e-Nq. The descriptor equality and selectors in `core.py` use lifted
variables chi=(e,w,q) and E=[-H_residual,-I,S]. N is NOT silently replaced with
P H_residual^T: the H18 accelerometer numerator masks held-bias terms, while S
retains bias uncertainty. Generic full-state information addition therefore
requires a separate effective-measurement/invariant bridge before it applies.
A21 identities also need the actual reset/projection/numerical branch bridge.

### Retain endogenous gains without linearizing the inverse

For two finite states with S0 q0=r0 and S1 q1=r1, exactly

```
S1 delta_q + delta_S q0 = delta_r,
delta_correction = N1 delta_q + delta_N q0.
```

`measurement_increment_lift` implements these equations with row-major
vectorized delta_N and delta_S. Nonzero nominal residuals are retained. This
is stronger algebra than freezing K or assuming r0=0 and avoids differentiability
at this step. Source/covariance/tuner constraints must still bind both N and S
to their actual states. The derivative limit contains the otherwise omitted
`dS*q` and `dN*q` terms. Floating-point solve/correction residuals need additional
supply ports; the real-arithmetic identity alone is not their enclosure.

Keep the implemented Joseph polynomial

```
P_next = P - K N^T - N K^T + K S K^T
```

and the actual covariance reset, PSD/floor handling and rejection branches.
An information-form identity may simplify a proved subcase; it cannot stand
in for these maps without equivalence. Structural S positivity alone supplies
neither a tight gain bound nor a complete-word contraction margin.

## 3. Same-history bias and nonlinear graphs

For EACH separately admitted BIAS0/BIAS1/BIAS2 history,

```
beta_next = phi_true beta + w_b,
e_ba_next = phi_hat e_ba + (phi_true-phi_hat) beta + w_b.
```

`bias_prediction_lift` retains the mismatch in the joint matrix and puts the
SAME driver column in e_ba and beta. H18 has phi_hat=1; A21 uses its actual
OU factor. Do not make two independent copies of w_b or independent per-sample
truth slots. Family definition validation is not startup or hardware admission.

The actual radial projection graph is

```
e_ba_next = beta - Pi_R(beta-e_ba_pre).
```

For increments u=delta_beta-delta_e_pre and y=delta_beta-delta_e_next,
firm nonexpansiveness gives `y^T(u-y)>=0`. The implemented joint IQC retains
beta on both sides and all relevant cross terms. This finite-increment graph
also handles the projection boundary without claiming an ordinary Jacobian
exists there. Finite-angle attitude injection, Cayley/reset, chord, floor,
clamps and acceptance branches need their own EXACT valid graph constraints;
no generic scalar sector is assumed merely because a map contains a sine.

Use hard finite-horizon IQCs or explicit filter-storage/terminal terms.
An infinite-horizon soft IQC does not automatically justify a three-second
word. Generic nonnegative scalar multipliers are valid for the stated sector;
arbitrary matrix multipliers need a separate sector theorem.

## 4. One master over the literal word

Stack all event variables and source continuation into chi. Let X0 and XN
select endpoints, E_xi chi=0 collect exact linear descriptor relations, and
`chi^T Q_i chi >= 0` be valid hard graph/source IQCs. For fixed rho and a
source edge xi->xi_next form

```
D = XN^T M_next XN - rho X0^T M_current X0
    - W^T Gamma W - C^T Beta C,
L = D + sum_i lambda_i Q_i + Y E_xi + E_xi^T Y^T,
lambda_i >= 0,   L <= 0.
```

On admissible graphs the equality term vanishes and the PLUS IQC sign implies
D<=0. `iqc_master` assembles this sign-correct expression. Strict margin can
be required on the motion/performance subspace, together with uniform storage
coercivity; descriptor redundancy need not itself have artificial strict loss.
A Schur complement or equality-constrained elimination is used when valid,
not a separately rectangularized inverse.

This is affine in the metric/multipliers for fixed coefficient data and rho.
Unknown source coefficients, coupled metric/coefficient decisions, and products
are not automatically a single convex LMI. A vertex-only check needs a proved
convex/affine dependence argument. No such universal reduction is claimed yet.

Search hierarchy: (1) one common joint metric with bounded neutral supply,
(2) low-order source-dependent metric, (3) piecewise metric on a certified
transition graph. Construct/test on analytic source-reachable families, never
fit a metric to captured trajectories. For a state-dependent Riemannian metric,
include its evaluation at the actual next state and a valid path integration
argument. Independent metrics picked anew for each word do not prove switching.

The word includes every due/accepted prediction, floor, S, accelerometer,
magnetometer, reset, projection and tuner event, plus rejected/not-due cases.
Individual events need not contract. Every literal prefix must have a valid
reachability bound that keeps the entire comparison/trajectory inside its
attitude and other domains. Endpoint loss alone is insufficient.

H18->A21 is a real edge: use its actual guard, estimate/covariance reset,
held-bias release and tuner state. Other allowed transitions need edges too.
The source state at one endpoint must be the next word's initial source state.
No hidden metric restart, new bias root or S re-zeroing is permitted.

## 5. Startup capture remains a separate theorem

For the actual Mahony/proxy/learning path, seek a scalar barrier/comparison
function U such that outside a proved landing set each flow/step decreases U
by a positive amount after its bounded disturbance allowance. For a discrete
step decrease at least d>0, the capture count is at most
`ceil((U_initial-U_target)/d)` provided the sublevel set really implies entry
and guards cannot evade progress. The corresponding continuous comparison
needs a genuine time decrease bound and non-Zeno/hybrid accounting.

Certify both measured-period takeover and prior-frequency timeout under the
padded physical family. A timeout proves a scheduler event occurs, not that
its state lies in the certified basin. Covariance consistency is not a hard
error bound. Do not import fresh-Live membership as an assumed success.

## 6. Executable evidence and current blockers

Run from the repository root:

```sh
python -m unittest discover -s tests/ou3_alt_contraction -v
python tools/stability/ou3_alt_contraction/diagnose.py --output /tmp/ou3-alt.json
```

With paired traces from the unchanged observer:

```sh
python tools/stability/ou3_alt_contraction/diagnose.py \
  --maps /tmp/ou3-alt-maps.bin --cov /tmp/ou3-alt-cov.bin \
  --output /tmp/ou3-alt-diagnostic.json
```

The workflow `ou3-alt-contraction` creates those traces and retains source
fingerprints, diagnostic JSON and signed operation ledgers. Captured words
are only falsification/exploratory evidence: no replay fitting or promotion.
The high-precision option re-evaluates the STORED binary32 word product at
80 and 120 digits; it does not reconstruct high-precision nonlinear events.

Current algebra tests check masked/unmasked solve graphs, finite endogenous
increments, shared bias recurrence, retained joint cross terms, IQC sign and
fail-closed status. An exact rational two-state masked-update analogue has

```
A = [[2/3,-1/3],[0,1]], M = [[1,1/4],[1/4,1]], rho=9/10,
A^T M A-rho M-diag(0,1)
  = [[-41/90,-101/360],[-101/360,-43/45]] < 0.
```

Its determinant is 1027/2880>0 and first diagonal negative, while det(M)=15/16.
This proves bounded-supply joint storage resolves the neutral-direction
obstruction in that analogue, NOT in the shipping filter.

The captured baseline has 29 H18 and 359 A21 retained same-mode words; worst
stored-map ratios are approximately 0.999513639919454 and 0.995985717839491.
H18 uses its existing 18-state homogeneous map with held error zero; A21 keeps
all 21 estimated-error states. The observer reconstructs H from P^-1 N and
freezes gain dependence. Neither result is full actual joint24 feasibility,
a source-uniform metric, a chart-retention proof or a hybrid-edge certificate.
The H18 margin of about 0.00048636 is the tighter observed frozen-map margin,
not an available budget for an unbound nonlinear master.

All three individual bias definition validators pass. The shared aggregate
supply builder currently fails fresh-entry hard-coordinate membership; ALT
reports that failure, leaves its source unchanged and does not claim entry.
The existing PE/vector domain consistency issue remains separately open.

### Next falsifiable work, in order

1. Attach a **source-uniform actual** finite-increment joint24 word, retaining
   N/S endogenous increments, true bias, frontend/tuner and guard ancestry.
   Reuse old validated source facts read-only, not old PASS labels.
2. Before rigorous enclosure, run a high-precision complete-word feasibility
   diagnostic of THAT graph. Report worst H/A directions, bias/source ports,
   every operation's margin consumption and hybrid edges. Current captured
   frozen-map ratios do not meet this requirement.
3. Search the common-storage/bounded-supply formulation first. If it fails,
   record the exact limiting direction and compare genuinely different metric
   or source-graph choices before refining. Use the repository two-strike rule.
4. Only after a feasible full master exists: outward-check its source cover,
   compatible endpoints, all prefix reachability, physical working tube,
   projection/bias bounds and finite precision. Then prove startup landing.

No rigorous full-word margin, uniform ultimate bound, certified basin radius,
finite capture time or H18->A21 landing is presently supplied by ALT. The
original proof can continue independently while these obligations are pursued.

## Reference

For differential-to-incremental contraction background: Forni and Sepulchre,
*A differential Lyapunov framework for contraction analysis*, arXiv:1208.2943.
The finite solve-increment and descriptor/IQC algebra above is written out
explicitly; no reference substitutes for its shipping implementation binding.
