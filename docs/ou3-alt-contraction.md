# Parallel OU-III whole-word dissipativity proof

> **Continuation:** start a NEW PR from latest `main` and read
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
IQC, structural anti-shortcut regressions, and actual paired shipping word
experiments are executable. The older frozen-map diagnostic remains separate.
A source-uniform shipping joint24 master has NOT been solved.
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

## 6. Executable actual-word attachment and its limits

The [ALT handover](ou3-alt-contraction-handover.md) gives the commands and current
results. `shipping_probe.py` creates a temporary header overlay that only adds
observations; a byte-recovery test rejects changes to existing shipping lines.
`shipping_word_probe.cpp` executes the real wrapper and Kalman code, including
natural startup and bias release, rather than reconstructing H or multiplying a
frozen homogeneous map. The probe is Linux/host-only and is not linked into
deployment. `--eigen PATH` supports an existing non-default Eigen installation.

A process fork at word entry copies the complete private state, not a hand-
written subset of covariance or frontend fields. The 21 initial estimated-error
perturbations include held bias. Three additional physical-bias-root variants
restart the complete source/estimator history at boot and retain their own
one-time Live S origin. They are not fictitious independent truth coordinates
inserted at a word boundary. There are natural held-mode, active-mode and
H18->A21 words for each of the three separate bias families.

The deterministic source uses a single eight-second harmonic for bounded
translation and roll/pitch motion, and separately defined slow BIAS0/1/2
histories. The generator's all-time triangle bounds are recorded. This probe
subfamily does not cover all physical BRMM, persistent excitation, sensor,
working-domain or startup premises. Conditional Live perturbations are not
hard-entry qualification. No proof-domain restriction is inferred from the
chosen probe amplitudes or diagnostic units.

Every accepted S/accelerometer/magnetometer solve records actual N, S, r, K,
P and the literal pre/post mean vectors. Thus the analyzer can distinguish:

```
q_i = solve(S_i,r_i)                   # ideal real graph on recorded operands
S1*delta_q + delta_S*q0 = delta_r
ideal_delta_correction = N1*delta_q + delta_N*q0
recorded_gain_defect = (K1*r1-K0*r0) - ideal_delta_correction
multiply_add_defect = delta(actual_post_mean-actual_pre_mean) - (K1*r1-K0*r0)
```

Here K*r in the diagnostic is a high-precision multiplication of recorded
binary32 K and r, not the deployed multiplication. The last line is essential:
the observed state increment includes the actual multiply/add rounding.
Innovation asymmetry and these defects are retained, not erased by replacing
the actual update with an information-form identity. Bounds are measured for
the probes only; no deployment-wide roundoff enclosure is claimed.

State tapes also retain prediction, covariance-floor, Joseph, finite reset,
projection, wrapper/tuner and mode events. A signed ledger telescopes the
predeclared diagnostic storage along every aligned literal prefix. Pending
attitude correction is interpreted through the corresponding finite-quaternion
chart extension; actual float injection/reset changes remain in the next event.
Changing a boot root can change scheduler timing. Such event-word mismatches
are reported with their first differing solve occurrence, not silently removed.
Guard completeness and prefix retention over the full family are still open.

The analyzer uses one predeclared joint24 metric, including bias/truth cross
terms, with no replay fitting. Its unsupplied ratios exceed one in all three
modes/families. This is a failed diagnostic metric, not a falsification of the
bounded-supply stability target. Finite coordinate secants are not Jacobian
columns or generalized-eigenvalue maxima. Absolute physical forcing does not
cancel merely because a paired incremental experiment is convenient; truth is
not an estimator trajectory. No supply coefficient or ultimate bound has been
certified by this experiment.

The 80/120-digit checks evaluate the chart and exact finite solve identities on
recorded operands. The nonlinear shipping word itself was executed in binary32.
These checks are not the requested source-uniform, high-precision nonlinear
master. The uninstrumented control option permits a same-compiler comparison
against the original header, separately from the byte-recovery audit. The
executed BIAS2 control matches all 50,775 recorded wrapper states across 75
words bitwise with the same compiler/options:

```sh
python tools/stability/ou3_alt_contraction/shipping_probe.py \
  --output /tmp/alt-control --family 2 --uninstrumented
python tools/stability/ou3_alt_contraction/word_diagnostic.py \
  --directory /tmp/alt/BIAS2 --output /tmp/alt/BIAS2/diagnostic.json \
  --solves --control-directory /tmp/alt-control/BIAS2
```

The independent `ou3-alt-contraction` workflow retains tapes, source fingerprints
and JSON. Twenty-three unit tests check the algebra, attachment, chart extension,
nonzero endogenous terms, roundoff-port distinction, mismatch handling and
fail-closed gates. The prior frozen-map workflow remains available without
changing or depending on any P2/P3/P4/P5 proof stage.

### Current blockers and decision before the next master

The actual probe improves implementation correspondence, but source-uniform
covariance/frontend/guard attachment and consecutive-word compatibility are not
closed. A fresh copy of zero internal increments at each word is not allowed.
Choose and justify an internal-increment storage formulation or a legitimate
same-trajectory physical-error identity before further enclosure work; retain
reference forcing, neutral bias supply and full motion/bias cross information.
The [research ledger](ou3-proof-research-state.md) records the failed metric,
limiting directions, critic and three alternatives.

Individual BIAS0/1/2 definition validators pass, while the aggregate fresh-Live
hard-entry builder still fails. Its gate is not bypassed. Uniform storage,
qualified supply, useful basin/ultimate bound, all nonlinear and hybrid edges,
every-prefix retention, actual startup capture and finite precision remain
open. All ALT theorem flags stay false; existing P4/P5 are untouched.

## Reference

For differential-to-incremental contraction background: Forni and Sepulchre,
*A differential Lyapunov framework for contraction analysis*, arXiv:1208.2943.
The finite solve-increment and descriptor/IQC algebra above is written out
explicitly; no reference substitutes for its shipping implementation binding.
