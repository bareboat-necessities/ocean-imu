# OU-III proof research state

## Current hypothesis

The canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional P3 is
closed and frozen at `delta=1e-18`; P4 is OPEN and P5 is BLOCKED. Production
filter code, the declared source/error domain, numerical quality gates, zero
lever arm, dormant-transparent vibration branch, H18/A21 semantics, and the
separate H18->A21 hybrid are fixed.

The canonical P4 target is the paper's **finite physical true-minus-estimated
state map**, not an estimator-pair shadow and not a state-dependent replacement
theorem. For each admitted complete SEA3 window it must close

`V_{k+N} <= rho V_k + forcing`, `rho<1`,

and, for every prefix,

`V_{k+ell} <= kappa_V V_k + forcing`,

with every prefix inside the certified chart/source domain. Outward
physical-state Jacobians/Clarke Jacobians are permitted as a rigorous
mean-value enclosure of this finite map; they do not redefine the theorem.

## Evidence

The corrected single shipping observer owns one source/tuner/Riccati history
and retains every operation. On the genuine PM+Stokes Hs=1.5 m source the
legal 600-sample linear point words are

- H18: `rho_linear=0.9998658024147671`, 600 predictions, 137 actual-R_S S
  updates, 600 accelerometer updates, 75 vector updates;
- A21: `rho_linear=0.9958536807113242`, 600 predictions, 108 actual-R_S S
  updates, 600 accelerometer updates, 75 vector updates.

The same-observer signed ledgers reproduce these maps within about `2.5e-6`
and telescope to numerical roundoff. H18 gives `rho=0.999863922596`; A21 gives
`rho=0.9958561068`. This validates the point linear observer but is not
universal source coverage.

The exact-Cayley estimator-pair diagnostic is final. On the same A21 word,
full-Phi estimator-pair storage first crosses one at scale 6.5
(`rho_Phi=1.00005078316` for the worse sign), reaches about `1.00215` at scale
8, and reaches `1.13413977623` at the declared limiting scale. Nominal state
reconstruction remains about `7.45e-9`.

This is a real result for the covariance-free **estimator-to-estimator** shadow,
but it is not a counterexample to the physical theorem map. The shadow stores
the difference of two estimates that both receive Kalman corrections. P4 uses
true-minus-estimated physical errors. The maps share the zero-error tangent but
need not share their finite nonlinear map.

## Retained theorem-facing machinery

The branch now contains exact/outward physical-state primitives for:

- Cayley true-minus-estimated prediction with the same committed `tau`, body
  rate and full v/p/S/a_w integrated-OU transition;
- S=0, accelerometer and vector Joseph events with K derived from the SAME
  source P/H/R cell;
- mandatory actual-applied SpectralMSE R_S provenance at every S event;
- deployed quaternion correction and exact physical attitude composition;
- A21 0.5 m/s^2 accelerometer-bias projection with Clarke generalized
  Jacobian when a cell crosses the projection boundary;
- the separate rectangular H18->A21 event;
- literal complete-word differential cocycle in shipping order;
- finite-map generalized mean-value endpoint and every-prefix matrix tests;
- exact whole-word endpoint transport and the joint accelerometer covariance
  channel, in which every later S/R_S event remains inside the suffix map.

P3 also supplies a source-uniform prefix information identity

`Psi_l^T P_l^-1 Psi_l <= P_0^-1`,

and exact event algebra preserves its established full-matrix margin through
prediction, every due S update, every accelerometer update, asynchronous vector
updates, covariance floors and immediate resets. This is a universal algebraic
fact over complete SEA3 and does not require a finite source-word list.

The reset-gauge attachment closes the **zero-error** P3-to-physical-P4 tangent:
P3's margin is invariant under every finite covariance reset congruence, while
for the homogeneous physical event `y(0)=0`, hence `dtheta=0` and `G=I`. This
justifies the physical zero-state tangent `I-KH`. It does **not** close finite
nonlinear reset-coordinate transport or nonzero forcing/noise reset attachment.

## Failure analysis / dead ends

Classification of the finite-shadow expansion: **diagnostic-attachment failure,
not theorem failure**. It invalidates optimizing the estimator-pair raw/Phi
endpoint storage as if it were the physical P4 map. It does not invalidate
frozen P3, the physical event equations, or the finite-state endpoint/prefix
P4 theorem.

A second experiment-design issue is now explicit: a point evaluator that drops
shipping covariance resets by selecting the congruent `G=I` representative can
verify the zero-state tangent, but finite-error energies are not automatically
the paper's quadratic physical-error storage. After a nonzero correction, the
physical post-correction Cayley/additive coordinates and a Joseph-posterior
covariance with its reset congruence omitted are not in the same finite
coordinate gauge unless the exact finite reset-coordinate transport is also
applied. The current reset-gauge lemma deliberately leaves that transport open.

Classification: **finite diagnostic metric-attachment gap**, not theorem
failure and not a P3 failure. It invalidates interpreting reset-normalized
finite point ratios as theorem-grade `V_N/V_0` before the finite coordinate
attachment is closed. It does not invalidate their zero-state Jacobian parity,
the physical state event map itself, or the source/event payload.

**DEAD_ENDS / forbidden rescues**

- duplicated host observers that reconstruct the scheduler independently;
- pure-`e_eta`, scalar Lipschitz, correction-radius, inverse-metric-floor, or
  packet-count-times-worst-remainder routes;
- selected-S replacement words or independent R_S/tuner schedules;
- Schur/elimination of `a_w` or any active state from the final P4 condition;
- 6/9-second optimization of the estimator-pair shadow;
- treating the estimator-pair finite rho as the physical theorem map;
- replacing the 601-sample SEA3 behavior by independent per-sample physical
  boxes;
- inventing an L2/spectral hard amplitude bound from JONSWAP. The theorem
  explicitly states that the spectrum alone does not provide the deterministic
  pathwise bound for `X^s_SEA3`;
- treating zero-error reset-gauge isometry as if it had already proved finite
  nonlinear reset-coordinate attachment.

## Current limiter

There is still **no theorem-grade low-pessimism point energy result for the
actual finite P4 storage**. The new same-observer physical payload and physical
state map are useful, but finite energy interpretation requires the exact
finite reset-coordinate attachment first. Going directly from a reset-normalized
finite ratio to a source-uniform interval proof would violate the research
protocol.

For universal closure, the complete SEA3 hard-window behavior `B^601_SEA3` is
compact and requires one common phase/parameter/response witness, but its
validated correlated outer-enclosure oracle is not implemented. Normal-Live
sample caps are explicitly not sufficient membership conditions. Conditional
P4 does not require the separate physical SEA0->SEA3 left inclusion, but it may
not replace the admitted SEA3 history by independent boxes.

A promising way to avoid materializing the 10,818-dimensional sampled behavior
is to use P3's universal event algebra and the complete-word endpoint
reduction. The joint accelerometer channel already gives

`d_acc^T P_N^-1 d_acc <= q^T R_acc^-1 q`

for the entire stacked nonlinear accelerometer history, with no packet-count
multiplier and with every later actual-R_S S event inside the suffix. The open
question is whether the stacked nonlinear graph and endpoint cross/boundary
terms can be dominated by the same complete-word information on a usable
finite cell.

## Independent critic pass

The strongest immediate reason to distrust the present point experiment is the
finite reset-coordinate mismatch, not a numerical rho. Zero-error parity alone
cannot justify evaluating finite post-correction physical errors in an
unreset Joseph covariance metric.

Three legitimate routes from this limiter are:

1. derive the exact finite reset-normalized physical coordinate and prove that
   the paper's chosen quadratic storage is represented without an error-dependent
   metric; then rerun the same-observer physical finite-map point test;
2. keep the paper's physical Cayley/additive coordinates and transport the
   source-indexed P3 metric through an exact source-only congruent representative
   at each boundary, proving that the representative does not depend on the
   finite error state;
3. if neither is compatible with the paper's stated source-indexed quadratic
   storage, use the whole-word endpoint identity directly and prove finite
   dissipation without introducing a state-dependent Riccati metric.

None permits state elimination, source replacement, independent R_S schedules,
or changing the theorem domain.

## Next falsifiable experiment

The same-observer payload path is implemented and should still be run because
it verifies event extraction, actual-R_S retention, zero-state H/K/Jacobian
parity, and the exact physical state map. Its finite ratios are **diagnostic
only** until reset-coordinate attachment closes.

Before using those finite ratios to select or reject a P4 cell, derive and test
the exact finite reset-coordinate identity at a single Joseph event:

- start from exact Cayley true-minus-estimated pre-correction error and the same
  P/H/R cell;
- apply the exact physical correction `E_plus=E Q(Ky)^-1` and additive
  true-minus-estimated corrections;
- relate the resulting finite physical coordinate to the covariance reset
  `P_R=G P_J G^T` without replacing `G` by `I` unless the corresponding state
  coordinate is transformed by the exact same gauge map;
- prove the energy identity in both H18 and A21, including the A21 bias
  projection branch;
- verify its derivative at zero reduces to the already-closed reset-gauge
  tangent.

Only after that single-event finite metric attachment passes should the full
same-observer physical word be assigned theorem-storage endpoint/prefix ratios.
If the corrected physical point result is strict, proceed to the source-uniform
generalized mean-value endpoint/prefix enclosure; if it is not strict, stop and
replan before adding an enclosure lemma.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
