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
legal 600-sample linear point words are approximately

- H18: `rho_linear=0.9998654`, 600 predictions, 137 actual-R_S S updates,
  600 accelerometer updates, 75 vector updates;
- A21: `rho_linear=0.9958522`, 600 predictions, 108 actual-R_S S updates,
  600 accelerometer updates, 75 vector updates.

The same-observer signed ledgers reproduce these maps within a few `1e-6` and
telescope to numerical roundoff. This validates the point linear observer but
is not universal source coverage.

The exact-Cayley estimator-pair diagnostic is now final. On the same A21 word,
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

## Failure analysis / dead ends

Classification of the finite-shadow expansion: **diagnostic-attachment failure,
not theorem failure**. It invalidates optimizing the estimator-pair raw/Phi
endpoint storage as if it were the physical P4 map. It does not invalidate
frozen P3, the physical event equations, or the finite-state endpoint/prefix
P4 theorem.

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
  pathwise bound for `X^s_SEA3`.

## Current limiter

There is still **no low-pessimism point feasibility result for the actual
physical finite P4 map**. Going directly from the estimator-pair result to a
source-uniform interval proof would violate the research protocol.

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

Strongest reason to abandon the present route would be a physical finite-map
point result with `rho>=1` well inside the candidate cell after exact event
parity is verified. That would show that interval sharpening cannot rescue the
controlling theorem inequality.

Qualitatively different outcomes after the physical-map probe are:

1. Physical endpoint and prefix ratios are strictly feasible: proceed to one
   same-history complete-word mean-value enclosure, exploiting P3 prefix
   information and the joint accelerometer channel rather than source boxes.
2. Endpoint is feasible but a prefix leaves the domain or has excessive gain:
   solve the actual prefix/domain-retention obligation; do not alter the
   endpoint metric or source.
3. Physical map is not feasible in the source-indexed quadratic metric but its
   full-state differential cocycle is uniformly stable: a different full-rank
   source-indexed/converse Lyapunov construction may be justified.
4. Physical finite map itself has an admissible expanding direction with no
   compatible full-state metric/path argument: report theorem-failure evidence
   rather than weaken the domain or eliminate states.

## Next falsifiable experiment

Extend the **existing single shipping operation observer** with an optional
selected-word source-payload trace. Do not create another estimator, scheduler,
or Riccati recursion. For every event record only the source quantities needed
by the already-tested theorem event functions: prediction `omega_hat`, `h`,
committed `tau`/A21 `tau_b`; pre-Joseph full P; actual R_S; accelerometer
`f_hat,R_hat`; accepted vector geometry; and the A21 true-bias/projection
coordinate required by the event contract.

Consume that trace with the existing Python physical prediction/Joseph event
functions and compose the finite true-minus-estimated map in exact Cayley
coordinates. On the same H18/A21 limiting words, test both signs and the
candidate-domain scales while reporting:

- finite endpoint `V_N/V_0` in the paper's source-indexed quadratic metric;
- maximum every-prefix `V_l/V_0` and chart/domain retention;
- zero-state Jacobian parity against the verified linear complete-word map;
- exact event counts and actual-R_S history;
- A21 projection branch/Clarke activity.

This is a non-promoting feasibility/falsification test only. If it is strict,
the next implementation is the source-uniform generalized mean-value endpoint
and prefix enclosure using the same complete-word event algebra. If it is not
strict, stop and replan before adding any enclosure lemma.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
