# OU-III proof research state

## Current hypothesis

The canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional P3 is
closed and frozen at `delta=1e-18`; P4 is OPEN and P5 is BLOCKED. The user has
authorized one production/proof-domain tightening: the accelerometer-bias
projection radius is now `0.4 m/s^2` instead of `0.5 m/s^2`; the Normal-Live
interior active-state bound is `0.35 m/s^2` so the previous `0.05 m/s^2`
projection margin is preserved, and the startup/handoff accelerometer-bias
error envelope is `0.4 m/s^2`. No other production filter tuning, source-domain
parameter, or numerical quality gate is changed. Zero lever arm,
dormant-transparent vibration branch, H18/A21 semantics, and the separate
H18->A21 hybrid remain fixed.

The canonical P4 target is the paper's finite physical true-minus-estimated
state map. Point/shadow/replay calculations are falsification and architecture
selection only; they cannot promote P4 or replace the complete SEA3 source.
Every valid accelerometer correction, every due S=0 correction with the actual
applied anisotropic SpectralMSE `R_S`, every applicable vector event, full Q,
covariance floors, and immediate resets remain in the literal word.

For one admitted proof window P4 must establish a strict finite endpoint
inequality, a finite every-prefix gain, and every-prefix chart/domain retention.
The proof window length is a theorem quantity, not intrinsically fixed to the
600-sample / 3 s point diagnostic. The paper's S-observability construction
explicitly states that a longer proof window may credit wider source-realizable
S-event separation. Any longer-window experiment must nevertheless be the
actual physical same-history complete word; the retired 6/9-second
estimator-pair-shadow optimization remains forbidden.

## Frozen linear evidence

The corrected single shipping observer owns one source/tuner/Riccati history
and retains every operation. On the genuine PM+Stokes Hs=1.5 m source the
pre-0.4 legal 600-sample point words were

- H18: `rho_linear=0.9998658024147671`, 600 predictions, 137 actual-R_S S
  updates, 600 accelerometer updates, 75 vector updates;
- A21: `rho_linear=0.9958536807113242`, 600 predictions, 108 actual-R_S S
  updates, 600 accelerometer updates, 75 vector updates.

The same-observer event ledgers reproduce these maps within a few `1e-6` and
telescope to roundoff. This establishes that the A21 obstruction is not a
missing linear-information or observer-parity problem. It is not universal
source coverage.

On the user-authorized 0.4 projection/domain head, the non-promoting physical
finite-map harness passes strict zero-state parity. Its reset-normalized point
linear ratios are approximately `0.9998358723128917` (H18) and
`0.9958491837932946` (A21). P3 remains frozen at `delta=1e-18`; this production
change does not authorize retuning or replacing P3.

P3 also supplies the source-uniform prefix information inequality

`Psi_l^T P_l^-1 Psi_l <= P_0^-1`

and exact event algebra preserves the established full-matrix margin through
prediction, every due S update, every accelerometer update, asynchronous vector
updates, covariance floors and finite reset congruences.

## Nonlinear point evidence

### Estimator-pair shadow

The exact-Cayley covariance-free estimator-pair shadow remains a historical
non-promoting diagnostic. H18 stays contractive. In A21, raw and full-Phi
storage first cross one at reliable scale 6.5; full-Phi reaches about
`1.13413977623` at the old declared limiting scale. This is real for the
estimator-to-estimator shadow but is not the paper's physical error map.

### Physical reset-normalized diagnostic

A single-observer physical source/event payload drives the exact
true-minus-estimated Cayley state equations. Host-reconstructed measurement H
is projected back only to exact shipping structural support (S selector,
skew vector block, SO(3) accelerometer a_w block, A21 identity b_a block), with
P/Q/schedule/timestamps and every actual-R_S S event retained bit-for-bit.
The original strict zero-state event and whole-word parity tolerances are not
relaxed.

With the 0.4 projection/domain change, H18 remains endpoint-contracting over all
retained tested point scales. A21 still first crosses one at scale `8.0`; the
bias projection is inactive in the problematic cases. The worst retained-domain
A21 endpoint ratio in this diagnostic is about `1.0860152320`. Therefore the
0.4 clamp tightens the declared A21 domain as requested but does not manufacture
or explain away the finite A21 obstruction.

Finite ratios from this reset-normalized point diagnostic remain non-promoting
until the finite source-indexed reset/metric attachment is closed. They are
nevertheless useful mechanism/falsification evidence because the physical state
equations and same-source Joseph data are retained.

The physical full-Phi endpoint was also tested. It improves A21 slightly but
does not solve it. A simple same-endpoint constant discrete-converse metric
obtained from a single A21 point map was also tested and remains expansive for
one sign at relevant finite scales. Neither is a proof route.

## Exact signed-information ledger

The branch contains the exact complete-word Joseph/reset identities. For one
same-cell Joseph event,

`y = H e + eta`, `S = H P H^T + R`, `K = P H^T S^-1`, `t=e-Ky`,

and

`V_J - V = - y^T S^-1 y + eta^T R^-1 eta`.

For the shipping reset `P_R=G P_J G^T`, if the exact finite physical error is
`e_R=G t+rho`, `b=G^-1 rho`, then

`V_R - V_J = 2 t^T J_J b + b^T J_J b`.

Thus a correction/reset contributes exactly

`Delta V = -I_y + E_eta + X_reset + E_reset`.

Every S=0 event has `eta=0` exactly and therefore contributes favorable
information with its actual applied R_S. Prediction with PSD Q and the PSD a_w
covariance floor are non-increasing in the corresponding moving linear
information energy. No packet-count remainder bound is needed or allowed.
Exact rational H18/A21 tests close these identities.

## A21 mechanism diagnosis

The user-authorized change from a `0.5` to `0.4 m/s^2` accelerometer-bias
projection does **not** remove the observed A21 finite point difficulty. At the
old scale-6.5 estimator-pair crossing the b_a component of the point worst
direction was only about `0.086 m/s^2`; in the 0.4 physical diagnostic the
projection is likewise inactive at the problematic scales. The change is a
legitimate production/proof-domain tightening requested independently of the
proof result, not a substitute for solving A21.

The signed point ledger localizes the A21 balance at scale +6.5 approximately
as follows (diagnostic values, not a certificate):

- initial information energy: about 42.25;
- prediction/process covariance contribution: about `-0.657`;
- covariance-floor contribution: about `-0.0011`;
- 600 accelerometer corrections together: about `+0.662`;
- vector and S corrections together: about `-0.0050`;
- total finite reset contribution: small and mildly favorable;
- nonlinear prediction excess over the linear prediction: only about `+1.8e-4`.

Inside the accelerometer contribution, the dominant adverse term is the
finite-angle nonlinear residual energy `eta^T R_acc^-1 eta`; the reset cross and
reset-defect energies are tiny by comparison. Therefore the controlling
mechanism is finite accelerometer curvature in an A21 direction where attitude,
latent acceleration, and accelerometer bias nearly cancel the first-order
accelerometer residual. Reset-radius tightening and inverse-metric-floor bounds
target the wrong mechanism.

The exact accelerometer model explains this cancellation. With lever arm off,

`y = (E-I) f_hat + E R_hat delta_a_w + delta_b_a`.

Its first-order row is

`H e = [c]_x f_hat + R_hat delta_a_w + delta_b_a`.

A single short word can therefore contain a direction in which the three
first-order pieces nearly cancel while second-order attitude curvature remains.
The long-lived residual bias (`tau_b` about 5000 s), the much faster latent OU
acceleration (order 2 s on the observed word), and changing attitude/specific-
force geometry provide the natural mechanism that can break this cancellation
over a longer source-correlated window.

## Theorem-facing A21 routes

The published active-bias theorem permits two full-state A21 routes:

1. full finite-window PE of attitude/gyro/accelerometer-bias coordinates; or
2. reduced attitude/gyro PE together with finite bounded `tau_b`, using the
   intrinsic residual-bias contraction/detectability fallback.

Route 2 retains all 21 states. It is not state elimination and does not turn
b_a into a removed coordinate. The current P4 machinery must not accidentally
over-constrain the paper by requiring the shipping P^-1 quadratic to contract
on every 3-second A21 word if a longer complete word or the finite-tau_b
full-state cascade is the theorem's valid construction.

## Retained theorem-facing machinery

The branch currently retains:

- exact/outward Cayley physical prediction with committed h/tau/body rate and
  the full v/p/S/a_w OU chain;
- source-cell S/accelerometer/vector Joseph events with K derived from the same
  P/H/R cell;
- actual applied R_S provenance on every S event;
- exact deployed quaternion correction and A21 `0.4 m/s^2` bias projection with
  Clarke generalized Jacobian;
- the separate H18->A21 rectangular hybrid event;
- literal complete-word differential cocycle and generalized mean-value bridge;
- exact whole-word endpoint transport with later actual-R_S suffix maps;
- joint accelerometer covariance channel;
- exact signed Joseph/reset information ledger;
- exact homogeneous Cayley residual-sector factorization, under which the A21
  delta-b_a term is exactly linear and contributes zero nonlinear eta.

These are structural components only. Source-uniform finite P4 endpoint,
every-prefix gain, every-prefix domain retention, and finite reset/source-metric
attachment remain open.

## Dead ends / forbidden rescues

- duplicated host observers or independent schedulers/Riccati histories;
- estimator-pair shadow promoted as the theorem map;
- raw or full-Phi estimator-pair endpoint optimization;
- physical full-Phi treated as an A21 solution after its point failure;
- arbitrary single-map converse-metric fitting;
- 6/9-second optimization of the estimator-pair shadow;
- selected-S words, independent R_S/tuner schedules, or independent per-sample
  source boxes;
- state elimination / a_w Schur final certificate;
- scalar Lipschitz, pure-e_eta, correction-radius, inverse-metric-floor, or
  packet-count-times-worst-remainder bounds;
- inventing an L2/spectral pathwise bound from JONSWAP;
- further proof-driven filter/domain tightening beyond the user-authorized 0.4
  projection change merely to make A21 easier.

## Current controlling blocker

The upstream complete-source contract still explicitly reports that the
correlated finite-window SEA3 realization is not materialized as an executable
outward source family. Parameter compactness, RAO/moment envelopes, hard
pathwise acceleration/body-rate caps, frontend parity, and adaptive-state
rate/jump bounds do not replace a same-history transition for the correlated
source state. P4 may not substitute independent per-sample acceleration/rate
boxes, independent tuner/R_S schedules, a finite harmonic/grid surrogate,
replay, or another source language.

There is also still no theorem-grade finite reset/source-metric attachment for
the actual finite P4 storage. Zero-error reset-gauge parity is valid, but finite
post-correction physical coordinates and the source-indexed quadratic metric
must be attached exactly before a finite replay ratio can be interpreted as the
paper's theorem storage ratio.

A promising way to avoid materializing a giant sampled behavior vector is to
use P3's universal event algebra, the exact residual sector, and the
complete-word endpoint/signed-information reductions. The joint accelerometer
channel already gives

`d_acc^T P_N^-1 d_acc <= q^T R_acc^-1 q`

for the entire stacked nonlinear accelerometer history, with no packet-count
multiplier and every later actual-R_S S event inside the suffix. The open
question is whether the stacked homogeneous nonlinear sector and endpoint
cross/boundary terms can be dominated by the same complete-word information on
a usable finite cell.

## Independent critic pass

The strongest reason to abandon a longer-window route would be a same-history
**physical** A21 point result that remains expansive on source-contiguous 6 s
and 9 s windows, especially after recomputing the linear maximizing direction
for each longer window. In that case extra source PE is not repairing the
finite-curvature mechanism and adding universal enclosure machinery would be
wasted effort.

Conversely, a strict longer physical point result is not proof. It only shows
that the paper-permitted longer complete-window route has enough margin to
justify the next universal enclosure. The universal proof must consume the
complete SEA3 source language, not the replay used by this falsification test.

If longer physical windows fail, the next qualitatively different route is the
paper's finite-tau_b full-state A21 cascade/detectability construction. That
route must produce a full-rank 21-state Lyapunov inequality and retain all
couplings; it may not eliminate b_a or weaken A21 into an 18-state result.

## Next falsifiable experiment

Finish the existing single-observer source-contiguous 1200-sample (6 s) and
1800-sample (9 s) **physical** finite-map diagnostic under the 0.4 projection
and proof-domain change. For each horizon:

1. recompute the worst linear H18/A21 direction from the exact shipping map and
   moving covariance boundaries for that same horizon;
2. require every valid accelerometer update and every due S update with actual
   applied R_S; reject hybrid/mode-changing windows exactly as for 3 s;
3. emit the selected source/event payload from the same observer;
4. canonicalize only the host-reconstructed measurement geometry, retaining
   strict zero-state event/whole-word parity;
5. execute the exact physical true-minus-estimated finite map for both signs
   over the declared finite cell;
6. report endpoint ratios, worst prefix gain, domain retention, event counts,
   projection activity, and the signed-information balance;
7. do not assert contraction as a CI infrastructure condition and do not
   promote P4 from this replay.

This is a theorem-facing window-length feasibility test, not the retired
estimator-shadow multiword optimization. If a longer physical window is strict,
the next proof implementation is the same-history complete-SEA3
signed-information/finite-map enclosure over that chosen finite window. If
neither 6 s nor 9 s is strict, stop this route and build the finite-tau_b
full-state A21 cascade certificate instead.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
