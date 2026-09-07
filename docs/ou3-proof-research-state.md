# OU-III proof research state

## Current hypothesis

The canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional P3 is
closed and frozen at `delta=1e-18`; P4 is OPEN and P5 is BLOCKED. Production
filter code, the declared source/error domain, numerical quality gates, zero
lever arm, dormant-transparent vibration branch, H18/A21 semantics, the
accelerometer-bias clamp, and the separate H18->A21 hybrid are fixed.

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
legal 600-sample point words are

- H18: `rho_linear=0.9998658024147671`, 600 predictions, 137 actual-R_S S
  updates, 600 accelerometer updates, 75 vector updates;
- A21: `rho_linear=0.9958536807113242`, 600 predictions, 108 actual-R_S S
  updates, 600 accelerometer updates, 75 vector updates.

The same-observer event ledgers reproduce these maps within a few `1e-6` and
telescope to roundoff. This establishes that the A21 obstruction is not a
missing linear-information or observer-parity problem. It is not universal
source coverage.

P3 also supplies the source-uniform prefix information inequality

`Psi_l^T P_l^-1 Psi_l <= P_0^-1`

and exact event algebra preserves the established full-matrix margin through
prediction, every due S update, every accelerometer update, asynchronous vector
updates, covariance floors and finite reset congruences. P3 is not being
retuned or reopened.

## Nonlinear point evidence

### Estimator-pair shadow

The exact-Cayley covariance-free estimator-pair shadow remains a historical
non-promoting diagnostic. H18 stays contractive. In A21, raw and full-Phi
storage first cross one at reliable scale 6.5; full-Phi reaches about
`1.13413977623` at the declared limiting scale. This is real for the
estimator-to-estimator shadow but is not the paper's physical error map.

### Physical reset-normalized diagnostic

A single-observer physical source/event payload now drives the exact
true-minus-estimated Cayley state equations. Host-reconstructed measurement H
is projected back only to exact shipping structural support (S selector,
skew vector block, SO(3) accelerometer a_w block, A21 identity b_a block), with
P/Q/schedule/timestamps and every actual-R_S S event retained bit-for-bit.
The original strict zero-state event and whole-word parity tolerances are not
relaxed.

Finite ratios from the reset-normalized point diagnostic remain non-promoting
until the finite source-indexed reset/metric attachment is closed. They are
nevertheless useful mechanism/falsification evidence because the physical
state equations and same-source Joseph data are retained.

The physical full-Phi endpoint was also tested. It improves A21 slightly but
does not solve it: around scale 6.5 the worse sign is still essentially at or
just above one, and the ratio rises above one at larger scales. Therefore
full-Phi is not the A21 solution and will not be promoted as one.

A simple same-endpoint constant discrete-converse metric obtained from a
single A21 point map was also tested. It remains expansive for one sign at the
relevant finite scales. Arbitrary pointwise metric fitting is therefore not a
proof route.

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

The first reliable A21 finite loss is **not** caused by the 0.5 m/s^2
accelerometer-bias projection. At scale 6.5 the b_a component of the point
worst direction is only about 0.086 m/s^2; the projection is inactive. The user
explicitly rejected lowering the clamp as a proof fix, and no clamp/domain
change will be made.

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
reset-defect energies are tiny by comparison. At scale 8 the same imbalance is
larger. Therefore the controlling mechanism is finite accelerometer curvature
in an A21 direction where attitude, latent acceleration, and accelerometer bias
nearly cancel the first-order accelerometer residual. Reset-radius tightening,
metric-floor bounds, and clamp changes target the wrong mechanism.

The exact accelerometer model explains this cancellation. With lever arm off,

`y = (E-I) f_hat + E R_hat delta_a_w + delta_b_a`.

Its first-order row is

`H e = [c]_x f_hat + R_hat delta_a_w + delta_b_a`.

A single short word can therefore contain a direction in which the three
first-order pieces nearly cancel, while the second-order attitude curvature
remains. The long-lived residual bias (`tau_b` about 5000 s), the much faster
latent OU acceleration (order 2 s on the observed word), and changing attitude /
specific-force geometry provide the natural mechanism that can break this
cancellation over a longer source-correlated window.

## Theorem-facing A21 routes

The published active-bias theorem already permits two full-state A21 routes:

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
- exact deployed quaternion correction and A21 0.5 m/s^2 bias projection with
  Clarke generalized Jacobian;
- the separate H18->A21 rectangular hybrid event;
- literal complete-word differential cocycle and generalized mean-value bridge;
- exact whole-word endpoint transport with later actual-R_S suffix maps;
- joint accelerometer covariance channel;
- exact signed Joseph/reset information ledger described above.

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
- changing the accelerometer-bias clamp or shrinking the declared domain to
  make A21 easier.

## Independent critic pass

The strongest reason to abandon a longer-window route would be a same-history
**physical** A21 point result that remains expansive on source-contiguous 6 s
and 9 s windows, especially after the linear maximizing direction is recomputed
for each longer window. In that case extra source PE is not repairing the
finite curvature mechanism and adding universal enclosure machinery would be
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

Use the existing **single shipping operation observer** and genuine coupled
source, but scan source-contiguous 1200-sample (6 s) and 1800-sample (9 s)
windows in addition to the existing 600-sample diagnostic. For each horizon:

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

This is a theorem-facing **window-length feasibility** test, not the retired
estimator-shadow multiword optimization. The paper explicitly allows longer
proof windows in its finite-window observability construction.

If a longer physical window is strict, the next proof implementation is the
same-history complete-SEA3 signed-information / finite-map enclosure over that
chosen finite window. If neither 6 s nor 9 s is strict, stop this route and
build the finite-tau_b full-state A21 cascade certificate instead.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
