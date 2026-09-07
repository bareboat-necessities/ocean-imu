# OU-III proof research state

## Current status

The canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional P3
is closed and frozen at `delta=1e-18`; P4 is OPEN and P5 is BLOCKED. Production
filter code, the declared source/error domain, numerical quality gates, zero
lever arm, dormant-transparent vibration branch, H18/A21 mode semantics, and
the H18->A21 hybrid separation are fixed.

The theorem objective is full-state nonlinear finite-window UES/ISS. A P4
metric must be full rank on every active coordinate and prove the same physical
true-minus-estimated H18/A21 error dynamics. State elimination, source
replacement, replay witnesses, independently bounded tuner/R_S sequences,
selected-S words, packet-count remainder sums, correction-radius fallbacks, and
domain shrinkage are not admissible.

## Frozen linear reference

The corrected single shipping observer owns one source/tuner/Riccati history and
retains every executed operation, including every due S=0 update with the actual
applied R_S matrix generated from the smoothed scalar `RS_applied` channel,
realized pseudo cadence, and deployed standard-deviation factors
`[0.72,0.72,1]`.

On the genuine PM+Stokes Hs=1.5 m same-history source, the legal 600-sample
point words are:

- H18: `rho_linear=0.9998658024147671`, with 600 predictions,
  600 accelerometer updates, 137 S updates, and 75 vector updates.
- A21: `rho_linear=0.9958536807113242`, with 600 predictions,
  600 accelerometer updates, 108 S updates, and 75 vector updates.

The selected-direction event ledgers reproduce those maps within about `2.5e-6`
and telescope to numerical roundoff. These are non-promoting feasibility checks,
not universal source coverage.

## Point-shadow diagnostics: what they do and do not establish

The covariance-free frozen-gain shadow uses one shipping estimator for the
source, covariance, gains, branch decisions, and actual-R_S word. Its nominal
reconstruction error is about `2.98e-8` in H18 and `7.45e-9` in A21. It was
useful for finding real harness defects: staged tuner/R_S commit ordering,
covariance-changing perturbation injection, and the finite attitude-chart
mismatch.

The final parity check injects the requested **exact Cayley** attitude coordinate
`c=2 tan(theta/2) u` by converting it to the equivalent rotation vector accepted
by the host shadow. Focused workflow run `34071140640` retained the verified A21
word `(600 prediction, 108 S, 600 accelerometer, 75 vector events)` and found,
for the estimator-pair full-Phi storage:

- scale 6: `rho_Phi=0.999459385872` / `0.999386787415`;
- scale 6.5: `1.00005078316` / `0.999995589256`;
- scale 7: `1.00067698956` / `1.00064611435`;
- scale 8: `1.00212538242` / `1.00215435028`;
- declared limiting scale `37.714067586719345`: worst
  `rho_Phi=1.13413977623`.

Those numbers are real for the **incremental estimator-to-estimator shadow**.
They are not a theorem-grade counterexample to the physical P4 map. The shadow
stores `shadow estimate - nominal estimate` and both estimates receive Kalman
corrections, whereas the theorem uses physical true-minus-estimated errors with
`R_true=E R_hat`, `delta a_w=a_true-a_hat`, and the true trajectory does not
receive estimator corrections. The two maps share the zero-error tangent, which
explains the linear agreement, but their finite nonlinear maps are not
identical.

Therefore the shadow result may falsify a proposed *diagnostic attachment* or
warn against optimizing a finite endpoint storage, but it may not be used to
claim that physical raw/Phi P4 is false. No P4 architecture may be promoted or
rejected solely from this estimator-pair rho.

## Exact physical-error convention retained for P4

The theorem-facing accelerometer coordinate is

`R_true = E R_hat`, `delta a_w = a_true-a_hat`,

with exact homogeneous residual

`y=(E-I) f_hat + E R_hat delta a_w + delta b_a`.

The shipping tangent is

`H z=[c]x f_hat + R_hat delta a_w + delta b_a`.

The exact measurement-linearizing shift remains

`Q_aw=R_hat^T E R_hat`,

`e_eta=R_hat^T((E-I)-[c]x)f_hat`,

`epsilon_aw=(Q_aw-I)delta a_w+e_eta`,

`Phi(z)=z+E_aw epsilon_aw`,

so `y=H Phi(z)` with the ORIGINAL shipping H/P/R/K/S. For a correction
`d=K y`, the true-minus-estimated additive error subtracts `d` and the exact
attitude error satisfies `E_plus=E Q(d)^-1`. Prediction retains the literal
full `F E_aw` transport, including v/p/S/a_w rows. H18->A21 remains a separate
21x18 physical/covariance hybrid.

## Dead ends and forbidden rescues

- duplicated host observers that reconstruct the scheduler independently;
- pure-`e_eta`, scalar Lipschitz, correction-radius, metric-floor, or
  packet-count-times-worst-remainder routes;
- selected-S replacement words or independent R_S schedules;
- reduced/end-point Schur elimination of `a_w` or any other active state from
  the final P4 coercivity condition;
- 6/9-second point-window optimization of the estimator-pair shadow;
- treating the estimator-pair finite rho as the physical theorem map.

The old `ou3-p4-a21-failure-mechanism` Schur/subdirection diagnostic is
historical debugging only and must not feed canonical P4.

## Canonical full-state differential architecture

The current canonical P4 architecture is the full-rank pullback differential
metric

`M(z,zeta)=D Phi(z,zeta)^T P(zeta)^-1 D Phi(z,zeta)`.

`D Phi` is block triangular and its a_w diagonal block is the orthogonal
`Q_aw`, so `det D Phi=1` on the finite Cayley chart. No active state is
eliminated. At zero error the metric reduces exactly to the frozen P3 moving
metric.

For one complete physical word `F_W` the controlling inequality is

`rho M_0 - D F_W(z)^T M_1 D F_W(z) > 0`, `0<rho<1`,

uniformly over every admitted complete SEA3 same-history word and every state
in one certified finite-error cell, separately in H18 and A21. The H->A
rectangular differential event is separate.

The existing differential-word cocycle composes literal event Jacobians in
shipping order and requires the same source token for every event. Every due
S event must carry actual-applied SpectralMSE R_S provenance. The terminal gate
is full interval LDLT, not a scalar packet norm.

## Current full-state differential obligations

The theorem-facing differential Joseph event module now:

- uses exact Cayley physical true-minus-estimated state;
- derives H/S/K from the SAME source-correlated P/H/R operation cell;
- forbids an independently supplied theorem K;
- requires the actual-applied-R_S provenance token on S events;
- differentiates the exact deployed quaternion correction outward;
- keeps all H18/A21 coordinates.

The following obligations remain OPEN and are controlling:

1. **Prediction differential.** Supply the exact physical true-minus-estimated
   prediction Jacobian for the same source cell, including attitude/gyro-bias
   propagation, the full v/p/S/a_w OU map, active b_a OU propagation, and the
   source-committed schedule. The covariance-floor event has identity state
   Jacobian but remains in the covariance/source word.
2. **A21 b_a projection hybrid.** Shipping projects the estimated b_a state to
   the 0.5 m/s^2 ball. Normal-Live declares nominal `||b_a_hat||<=0.45`, giving
   only 0.05 m/s^2 nominal interior margin; the declared nonlinear error set can
   cross the projection boundary. Projection inactivity is therefore NOT
   source-uniformly proved. P4 must either derive a stronger invariant from
   existing assumptions or include the exact nonsmooth projection hybrid /
   generalized Jacobian. The domain must not be shrunk to avoid it.
3. **Source-correlated finite event cells.** P/H/R, geometry, schedule and state
   cells must be generated from one complete SEA3 history. Independent boxes
   may not be multiplied into a fake word.
4. **H18->A21 differential hybrid.** Supply the actual 21x18 lift and covariance
   seed/reset event; do not reuse an H18 ceiling as A21.
5. **Endpoint pullback cells.** Enclose P0^-1/P1^-1 and D Phi at the same source
   endpoints and compose one full interval Jacobian through every prediction,
   actual-R_S S event, accelerometer event, vector event, floor and reset.
6. **Terminal full-matrix test.** Run interval LDLT for the candidate cells
   `[30,25,20,15]` degrees, selecting only the widest cell for which the full
   same-history H18 and A21 inequalities close. A green structural test alone
   cannot promote P4.

## Independent critic pass

The strongest reason this architecture could still fail is not the point-shadow
rho. It is source-uniform interval dependency combined with the nonsmooth A21
bias projection and source-correlated Kalman gains. A proof that boxes P/H/R/K,
geometry, or schedule independently would almost certainly manufacture excess
width and would violate complete SEA3 even if it numerically passed.

Three legitimate alternatives if the pullback metric fails after a faithful
same-history enclosure are:

1. a different **full-rank** source/state-dependent Riemannian metric;
2. a full-state converse/path-memory Lyapunov metric built from the same
   differential cocycle;
3. theorem-failure analysis if no uniformly coercive full-state metric can
   close on the declared domain.

None permits state elimination, source replacement, or domain shrinkage.

## Next falsifiable experiment

Complete the theorem-facing event library before constructing a universal word:

- verify every nonlinear Joseph event has the correct zero-cell tangent against
  the literal shipping H/K/reset map;
- implement and test the exact physical prediction Jacobian;
- implement the A21 projection hybrid/generalized Jacobian (or prove its
  inactivity from existing declarations without changing those declarations);
- implement the H18->A21 rectangular differential event.

Only after those event maps pass exact point/interval algebra tests should one
bind them to one source-correlated complete SEA3 word and evaluate the master
full-matrix inequality. Do not return to finite endpoint-storage optimization.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
