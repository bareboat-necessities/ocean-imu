# OU-III proof research state

## Current hypothesis

The canonical source remains `COMPLETE_SEA3_NORMAL_LIVE_WORD`. Conditional P3 is
closed and frozen at `delta=1e-18`; P4 is open and P5 is blocked. Production
filter code, the declared source/error domain, quality gates, H18/A21 mode
semantics, zero lever arm, and dormant-transparent vibration branch are fixed.

The theorem objective is full-state nonlinear finite-window UES/ISS. Existing
moving-Riccati storage constructions are proof hypotheses, not the theorem.
A replacement P4 metric is admissible only if it is full rank on all active
coordinates and proves the same full H18/A21 error dynamics. Eliminating a
state, replacing the source, selecting favorable S events, replay fitting, or
shrinking the declared domain is not admissible.

## Evidence

The corrected single shipping observer is the point-diagnostic reference. It
uses one source/tuner/Riccati history and retains every shipping operation,
including every due S=0 update with the actual applied R_S matrix generated from
the smoothed scalar `RS_applied` channel, realized pseudo cadence, and deployed
standard-deviation factors `[0.72,0.72,1]`.

For the genuine PM+Stokes Hs=1.5 m same-history source, legal 600-sample words
give:

- H18: complete-word linear `rho=0.9998658024147671`; 600 predictions,
  600 accelerometer updates, 137 S updates, 75 vector updates.
- A21: complete-word linear `rho=0.9958536807113242`; 600 predictions,
  600 accelerometer updates, 108 S updates, 75 vector updates.
- The selected-direction event ledgers reproduce those maps within about
  `2.5e-6` and telescope to numerical roundoff.
- The covariance-free frozen-gain nonlinear shadow reconstructs the nominal
  shipping state to `2.98e-8` in H18 and `7.45e-9` in A21, so its finite-error
  result is not a second-Riccati artifact.

Two full-state storage hypotheses have now been falsified on that legal A21
word:

1. Raw physical-error moving-Riccati storage `z^T P^-1 z` first crosses
   `rho=1` at reliable scale 6.5 and reaches `rho=1.13702642918` on the tested
   declared direction.
2. Full measurement-linearizing storage
   `Phi(z)^T P^-1 Phi(z)`, with
   `epsilon_aw=(Q_aw-I)delta_a_w+e_eta` and all mixed transport retained, also
   first crosses `rho=1` at scale 6.5 and reaches `rho=1.1339495182`.

At scale 6.5 the limiting A21 direction has approximately
`||delta theta||=0.0235311 rad = 1.348 deg`,
`||delta b_a||=0.0861747 m/s^2`, and
`||delta a_w||=0.0121997 m/s^2`. The failure is therefore not caused by the
30-degree attitude candidate and cannot be repaired by selecting 25/20/15
degrees. H18 remains pointwise contractive throughout its tested declared
limiting direction; its worst reliable full-Phi ratio is `0.999973833561`.

Binary32 subtraction below about scale `0.008` is not used to judge the tangent
map. At resolved scale `0.125`, both raw and full-Phi central ratios recover the
linear complete-word tangent to about `1.3e-5` in H18 and `1e-4` in A21.

## Failure analysis

Classification: **proof-method failure for the current moving-Riccati storage
architectures**, not yet a theorem failure and not a source/CI/filter failure.

Failed controlling inequalities on one legal A21 complete word are

`V_raw(F_W(z)) < V_raw(z)`

and

`V_Phi(F_W(z)) < V_Phi(z)`.

Both are false for admissible finite errors beginning at scale 6.5 along the
reported full 21-state direction.

This invalidates the strategy of proving P4 by sharpening remainder bounds
around either of those two storages. It does **not** invalidate frozen P3, the
linear UES certificate, the exact residual/Joseph/reset identities, complete
SEA3 source generation, or the possibility of nonlinear UES in another
full-rank metric or over a longer finite window.

**DEAD_ENDS**

- duplicated pre-commit host observers;
- standalone/pure `e_eta` or packet-count remainder budgets;
- raw moving-Riccati endpoint storage as the nonlinear P4 certificate;
- full-Phi moving-Riccati endpoint storage as the nonlinear P4 certificate;
- Schur/elimination of `a_w` or any other active state from the final P4
  coercivity condition. Such a reduced quantity is not the required full-state
  theorem and must not be pursued as a rescue.

## Independent critic pass

Strongest reason to abandon the current architecture: the A21 obstruction is
not interval pessimism. A covariance-free nonlinear shadow using the exact
shipping source/gain word produces `rho>1` in both candidate storages well
inside the declared domain. Tightening local bounds around either storage is
therefore optimizing a false controlling inequality.

Qualitatively different full-state alternatives are:

1. **Longer complete-SEA3 differential window in the retained moving metric.**
   UES requires contraction over some bounded source-contiguous window, not
   necessarily the 600-sample diagnostic window. Test the full nonlinear
   Jacobian over 3 s / 6 s / 9 s while retaining all intervening shipping
   operations and the actual endpoint Riccati metrics.
2. **Full-state differential contraction / Riemannian metric.** If no longer
   moving-Riccati window is strict but the differential cocycle is uniformly
   stable, construct a full-rank state/source-dependent metric. No state may be
   eliminated and all actual-R_S operations remain in the word.
3. **Full-state converse/path-memory Lyapunov construction.** Build a
   source/state-indexed finite-horizon pullback metric from the complete
   Jacobian cocycle if a simple one-window metric is insufficient. This remains
   full rank and cannot be replaced by endpoint block elimination.
4. **Theorem falsification.** If the full-state nonlinear differential cocycle
   cannot be made uniformly contractive over bounded complete source windows
   while retaining uniformly coercive/bounded metrics, report that the
   declared nonlinear P4/UES formulation is unsupported on the current domain.

## Master inequality before new proof code

For a source state/history `s`, nonlinear active-state error `z`, and complete
word map `F_{W,s}`, a full-rank differential metric `M(s,z)>0` must satisfy

`DF_{W,s}(z)^T M(s+,F_{W,s}(z)) DF_{W,s}(z) <= rho M(s,z)`

with one uniform `rho<1` over every admitted complete SEA3 word and every state
in the certified nonlinear region, separately for H18 and A21. The H->A
21x18 dimension-changing event remains a separate hybrid obligation.

Because the system and admissible Lyapunov metric are time/source varying, the
spectral radius of one finite transition matrix is **not** a standalone no-go
criterion. The immediate feasibility quantity for the retained moving metric is

`gamma_W(z)^2 = || P_1^(-1/2) DF_W(z) P_0^(1/2) ||_2^2`.

If a bounded source-contiguous window has `sup gamma_W(z)^2 < 1`, then that
window is a valid full-state differential-contraction candidate in the retained
moving metric. If 3 s fails but 6 s or 9 s is strict, no new metric architecture
is needed. If all tested bounded windows fail, a different uniformly coercive
full-rank metric may still exist; its feasibility must be assessed through the
full differential cocycle, not through the spectral radius of one time-varying
transition.

No new nonlinear lemma is authorized unless it enters this matrix inequality
and has a quantified path to moving its largest generalized eigenvalue below 1.

## Retained facts

- canonical complete SEA3 same-history source and scheduler;
- frozen conditional P3 at `1e-18`;
- exact shipping H/P/R/K/S Joseph identities and full posterior precision;
- exact full shift `epsilon_aw=(Q_aw-I)delta_a_w+e_eta`;
- exact correction/reset and prediction shift transport, including all
  `F E_aw` v/p/S/a_w rows;
- actual applied R_S on every due S update;
- separate H18->A21 hybrid event;
- existing outward interval-AD Cayley/deployed-quaternion primitives.

## Next falsifiable experiment

Build one host-only **full-state differential Jacobian diagnostic**, driven by
the same single shipping estimator and complete SEA3 source. It must not create
a second covariance, tuner, source, or acceptance history. Propagate central
perturbation pairs for all active coordinates in parallel using higher-precision
shadow states and the frozen shipping gains/branch sequence, then recover the
full Jacobian.

For the A21 limiting source phase evaluate the origin and finite states on the
same limiting direction at scales 4, 6, 6.5, and 8. For each center evaluate
source-contiguous 3 s, 6 s, and 9 s windows and report:

- the complete 21x21 Jacobian and all active coordinates;
- `gamma_W^2`, the largest generalized moving-metric differential ratio;
- ordinary singular/eigen diagnostics only as conditioning information, not as
  theorem criteria;
- central-difference consistency under step refinement;
- endpoint P0/P1 conditioning;
- exact prediction/accelerometer/S/vector counts and a deterministic hash of
  the actual-R_S event sequence.

Also check H18 on its legal 3 s limiting word so both active dimensions remain
attached to the same diagnostic machinery. If a longer A21 window makes
`gamma_W^2` clearly below one, rigorous complete-SEA3 interval-AD enclosure of
that **full state and full window** is quantitatively justified. If no tested
window is strict, perform another architecture review before introducing a new
metric; do not return to endpoint storage transformations or eliminated-state
certificates.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
