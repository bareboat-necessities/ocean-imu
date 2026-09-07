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

1. **Full-state differential contraction / Riemannian metric.** Enclose the
   exact 21x21 Jacobian of the nonlinear complete word and use a full-rank
   state/source-dependent metric. No state is eliminated and all actual-R_S
   operations remain in the word.
2. **Longer complete-SEA3 finite window.** UES requires contraction over some
   bounded window, not necessarily the 600-sample diagnostic window. Compose
   consecutive complete source words with all intervening shipping events and
   test 6 s / 9 s full-state Jacobians before changing metric architecture.
3. **Full-state converse/path-memory Lyapunov construction.** If longer exact
   Jacobian products are uniformly stable but no simple one-window metric is,
   construct a source/state-indexed finite-horizon pullback metric from the
   full Jacobian products. This remains full rank and is not an endpoint state
   elimination.
4. **Theorem falsification.** If an admissible full-state nonlinear Jacobian has
   spectral radius >=1 persistently over bounded longer windows, stop proof
   construction and report that the declared nonlinear P4/UES formulation is
   not supported on the current domain.

## Master inequality before new proof code

For a source state/history `s`, nonlinear active-state error `z`, and complete
word map `F_{W,s}`, a full-rank differential metric `M(s,z) > 0` must satisfy

`DF_{W,s}(z)^T M(s+,F_{W,s}(z)) DF_{W,s}(z) <= rho M(s,z)`

with one uniform `rho<1` over every admitted complete SEA3 word and every state
in the certified nonlinear region, separately for H18 and A21. The H->A
21x18 dimension-changing event remains a separate hybrid obligation.

The immediate feasibility quantity is the spectral radius of the **full**
Jacobian `J_W(z)=DF_W(z)`. If `spectral_radius(J_W(z))>1`, no equivalent norm can
make that same complete-word map locally contractive at that state; metric
search for that word is then futile. If the 3 s Jacobian is locally expansive
but a longer full-source product is Schur stable, the theorem may still close
with a longer finite window.

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

Build one host-only **full-state Jacobian diagnostic**, driven by the same single
shipping estimator and complete SEA3 word. It must not create a second covariance
or source history. Propagate central perturbation pairs for all active
coordinates in parallel using higher-precision shadow states and the frozen
shipping gains/branch sequence, then recover the complete-word Jacobian.

For the A21 limiting word evaluate at least the origin and the finite states on
the same direction at scales 4, 6, 6.5, and 8. Report:

- all 21 singular/eigen directions, not a reduced block;
- spectral radius of the full Jacobian;
- largest singular value in the current shipping metric only as a diagnostic;
- determinant/conditioning and central-difference consistency under step
  refinement;
- exact event counts and actual-R_S history.

Then compose the same full-state diagnostic over 6 s and 9 s source-contiguous
windows beginning at the same phase. If spectral radius is below one with
numerical margin, a full-rank differential/converse metric is quantitatively
justified and only then should interval-AD enclosure be implemented. If it is
not below one on any bounded tested window, treat that as theorem-falsification
evidence rather than introducing another storage shortcut.

P4 remains OPEN. P5 remains BLOCKED. No merge is authorized.
