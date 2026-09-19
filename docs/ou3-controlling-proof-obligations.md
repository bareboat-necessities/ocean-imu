# OU-III controlling proof obligations — PR #558

This list is fail-closed. “Infrastructure closed” means the theorem machinery
needed to attack the obligation exists; it does not mean the theorem obligation
is discharged.

| # | Obligation | Current state | Controlling next certificate |
|---|---|---|---|
| 1 | Root covariance floor p_min | OPEN; generic interval layer + event order + conservative maximal-correction order closed | Numerically verify a self-containing 21-state box; subdivide parameter/attitude branch cells if entrywise dependency prevents innovation inclusion |
| 2 | Covariance-normalized mu/rho0 | BLOCKED by 1 | mu_cov=p_min*2.04e-3; rho0=1/(1+mu_cov) only after 1 |
| 3 | Explicit nonlinear r_* | BLOCKED by 2 | Certified L2 on retained finite-error domain and r*=min(r_guard,(1-sqrt(rho0))/L2) |
| 4 | Whole-word float32 supply | OPEN; local-kernel composition closed | Obtain magnitude/prefix-gain boxes from the verified Riccati/tail enclosure, then insert literal per-kernel operation counts using gamma_n |
| 5 | Finite startup/capture | OPEN | History-dependent capture into the explicit retained domain; no common deadline is assumed |
| 6 | Finite H18 bridge retention | OPEN; capture→bridge→release set-inclusion composition closed | Insert actual capture storage, history-dependent bridge length, supply/gain and release map bounds |
| 7 | Actual refinement/release retention | OPEN | Compose captured-domain tuner completion, 250-update/1-s gate, covariance release and state continuity |
| 8 | Capture/release into explicit A21 region | BLOCKED by 2–7 | Compare certified release set with explicit r_* / projection-sector retained set |
| 9 | Every-prefix tail retention | OPEN; prefix-retention composition closed | Obtain certified prefix gains/additive bounds from the explicit tail enclosure |
| 10 | Recurring source-uniform magnetic service | EXTERNAL/ALL-TIME; continuation schema closed and tested | Supply one-history every-window certificate using actual applied-event/innovation semantics; finite replay is insufficient |
| 11 | Physical sensor/bias/marine qualification | EXTERNAL/ALL-TIME; same-history composition closed and tested | Supply MARINE MOTION + IMU BIAS + MAGNETIC SERVICE continuations plus assembled sensor/mount/calibration qualification |
| 12 | Final implementation/arithmetic totality | OPEN | Close every finite branch, innovation solve, projection/reset, scheduler counter and arithmetic enclosure |

## Interval-Riccati progress

The full generic 21-state interval layer now exists: source-uniform prediction
F/Q boxes, accelerometer/S/magnetometer H/R boxes, verified innovation inverse
and gain, literal Joseph update, all covariance cross terms, periodic a_w
covariance-floor event, H18-to-A21 bias covariance-floor event, and
self-containing recurring-box iteration.

The remaining work for obligation 1 is not another algebraic kernel. It is the
literal hybrid event-word assembly: encode the shipping sample/pseudo-update/
magnetic-service cadence and branch cells tightly enough that interval
dependency does not destroy the innovation inverse and recurring-box
inclusions. The certificate remains OPEN until an actual box verifies.
