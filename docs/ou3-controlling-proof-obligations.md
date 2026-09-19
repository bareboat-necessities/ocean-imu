# OU-III controlling proof obligations — PR #558

This list is fail-closed. “Infrastructure closed” means the theorem machinery
needed to attack the obligation exists; it does not mean the theorem obligation
is discharged.

| # | Obligation | Current state | Controlling next certificate |
|---|---|---|---|
| 1 | Root covariance floor p_min | OPEN; interval infrastructure closed | Assemble literal shipping event cadence and verify a self-containing 21-state recurring covariance box |
| 2 | Covariance-normalized mu/rho0 | BLOCKED by 1 | mu_cov=p_min*2.04e-3; rho0=1/(1+mu_cov) only after 1 |
| 3 | Explicit nonlinear r_* | BLOCKED by 2 | Certified L2 on retained finite-error domain and r*=min(r_guard,(1-sqrt(rho0))/L2) |
| 4 | Whole-word float32 supply | OPEN | Literal per-kernel operation counts/magnitude boxes composed over the proof word using gamma_n |
| 5 | Finite startup/capture | OPEN | History-dependent capture into the explicit retained domain; no common deadline is assumed |
| 6 | Finite H18 bridge retention | OPEN | Prefix bound from capture through actual refinement/release time |
| 7 | Actual refinement/release retention | OPEN | Compose captured-domain tuner completion, 250-update/1-s gate, covariance release and state continuity |
| 8 | Capture/release into explicit A21 region | BLOCKED by 2–7 | Compare certified release set with explicit r_* / projection-sector retained set |
| 9 | Every-prefix tail retention | OPEN | Prefix gains for every operation inside the recurring A21 word |
| 10 | Recurring source-uniform magnetic service | EXTERNAL/ALL-TIME; schema closed | One-history continuation certificate for every service window; finite replay is insufficient |
| 11 | Physical sensor/bias/marine qualification | EXTERNAL/ALL-TIME; composition closed | Same-history MARINE MOTION + IMU BIAS + MAGNETIC SERVICE plus assembled sensor/mount/calibration qualification |
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
