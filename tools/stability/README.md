# OU-III stability tooling

The repository has one OU-III stability theorem architecture. The shipping
implementation in `src/kalman_ou_iii/` is authoritative for estimator
behavior.

MARINE MOTION, IMU BIAS, and MAGNETIC SERVICE apply simultaneously to one
persistent physical execution. The executable theorem package is
`tools/stability/ou3_theorem/`.

The accelerometer-bias prediction relation is shared by held and active modes:
the estimator multiplier is one while held and the literal OU multiplier while
active. Measurement correction, estimate projection, release, and attitude
operations are composed separately using the shipping maps.

The proof path is
`construction -> startup/capture -> magnetically informed Live/H18 -> H18-to-A21 release -> magnetically informed A21 -> regional practical stability`.

The committed theorem status is fail-closed. Capture, finite-error
dissipativity, release retention, recurring magnetic information, every-prefix
retention, physical constant qualification, and finite-precision closure remain
open until certified.
