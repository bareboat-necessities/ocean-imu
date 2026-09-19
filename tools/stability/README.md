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

## Measuring one H18 service superword

    make -C tools/stability h18-superword OUTPUT_DIR=/tmp

`ou3_theorem/h18_superword_export.cpp` drives the shipping filter through its
own startup and deployed handoff, keeps the accelerometer bias held with the
shipping external hold, and exports one magnetically informed service superword
rooted at an actually reached state: the covariance at every prefix, central
differences of the complete closed-loop finite-error map taken through the
shipping code, and the actual innovation covariance and sensitivity of every
correction the update really applied.

`ou3_theorem/h18_superword.py` evaluates the candidate storage `V(e)=e^T P^-1 e`
on that export in `decimal`, reporting the worst admissible ratio, its limiting
direction, every-prefix retention, the applied magnetic information and the
held-bias non-contraction obstruction. It is a feasibility diagnostic: it
measures one execution, sets no obligation, and cannot become a certificate.
The measured numbers are recorded in `docs/ou3-proof-research-state.md`.

The committed theorem status is fail-closed. Capture, finite-error
dissipativity, release retention, recurring magnetic information, every-prefix
retention, physical constant qualification, and finite-precision closure remain
open until certified.
