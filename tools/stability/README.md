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

    make -C tools/stability h18-iss OUTPUT_DIR=/tmp

The target runs one shipping-history export, the full-state obstruction
measurement, and the H18 complement/input diagnostic on that **same export**.
Use `h18-superword` to run only the first two steps.

`ou3_theorem/h18_superword_export.cpp` starts the shipping filter through its
own startup and handoff, using the existing external bias hold. The V2 export
records covariance and central differences at every prefix. Applied magnetic
sensitivities and transported error responses are captured before correction;
the whitening covariance is the actual one factored by that correction. A
changed applied-event sequence in the perturbed runs refuses the measurement.

`h18_superword.py` evaluates the full local incremental map at 60 Decimal
digits. Its prefix profile follows only the endpoint-maximizing direction and
is named accordingly. `h18_iss.py` separately evaluates **all directions at
every prefix** on the 18-dimensional complement using a complete symmetric
Jacobi eigensystem. It reports the fine/coarse discrepancy, available norm
margin, conditional Young multiplier and isolated linear bias Schur gain.

The exact conditional relation is `x_plus=A x+B e_b+r`. Local central
differences do not bound `r`, which includes reference forcing and nonlinear,
model and arithmetic remainders. Increasing evaluation precision does not
restore information lost in single-precision shipping differences. An observed
fine/coarse discrepancy is not a rigorous uncertainty bound.

Old post-correction magnetic exports, malformed measurements and incomplete
held histories are rejected. CI requires valid measurements but does not gate
on their ratios. Reports keep certificate, source-uniform and obligation flags
false, whether the numerical ratio is below or above one. Research findings
are recorded in `docs/ou3-proof-research-state.md`.

Capture, finite-error dissipativity, release retention, recurring magnetic
information, every-prefix retention, physical constant qualification, and
finite-precision closure remain open until certified.
