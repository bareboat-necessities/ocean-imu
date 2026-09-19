# OU-III stability handoff

## Architecture

There is one theorem under simultaneous MARINE MOTION, IMU BIAS, and MAGNETIC
SERVICE assumptions on one persistent physical execution. Shipping code is
authoritative.

The proof route is
`construction -> capture -> finite H18 bridge -> finite reference refinement
and bias release -> recurring magnetically informed A21 -> regional practical
stability`.

## Key correction in PR #558

H18 is no longer asked to provide asymptotic contraction. PR #557 showed why a
particular full-state held-bias contraction formulation cannot work. That
negative result is retained in prose, while its experimental export and
finite-difference scripts have been removed.

The shipping wrapper itself supplies the reason to move on: accelerometer-bias
learning is externally held while the magnetic reference is provisional, and
the hold is released when refinement completes. The internal gate additionally
requires its accepted-magnetometer-update threshold and a one-second guard.
Under recurring MAGNETIC SERVICE the count/guard part is finite once reference
refinement is finite.

## Productive proof obligations

1. Qualify the assembled sensor/bias limits and all-time marine-motion
   membership.
2. Prove finite history-dependent startup/capture.
3. Prove finite magnetic-reference refinement and hence finite H18 release.
4. Bound the finite H18 bridge and prove the literal release retains the A21
   tail domain.
5. Establish a source-uniform complete A21 information floor in normalized
   covariance coordinates. MAGNETIC SERVICE supplies the heading/axial
   gyro-bias component; gravity/accelerometer and integral pseudo-updates must
   close the remaining directions. A full floor mu gives the linear comparison
   rho0<=1/(1+mu).
6. Retain A21 inside the inactive accelerometer-bias projection region
   (bias-error radius <0.1748334 m/s^2), then bound the remaining nonlinear
   reset/tuner/arithmetic remainder. Close finite-error contraction with
   (sqrt(rho0)+eta)^2 < 1.
7. Close every-prefix retention, recurring service, and finite-precision
   arithmetic.

## Reproduction

`cd tests/validation && python3 -m unittest -v test_ou3_architecture_cleanup test_ou3_imu_bias test_ou3_magnetic_service test_ou3_marine_motion test_ou3_no_mag_obstruction test_ou3_same_execution test_ou3_theorem_status test_ou3_tail_stability`

`python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json`

`make -C tests/kalman_ou_iii shipping_contract-test shipping_transition-test && ./tests/kalman_ou_iii/shipping_contract-test && ./tests/kalman_ou_iii/shipping_transition-test`

No H18 numerical diagnostic is a proof gate. The next implementation work
should target finite reference refinement/release and the analytic A21
dissipativity inequality.
