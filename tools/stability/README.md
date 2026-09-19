# OU-III stability tooling

The executable package contains only contracts and algebra used by the single
OU-III theorem. Shipping estimator behavior remains authoritative.

The proof route is
`construction -> capture -> finite H18 bridge -> finite reference refinement
and release -> recurring magnetically informed A21 -> regional practical
stability`.

PR #557 supplied a negative research result: strict full-state H18 incremental
contraction is obstructed by the held accelerometer-bias identity block. Its
trajectory export and finite-difference diagnostics are intentionally not
retained as proof machinery.

Current theorem tooling covers MARINE MOTION, IMU BIAS, MAGNETIC SERVICE,
same-execution linkage, the no-heading-service necessity result, and analytic
finite-bridge/A21 small-gain lemmas. No sampled trajectory can promote a proof
obligation.

Run the theorem contracts with:

    cd tests/validation
    python3 -m unittest -v test_ou3_architecture_cleanup test_ou3_imu_bias       test_ou3_magnetic_service test_ou3_marine_motion       test_ou3_no_mag_obstruction test_ou3_same_execution       test_ou3_theorem_status test_ou3_tail_stability

The controlling open work is finite magnetic-reference refinement/release,
source-uniform A21 linear dissipation, and a nonlinear finite-error remainder
bound small enough to satisfy `sqrt(rho0)+eta < 1`.
