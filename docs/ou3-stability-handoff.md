# OU-III stability handoff

## Current theorem

There is one architecture with three simultaneous principal assumptions: MARINE MOTION, IMU BIAS, and MAGNETIC SERVICE. Shipping implementation is authoritative and pinned by `reports/results/ou3_stability/provenance.json`.

## Closed structural facts

- Permanent nonzero wave-displacement DC is excluded by the all-time bounded-primitive contract; quiet water is admitted.
- Motion and bias samples are predecessor-linked parts of one persistent physical history.
- Physical accelerometer and gyro biases are independent bounded/rate-bounded total residuals.
- One accelerometer-bias prediction equation covers both estimator modes: `phi_e=1` in H18 and `phi_e=phi_OU` in A21. Correction and projection are separate literal shipping operations.
- An attempted magnetic callback is not the same as an applied correction. Informative service consumes only actually applied, valid, unsaturated, gauged corrections and their transported/whitened information.
- The indefinitely ungauged heading/axial-gyro-bias block is unipotent with spectral radius one and survives only as a necessity result.
- Live handoff and H18-to-A21 release inherit the existing execution and are not proof resets.

## Open controlling obligations

1. Qualify assembled-sensor bias/residual limits, including the accelerometer bias-rate limit.
2. Supply an all-time marine bounded-primitive certificate.
3. Prove finite history-dependent startup/capture into the finite-error tail domain.
4. Build a complete same-history magnetically informed H18 finite-error storage inequality with strict service-superword dissipation and every-prefix retention.
5. Prove the actual H18-to-A21 transition retains the certified domain.
6. Close the corresponding magnetically informed A21 finite-error inequality.
7. Prove recurring MAGNETIC SERVICE using actually applied sensitivities.
8. Close shipping finite-precision/arithmetic totality and all hard-event retention.

No end-to-end stability claim is authorized while any item remains open.

## Reproduction

`cd tests/validation && python3 -m unittest -v test_ou3_architecture_cleanup test_ou3_imu_bias test_ou3_magnetic_service test_ou3_marine_motion test_ou3_no_mag_obstruction test_ou3_theorem_status`

`python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json`

`make -C tests/kalman_ou_iii shipping_contract-test shipping_transition-test && ./tests/kalman_ou_iii/shipping_contract-test && ./tests/kalman_ou_iii/shipping_transition-test`

## Latest continuation

Keep all work on PR #553, `codex/ou3-single-stability-architecture`. The duplicate
PR #554 is closed. Do not reopen it or introduce another working branch.

The failing cleanup gate was reproduced locally and its seven file-level
causes removed without an exclusion. The performance replay fingerprint is
explicitly invalidated; only a new validated full evidence run can replace it.
New finite-error lemmas retain the joint true-bias/error prediction, separate
projection defect, conditional Euclidean projection sector, fixed-origin
integral innovation, and rectangular information transport. None closes the
complete source-uniform H18/A21 storage inequality.

Run the shared workflow and fingerprint regressions alongside the theorem tests:
`cd tests/validation && python3 -m unittest -v test_workflow_contract test_ou_replay_fingerprint`.
The native shipping-transition test also checks actual gain correction,
nonpositive projection radius, covariance preservation at projection, and the
invalid-injection branch that bypasses projection. Source-pinned CI archives
and logs are retained even when a contract fails.
