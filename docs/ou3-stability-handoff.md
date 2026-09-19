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
- The H18 hold is exact: the held accelerometer-bias estimate, its covariance block and its zeroed cross-covariances are bit-for-bit unchanged across a magnetically served window. The held bias is therefore an invariant coordinate of the H18 superword map, and strict full-state contraction is unavailable on that leg. It is bounded instead, by the estimate-projection lemma.

## Open controlling obligations

1. Qualify assembled-sensor bias/residual limits, including the accelerometer bias-rate limit.
2. Supply an all-time marine bounded-primitive certificate.
3. Prove finite history-dependent startup/capture into the finite-error tail domain.
4. Restate the H18 obligation off the held accelerometer-bias coordinates -- strict dissipation on the complement, with the held bias as a bounded input -- and prove the restated same-history inequality with every-prefix retention. The full-state form is refuted, not merely unproved.
5. Prove the actual H18-to-A21 transition retains the certified domain.
6. Close the corresponding magnetically informed A21 finite-error inequality.
7. Prove recurring MAGNETIC SERVICE using actually applied sensitivities.
8. Close shipping finite-precision/arithmetic totality and all hard-event retention.

No end-to-end stability claim is authorized while any item remains open.

## Reproduction

`cd tests/validation && python3 -m unittest -v test_ou3_architecture_cleanup test_ou3_imu_bias test_ou3_magnetic_service test_ou3_marine_motion test_ou3_no_mag_obstruction test_ou3_theorem_status test_ou3_finite_error test_ou3_h18_superword`

`python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json`

`make -C tests/kalman_ou_iii shipping_contract-test shipping_transition-test && ./tests/kalman_ou_iii/shipping_contract-test && ./tests/kalman_ou_iii/shipping_transition-test`

## Measuring one superword

    make -C tools/stability h18-superword OUTPUT_DIR=/tmp

`tools/stability/ou3_theorem/h18_superword_export.cpp` runs the shipping filter
through its own startup and deployed handoff, keeps H18 with the shipping
external hold, and exports one magnetically informed service superword: the
covariance at every prefix, central differences of the complete closed-loop
finite-error map taken through the shipping code, and the actual innovation
covariance and sensitivity of every correction the update really applied.
`h18_superword.py` evaluates the candidate storage on that export in `decimal`
at 60 digits and reports the worst admissible ratio, its limiting direction,
every-prefix retention, applied magnetic information and the held-bias
obstruction. Neither step gates CI and neither can discharge an obligation; the
measured numbers live in `docs/ou3-proof-research-state.md`.

## Latest continuation

The full-state H18 dissipation target is refuted by the held-bias obstruction,
so the next measurement is the complement map described under Next falsifiable
experiment in the research state. Do not spend further enclosure effort on the
full-state form.

Run the shared workflow and fingerprint regressions alongside the theorem tests:
`cd tests/validation && python3 -m unittest -v test_workflow_contract test_ou_replay_fingerprint`.
The native shipping-transition test also checks actual gain correction,
nonpositive projection radius, covariance preservation at projection, the
invalid-injection branch that bypasses projection, and the exactness of the H18
hold across a served window. Source-pinned CI archives and logs are retained
even when a contract fails.
