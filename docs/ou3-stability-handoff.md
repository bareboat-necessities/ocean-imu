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
- The H18 hold is exact: the held accelerometer-bias estimate, its covariance block and its zeroed cross-covariances are bit-for-bit unchanged across a magnetically served window. On feasible held segments without projection-changing or frame/relock events, the same-history incremental bias coordinate is invariant. Absolute physical bias error still carries its predecessor increment. Strict full-state incremental contraction is unavailable there; the absolute bias is bounded under the projection lemma premises.

## Open controlling obligations

1. Qualify assembled-sensor bias/residual limits, including the accelerometer bias-rate limit.
2. Supply an all-time marine bounded-primitive certificate.
3. Prove finite history-dependent startup/capture into the finite-error tail domain.
4. Certify the H18 complement matrix inequality and total finite-error supply, including held bias and nonlinear/reference forcing, with every-prefix retention. The conditional supply algebra and complement feasibility diagnostic are implemented; their source-uniform premises remain open.
5. Prove the actual H18-to-A21 transition retains the certified domain.
6. Close the corresponding magnetically informed A21 finite-error inequality.
7. Prove recurring MAGNETIC SERVICE using actually applied sensitivities.
8. Close shipping finite-precision/arithmetic totality and all hard-event retention.

No end-to-end stability claim is authorized while any item remains open.

## Reproduction

`cd tests/validation && python3 -m unittest -v test_ou3_architecture_cleanup test_ou3_imu_bias test_ou3_magnetic_service test_ou3_marine_motion test_ou3_no_mag_obstruction test_ou3_theorem_status test_ou3_finite_error test_ou3_h18_superword test_ou3_h18_iss`

`python3 tools/stability/ou3_theorem/build_evidence.py --output /tmp/ou3-stability-evidence.json`

`make -C tests/kalman_ou_iii shipping_contract-test shipping_transition-test && ./tests/kalman_ou_iii/shipping_contract-test && ./tests/kalman_ou_iii/shipping_transition-test`

## Measuring one superword

    make -C tools/stability h18-iss OUTPUT_DIR=/tmp

The exporter runs shipping startup and the deployed handoff, preserves one
physical history, and holds bias using the existing external control. V2
records each applied magnetic sensitivity and error response immediately
**before** correction, with the actual innovation covariance. Perturbed runs
must have the same applied-event sequence. Old post-correction exports are
refused, not silently reinterpreted.

`h18_superword.py` measures the full incremental map. Its endpoint-direction
prefix profile is labeled accordingly. `h18_iss.py` uses the same export for
all-direction, every-prefix complement gains, fine/coarse conditioning, the
conditional Young multiplier and the isolated linear bias Schur gain. Decimal
Jacobi evaluation resolves the whole symmetric spectrum, not a single power
iteration direction. None is a finite-error or source-uniform certificate.
CI gates valid production of measurements and contracts, never the numerical
contraction result. Raw exports and both reports are uploaded together.

## Latest continuation

At the one-second word rooted 30 seconds after Live, the corrected magnetic
information minimum is 1.707892 (floor 1). The complement endpoint ratio is
0.999803093, but the whitened fine/coarse map discrepancy 0.0140108 is about
142 times the available perturbation margin. Both scales have endpoint ratios
below one; their disagreement is not a rigorous error bound or an instability
result. The all-direction prefix maximum is 1.001023761 at sample 5.

The exact conditional supply lemma is stated in `docs/ou3-stability-proof.md`.
For the diagnostic rho 0.999901546, the isolated linear bias gain is about
1.690e6; the total-supply Young multiplier is about 1.016e4. The nonlinear and
reference forcing remain unbounded in this storage. These are not practical
error bounds and do not authorize promotion of any theorem obligation.

The next falsifiable experiment is a branch-consistent sensitivity/remainder
calculation about the same inherited execution, with all auxiliary states
carried. Its error budget must be compared quantitatively with 9.846e-5 in the
whitened operator norm. Do not pursue blind step-size or interval refinement;
change the storage construction if a reliable matrix bound is not below one.

Run the shared workflow and fingerprint regressions alongside the theorem tests:
`cd tests/validation && python3 -m unittest -v test_workflow_contract test_ou_replay_fingerprint`.
The native transition regressions remain unchanged. Mainline estimator source,
quality gates and physical certification constants are unchanged.
