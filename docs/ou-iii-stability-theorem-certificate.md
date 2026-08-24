# OU-III analytical stability theorem certificate

This note records the machine-checkable closure of the OU-III analytical stability argument. It is intentionally separate from the stronger numerical source-funnel/deployment certificate described in `docs/ou-iii-numerical-certificate.md`.

## Result

For the configured 200 Hz validation/runtime timing contract, the current adaptive OU-III implementation has a source-bound **conditional normal-Live local ISS certificate** in both fixed-dimensional Live modes:

- held accelerometer-bias mode `H`, dimension 18;
- active accelerometer-bias mode `A`, dimension 21.

The result is conditional on the explicit persistent-excitation operating envelope required for full-heading observability: accepted gravity/specific-force and magnetic vector packets must retain nonzero vector magnitude/separation and a bounded accepted-packet gap. This qualification is mathematical, not empirical. Arbitrarily long magnetometer rejection cannot support an unconditional full-heading theorem; gravity-only operation is instead interpreted on the yaw quotient as in the manuscript.

The producer is `tools/ou3_stability_theorem_certificate.py`. It regenerates every upstream source-bound lemma from the current implementation and returns `PASS_CONDITIONAL_LOCAL_ISS` only when all obligations pass.

## Closed proof chain

The certificate composes the following current-source obligations rather than replay-derived minima:

1. **Compact source domain.** `ou3_source_domain_contract.py` extracts deployed binary32 constants, adaptation clamps, timing constants, source branches, and an outward-rounded continuous parameter box.
2. **Exact scalar OU transition.** `ou3_scalar_ou_enclosure.py` encloses the shipping one-step OU/integrator transition with validated transcendental arithmetic.
3. **Translational UCO/UCC.** `ou3_translational_uco_ucc.py` supplies a strict bounded-gap observation/detectability certificate for `(v,p,S,a_w)` and a source-uniform process-excitation lower bound.
4. **Attitude/gyro-bias vector UCO.** `ou3_vector_uco_certificate.py` supplies a strict information lower bound under the explicit vector persistent-excitation envelope. Collinear or indefinitely rejected vector packets are not silently certified.
5. **Complete process UCC.** `ou3_full_process_ucc.py` gives strict one-step process-covariance lower bounds for both H and A modes.
6. **Stable Gauss-Markov tails.** The OU acceleration state and active accelerometer-bias state have strict discrete decay factors below one; the latter is independently enclosed with the same validated exponential backend.
7. **Source Gaussian primitive model.** `ou3_source_noise_certificate.py` derives the finite standardized primitive-noise covariance from source/configuration rather than truth-error replay.
8. **Periodic `a_w` covariance synchronization.** `ou3_hybrid_aw_sync_proof.py` proves the shipping synchronization is PSD covariance inflation. Therefore `P+ >= P- > 0`, `P+^{-1} <= P-^{-1}`, and the information energy is nonexpansive across that jump.

These establish the hypotheses used by the manuscript's discrete LTV Kalman/Riccati stability argument: bounded source coefficients, uniform detectability on the stated operating envelope, uniform complete controllability/stabilizability, and positive finite measurement/process covariances. Hence the stabilizing Riccati family is uniformly bounded and the homogeneous normal-Live linearized H/A recursions are uniformly exponentially stable.

On any sufficiently small geodesic chart that does not cross a source branch boundary, the fixed-branch MEKF prediction/correction/reset maps are piecewise `C^1`; their error relative to the linearization is uniformly quadratic on the compact source box. UES plus the discrete variation-of-constants/small-gain argument therefore gives a nonzero invariant local practical-ISS neighborhood for the nonlinear normal-Live recursion.

## What `PASS_CONDITIONAL_LOCAL_ISS` means

A PASS means all of the following are established without fitting to the eight wave trajectories:

- H-mode linearized normal-Live UES;
- A-mode linearized normal-Live UES;
- existence of a nonzero nonlinear local ISS neighborhood in each fixed mode;
- explicit dependence of full-heading stability on the persistent-excitation hypothesis;
- exact nonexpansiveness of the periodic `a_w` covariance synchronization jump;
- source-code binding through hashes of the implementation headers used by the proof.

The certificate deliberately does **not** fabricate a numerical basin radius from sampled perturbation searches.

## Separate stronger deployment-funnel claim

The analytical local-ISS theorem is now closed. A different, stronger question remains available in the numerical-certificate pipeline: produce an explicit, outward-rounded source-funnel radius together with startup capture, all hard reset/regauge inequalities, finite-horizon stochastic concentration, and finite-capture constants.

That deployment claim requires physical operating-envelope quantities that the estimator source cannot bound by itself, including startup/re-lock non-gravitational specific force, accepted magnetic excitation/gap, heading-regauge magnitude, and bounded physical world-vector magnitudes. Those assumptions must be supplied explicitly by the intended deployment environment; treating eight simulated seas as universal bounds would invalidate the theorem.

Accordingly the repository keeps these two statuses separate:

- **analytical normal-Live stability theorem:** `PASS_CONDITIONAL_LOCAL_ISS` when the composed source-bound proof passes;
- **explicit numerical source-funnel/deployment theorem:** remains governed by `ou3_validate_enclosure.py` and `ou3_deployment_gate.py` and may report `NOT_ESTABLISHED` until a declared physical deployment envelope is supplied and validated.

This separation is intentional: the second claim is strictly stronger than the first and is not an unfinished algebraic step in the local-ISS proof.

## CI

`ou3-proof-fast` contains a dedicated `stability-theorem` job. It recompiles the producer and its tests, regenerates the certificate against the current implementation, and fails if any upstream proof obligation, persistent-excitation qualification, source binding, strict stable tail, H/A UCC bound, or claim-separation invariant is lost.
