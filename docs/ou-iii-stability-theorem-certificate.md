# OU-III analytical stability theorem certificate

This note records the machine-checkable closure of the OU-III analytical stability argument. It is intentionally separate from the stronger numerical source-funnel/deployment certificate described in `docs/ou-iii-numerical-certificate.md`.

## Result

For the configured 200 Hz validation/runtime timing contract, the current adaptive OU-III implementation has a source-bound **conditional branch-regular normal-Live local ISS certificate** in both fixed-dimensional Live modes:

- held accelerometer-bias mode `H`, dimension 18;
- active accelerometer-bias mode `A`, dimension 21.

The full-heading result is conditional on the explicit persistent-excitation packet envelope proved by `ou3_vector_uco_certificate.py`: the proof packet uses accepted accelerometer and magnetometer vectors, with two consecutive accepted magnetic packets at the configured 25 Hz spacing, nonzero vector magnitudes/separation, and the stated rate bound. Arbitrarily long rejection cannot support an unconditional full-heading theorem.

The nonlinear lift is also explicitly **branch regular**. The nominal word must have positive margin to innovation/gating thresholds so a sufficiently small error neighborhood follows the same finite source branch word. A nominal point exactly on a gate discontinuity is outside this local theorem. Gravity-only operation remains a separate yaw-quotient result in the manuscript.

These qualifications are mathematical hypotheses, not values inferred from the eight replay trajectories.

The producer is `tools/ou3_stability_theorem_certificate.py`. It regenerates every upstream source-bound lemma from the current implementation and returns `PASS_CONDITIONAL_LOCAL_ISS` only when all machine-checkable obligations pass and the non-source hypotheses are stated explicitly.

## Closed proof chain

The certificate composes the following current-source obligations rather than replay-derived minima:

1. **Compact source domain.** `ou3_source_domain_contract.py` extracts deployed binary32 constants, adaptation clamps, timing constants, source branches, and an outward-rounded continuous parameter box.
2. **Exact scalar OU transition.** `ou3_scalar_ou_enclosure.py` encloses the shipping one-step OU/integrator transition with validated transcendental arithmetic.
3. **Translational UCO/UCC.** `ou3_translational_uco_ucc.py` supplies a strict bounded-gap observation/detectability certificate for `(v,p,S,a_w)` and a source-uniform process-excitation lower bound.
4. **Attitude/gyro-bias vector UCO.** `ou3_vector_uco_certificate.py` supplies a strict information lower bound under the explicit accepted-packet persistent-excitation envelope. Collinear or indefinitely rejected vector histories are not silently certified.
5. **Complete process UCC.** `ou3_full_process_ucc.py` gives strict one-step process-covariance lower bounds for both H and A modes.
6. **Stable Gauss-Markov tails.** The OU acceleration state and active accelerometer-bias state have strict discrete decay factors below one; the latter is independently enclosed with the same validated exponential backend.
7. **Source Gaussian primitive model.** `ou3_source_noise_certificate.py` derives the finite standardized primitive-noise covariance from source/configuration rather than truth-error replay.
8. **Periodic `a_w` covariance synchronization.** `ou3_hybrid_aw_sync_proof.py` proves the shipping synchronization is PSD covariance inflation. Therefore `P+ >= P- > 0`, `P+^{-1} <= P-^{-1}`, and information energy is nonexpansive across that jump.

These establish the hypotheses used by the manuscript's discrete LTV Kalman/Riccati stability argument: bounded source coefficients, uniform detectability on the stated accepted-packet operating envelope, uniform complete controllability/stabilizability, and positive finite measurement/process covariances. Hence the stabilizing Riccati family is uniformly bounded and the homogeneous normal-Live H/A linearized recursions are uniformly exponentially stable.

On a branch-regular nominal source word and a sufficiently small geodesic chart, the selected MEKF prediction/correction/reset branches are `C^1`; their error relative to the linearization is uniformly quadratic on the compact source box. Uniform exponential stability plus the discrete variation-of-constants/small-gain argument therefore gives a nonzero invariant local practical-ISS neighborhood.

## What `PASS_CONDITIONAL_LOCAL_ISS` means

A PASS establishes, without fitting to the eight wave trajectories:

- H-mode linearized normal-Live UES;
- A-mode linearized normal-Live UES;
- existence of a nonzero nonlinear local ISS neighborhood on branch-regular source words;
- explicit dependence of full-heading stability on accepted consecutive vector-packet persistent excitation;
- explicit exclusion of innovation-gate boundary points from the local nonlinear theorem;
- exact nonexpansiveness of periodic `a_w` covariance synchronization;
- source-code binding through hashes of the implementation headers used by the proof.

The certificate deliberately does **not** fabricate a numerical basin radius from sampled perturbation searches.

## Separate stronger deployment-funnel claim

The analytical local-ISS theorem is now closed. A different, stronger question remains available in the numerical-certificate pipeline: produce an explicit outward-rounded source-funnel radius together with startup capture, all hard reset/regauge inequalities, finite-horizon stochastic concentration, and finite-capture constants.

That deployment claim requires physical operating-envelope quantities that estimator source cannot bound by itself, including startup/re-lock non-gravitational specific force, accepted magnetic excitation/rejection gaps, heading-regauge magnitude, and bounded physical world-vector magnitudes. Those assumptions must be supplied by the intended deployment environment; treating eight simulated seas as universal physical bounds would invalidate the theorem.

Accordingly the repository keeps two statuses separate:

- **analytical branch-regular normal-Live stability theorem:** `PASS_CONDITIONAL_LOCAL_ISS` when the composed source-bound proof passes;
- **explicit numerical source-funnel/deployment theorem:** remains governed by `ou3_validate_enclosure.py` and `ou3_deployment_gate.py` and may report `NOT_ESTABLISHED` until a declared physical deployment envelope is supplied and validated.

The second claim is strictly stronger than the first and is not an unfinished algebraic step in the local-ISS proof.

## CI

`ou3-proof-fast` contains a dedicated `stability-theorem` job. It recompiles the producer and its tests, regenerates the certificate against the current implementation, and fails if an upstream proof obligation, accepted-packet/branch-regular qualification, source binding, strict stable tail, H/A UCC bound, or claim-separation invariant is lost.
