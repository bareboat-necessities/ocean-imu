# AGENTS.md

## Repository expectations

- Keep changes minimal and targeted.
- Do not change public estimator APIs or deployed behavior merely to simplify analysis.
- Treat the shipping OU-III implementation as authoritative for estimator behavior.

## Build and validation

- Build data comes from the oceanography-waves-lib release used by .github/workflows/build.yml.
- Primary validation command: `make all`.
- After changing C/C++ headers, sources, Makefiles, or deployed estimator behavior, run `make all`.
- If the build fails, report the exact failing command and error.
- Prefer vendored Eigen when present; otherwise preserve the `/usr/include/eigen3` fallback.

## Pull-request hygiene

Remove files, tests, workflows, generated evidence, and documentation made obsolete by the change. Documentation describes current behavior; chronological change history belongs in PR metadata.

## OU-III stability research protocol

There is one stability architecture. The three principal physical/deployment assumptions are MARINE MOTION, IMU BIAS, and MAGNETIC SERVICE, applied simultaneously to one persistent physical execution. Sensor-noise, model-error, arithmetic, and finite-precision premises may be stated separately when required.

The proof path is

`construction -> startup/capture -> magnetically informed Live/H18 -> H18-to-A21 release -> magnetically informed A21 -> regional practical stability`.

A certified tail inherits the actual estimator state, covariance, bias estimates, physical bias histories, frontend/tuner state, committed parameters, magnetic reference/state, scheduler, clocks, and physical source state. Do not reseed the estimator or restart physical-history coordinates at proof-word boundaries.

### Physical-history rules

- Vessel displacement is the wave coordinate about a local equilibrium/reference. Position, velocity, acceleration, attitude, and angular rate come from one continuous history.
- The displacement primitive is globally bounded on every admitted continuation. A permanent nonzero displacement DC component is inadmissible; quiet water is admissible.
- A local-equilibrium decomposition must not remove physical acceleration seen by the IMU.
- Accelerometer and gyroscope residual biases are separate bounded, rate-bounded physical histories. Their sampled successors are constrained by their predecessors.
- The estimator bias prior is not a law imposed on physical truth. Carry exact estimator/model mismatch.
- For accelerometer-bias prediction use one mode-dependent shipping relation: phi_e=1 in H18/held prediction and phi_e=phi_OU in A21/active prediction. Carry measurement correction and estimate projection as separate literal operations.
- Magnetic service is based only on actually applied, informative magnetic corrections. Calls, due events, packets, rejected measurements, invalid values, or saturated values do not establish service.

### Failure analysis

After a mathematical, enclosure, conditioning, numerical, implementation, CI, or infrastructure failure, update `docs/ou3-proof-research-state.md` with the exact failed quantity, failure classification, invalidated hypothesis, retained facts, current limiter, and next falsifiable experiment.

A tactic gets one implementation and at most one mathematically motivated refinement. A second failure of the same mechanism requires an architecture review. Do not repeat subdivision, interval refinement, tighter scalar norms, or deeper search without a quantitative argument that the change can cross the controlling threshold.

### Controlling theorem inequality

Before implementing a new lemma, state where it enters the finite-error tail inequality. A typical service-superword target is

`V_(j+1) <= rho V_j + c_d ||d||^2, rho < 1`,

together with coercivity and every-prefix retention on the same physical history. Local covariance, reset, measurement, or arithmetic lemmas are subordinate to that complete finite-error statement.

Run a non-promoting high-precision feasibility diagnostic before rigorous enclosure of a new contraction/storage construction. If its worst admissible ratio is above one, change the formulation instead of sharpening unrelated bounds.

### Research ledger

Keep `docs/ou3-proof-research-state.md` concise and current under Current hypothesis, Evidence, Current limiter, Failed approaches / DEAD_ENDS, Retained facts, Alternatives, and Next falsifiable experiment.

Do not claim the end-to-end theorem until finite capture, finite-error dissipativity, retention, and implementation/arithmetic premises are all closed.
