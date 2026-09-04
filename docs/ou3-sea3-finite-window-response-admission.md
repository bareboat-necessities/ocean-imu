# OU-III SEA3 finite-window response admission

This follow-up starts from merged PR #482 and advances the remaining SEA0 left-inclusion obligation without changing the filter, the declared Normal-Live P1 caps, the frozen P2 interface, or the canonical P3 gate.

## What is now mechanical

`tools/ou3_sea3_finite_window_response_admission.py` defines the exact finite-window admission predicate between the coupled physical SEA3 sea/RAO family and `Lhat_SEA3`.

A finite sea/ship response window can be admitted only when a separate validated realization producer supplies evidence that:

- uses validated arithmetic and outward enclosures;
- encloses the post-RAO response for every valid IMU sample in the window;
- is bound to the exact certified continuum RAO parameter box;
- does not use trajectory replay as proof evidence;
- does not substitute a Gaussian spectrum, RMS value, or PSD moment for a pathwise bound;
- keeps the post-RAO non-gravitational CoG acceleration norm at or below the existing `4 m/s^2` Normal-Live cap; and
- keeps body rate at or below the existing `30 deg/s` Normal-Live cap.

The response-box binding is recorded by a canonical SHA-256 digest so a finite-window producer cannot silently certify a different vessel-response family and feed the result into the SEA3 source bridge.

The predicate is deliberately fail-closed. A spectrum or response moment may be useful for constructing a validated realization enclosure, but it cannot by itself set the admission decision to PASS.

## What is still open

This does **not** promote

`L_actual_sea subset Lhat_SEA3`.

The missing producer is still the hard part: a replay-free oscillator/IQC, interval realization, or equivalent finite-window construction must prove that every physical sea/ship window claimed by the deployment theorem can furnish the required pathwise evidence. Until that producer exists, the global left inclusion remains false.

This also does not prune the 800-state P2 language. If the current full-P2 canonical H/A calculation passes the unchanged `1e-18` gate, the stronger full-P2 certificate transfers to SEA3. If it fails, the next numerical step remains response/estimator-history pruning followed by the same frozen H/A theorem interface.

## Regression found after #482 merge

The `source-foundation` job on the merged commit exposed one test-construction bug in `test_sea3_validator_rejects_mismatched_consumed_rao_box`: the in-memory inclusion object shares the response parameter-box object, so mutating the consumed box also mutated the source box before validation. The test now detaches the consumed candidate first, matching the serialized/stale-artifact mismatch that the validator is intended to reject. No validator condition is relaxed.

## Next certificate increment

1. Construct a deterministic finite-window realization producer against the new admission interface.
2. Prove coverage of the claimed physical SEA3 population, rather than merely checking individual supplied windows.
3. Promote the left inclusion only if that coverage proof closes.
4. In parallel, consume the clean current-head canonical P3 H/A result. If full P2 passes, proceed to source-complete P4 endpoint/prefix dissipation. If it does not, build SEA3-specific estimator-history pruning before rerunning H/A.
5. Attempt P5 finite capture only after canonical P3 and P4 are unambiguous.
