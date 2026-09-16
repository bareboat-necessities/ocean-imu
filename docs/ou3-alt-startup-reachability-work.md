# OU-III ALT physical startup reachability

## Research question

Determine, for the current shipping implementation and the existing admitted corrected COMPLETE-BRMM / startup sensor / BIAS0/1/2 source class, whether every source history reaches Live in finite history-dependent time:

`for every admitted h, there exists finite T(h) such that startup reaches Live.`

A uniform `T_max` is useful only if it follows from the existing physical source class. The falsified 30,002-sample deadline must not be replaced by an enlarged guessed deadline.

## Mandatory physical relation

The reachability argument must retain one continuous physical history. Sampled accelerometer directions are not independent controls: they are restrictions of the same corrected COMPLETE-BRMM trajectory and retain coupled `p -> v -> a`, frequency support `[0.018,0.88] Hz`, amplitude/velocity/acceleration and centered-primitive restrictions, physical attitude/rate, sample clock, sensor temporal contracts, and one BIAS history. An optimization over arbitrary per-sample vectors satisfying only `|a| <= 8.8` is inadmissible evidence.

## Controlling alternatives

1. Prove actual construction remains separated from the conditional antipodal/bad Mahony set and physical temporal cancellation forces finite capture.
2. Prove history-dependent finite capture without a uniform deadline.
3. Produce one all-time continuous corrected-COMPLETE-BRMM/sensor history from actual construction that remains non-Live forever; this falsifies universal eventual reachability.

## Required analysis sequence

Reconstruct and characterize the exact 33,447-sample witness, including physical margins and apparent-gravity/Mahony/frontend evolution. Derive a BRMM-to-sampled apparent-gravity temporal lemma strong enough to enter the actual discrete recurrence. Classify the conditional bad equilibria against actual reset/first-sample reachability. If a uniform deadline does not follow, use the theorem `forall h exists finite T(h)` and carry the source-dependent startup prefix forward.

The constructive circular docking candidate already present on main is the strongest critic of eventual capture. It must either be extended to an all-time invariant with source-uniform residual/rounding/guard bounds, or be shown to escape by a rigorous physical mechanism. A finite native non-Live prefix is not an infinite counterexample.

## Anti-dead-end guards

Do not enlarge the timeout, optimize independent sample accelerations, use replay extrema as universal bounds, infer sampled cancellation from continuous means without a sampling lemma, assume algebraic bad equilibria are reachable/unreachable without source ancestry, change shipping constants, or narrow BRMM. Keep `storage_search_allowed=false`. Preserve the independent P2/P3/P4/P5 route and P3=`1e-18`.

## Current status

The 30,002-sample universal deadline is falsified. The 33,447-sample recovery proves neither a larger universal deadline nor eventual capture. Universal eventual startup reachability remains open pending the physical temporal lemma and all-time classification of the docking/bad set.
