# OU-III ALT physical startup reachability

## Research question

Determine, for the current shipping implementation and the existing admitted corrected COMPLETE-BRMM / startup sensor / BIAS0/1/2 source class, whether every source history reaches Live in finite history-dependent time:

`for every admitted h, there exists finite T(h) such that startup reaches Live.`

A uniform `T_max` is useful only if it follows from the existing physical source class. The falsified 30,002-sample deadline must not be replaced by an enlarged guessed deadline.

## Mandatory physical relation

The reachability argument must retain one continuous physical history. In particular, sampled accelerometer directions are not independent controls. They are restrictions of the same corrected COMPLETE-BRMM trajectory and retain the coupled `p -> v -> a` relation, frequency support `[0.018, 0.88] Hz`, amplitude/velocity/acceleration and centered-primitive restrictions, physical attitude/rate history, sample clock, sensor residual temporal contract, and one BIAS history.

Any optimization used to characterize slow capture must parameterize or outer-approximate that continuous admissible history. An optimization over arbitrary per-sample vectors satisfying only `|a| <= 8.8` is inadmissible evidence for either capture or noncapture.

## Controlling alternatives

The proof must distinguish these mutually meaningful outcomes:

1. **Invariant separation.** Actual construction plus the physical source relation remains separated from the conditional antipodal/bad Mahony invariant set, and BRMM temporal cancellation forces finite capture.
2. **History-dependent capture without uniform deadline.** A progress/storage quantity has source-dependent finite excursion but every admitted history accumulates enough restoring apparent-gravity information to cross the shipping alignment predicate in finite time.
3. **Physical noncapture.** A single continuous corrected-COMPLETE-BRMM history, with actual sensor and startup machine ancestry, remains outside Live forever. Such a history falsifies universal eventual reachability and must be carried as a theorem counterexample rather than repaired by changing the filter or source class.

## Required analysis sequence

First reconstruct the exact 33,447-sample witness as a continuous physical trajectory and report its frequency, displacement/velocity/acceleration extrema, centered primitive, attitude/rate, residual histories, apparent-gravity directional bias over moving windows, Mahony attitude/integral trajectory, frontend alignment statistic, and the exact event that first permits Live. Compare every quantity to the existing source margins; do not call it representative merely because it is admitted.

Second derive a BRMM-to-sampled apparent-gravity temporal lemma. The target is a bound on accumulated or windowed directional forcing that uses the coupled continuous dynamics and oscillatory/primitive constraints. The lemma must be strong enough to enter the actual discrete Mahony/alignment recurrence; a continuous mean bound with no sampled-data implication is insufficient.

Third classify the conditional bad equilibria against the actual reset/first-sample seed and the reachable source product. Prove a separating invariant/barrier if one exists. If separation fails, construct the corresponding source history and test complete BRMM/sensor ancestry rather than injecting the equilibrium state.

Only after this classification should a universal finite-time argument be attempted. If no source-uniform deadline follows, state the theorem with finite `T(h)` and carry the source-dependent startup prefix into the later finite-word theorem.

## Anti-dead-end guards

- Do not enlarge the old timeout and call that a proof.
- Do not optimize independent sample accelerations, directions, gains, or Mahony states.
- Do not use replay maxima/minima as universal source bounds.
- Do not infer sampled cancellation from a continuous mean/primitive premise without a sampling lemma.
- Do not assume the antipodal equilibrium is reachable because it is an algebraic equilibrium, or unreachable because ordinary simulations avoid it.
- Do not change shipping constants or narrow corrected COMPLETE-BRMM to obtain capture.
- Keep `storage_search_allowed=false` until the existing finite-master guard closes.
- Preserve the independent P2/P3/P4/P5 route and frozen P3 threshold `1e-18`.

## Current status

The 30,002-sample universal deadline remains falsified. The known admitted witness reaching Live at sample 33,447 establishes neither a larger universal deadline nor eventual capture. Universal eventual startup reachability remains open until the physical temporal lemma and bad-set reachability classification are complete.
