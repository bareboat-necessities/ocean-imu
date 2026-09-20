# OU-III stability handoff

Continue only in PR #560, branch `ou3-explicit-regional-practical-stability`.
Read AGENTS.md and `ou3-proof-research-state.md`. The same-history construction
-> capture -> H18 -> release -> A21 architecture is retained.

**The declared universal six-degree physical capture is refuted on the current
point-sample domain.** Read `ou3-sampled-capture-obstruction.md` before any new
positive construction. Its two nonzero continuous histories have identical
stationary measurements, meet all motion/bias bounds and all-time actual
magnetic service, and retain 53.13-degree physical tilt error through literal
startup and A21 release. Neither quiet-water exclusion nor tighter nominal
covariance estimates can remove this sampling ambiguity.

The exact rational artifact includes all-time envelope bounds, an invariant
source 2x2 Riccati comparison, and actual innovation service >1 for either
nominal or true down-axis convention. The native `sampled_capture-test` checks
the shipping float startup/refinement/release and stationary invariance;
it is not used to extrapolate all-time service or float32 totality.

A quantitative cell-integral sampling-fidelity condition is derived but is
not adopted or physically qualified. AGENTS.md does not authorize silently
changing physical assumptions. A positive end-to-end theorem needs a justified
revised domain or a different capture target, as well as the remaining loss,
nonlinear, retention and arithmetic work.

Retain the exact word energy and factor/range algebra, corrected LIN path
lower comparison, joint full-state covariance floor after 16 s, nuisance
covariance upper comparison after 17 s, and stationary A21 detectability.
Their scopes are conditional and do not imply physical capture. Do not revive
restricted-information lifting, endpoint-batch floors or entrywise Riccati
subdivision. The applicability map and finite quiet-word diagnostic remain
removed.

Reproduce with:

```
python3 -m unittest discover -s tests/validation -p 'test_ou3_*.py'
python3 tools/stability/ou3_theorem/build_evidence.py
make -C tests/kalman_ou_iii sampled_capture-test
./tests/kalman_ou_iii/sampled_capture-test
```
