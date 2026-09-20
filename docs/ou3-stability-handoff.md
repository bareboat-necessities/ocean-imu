# OU-III stability handoff

Continue only in PR #560, branch `ou3-explicit-regional-practical-stability`.
Read AGENTS.md and `ou3-proof-research-state.md`. Keep the single same-history
construction -> capture -> H18 -> release -> A21 architecture.

MARINE MOTION now includes locally absolutely continuous physical acceleration
with jerk <=100 m/s^3. The user authorized this domain revision. Quiet water,
all other physical limits and shipping behavior are retained.

Read `ou3-sampling-fidelity.md`: it proves the sharp nonuniform trapezoidal
mean bound, exclusion of every constant-attitude stationary-sample alias of
at least six degrees over 32 s, a positive joint 3-D measured-vector Gramian
after 64 s (both magnetic residual envelopes charged), finite-rotation
coercivity, and physical LIN prediction supply in the full covariance metric.
The exact constants are reproduced by `sampling_fidelity.py`. The physical
Gramian uses true world transport; it is not the actual corrected full loss.

The 200-Hz counterexample and native `sampled_capture-test` now have the
explicit scope of the previous domain without a jerk bound. They explain why
the new premise is necessary and remain useful regression checks. They do
not refute capture under the revised domain.

Retain the exact complete-word loss/factor/range algebra, corrected LIN path
lower comparison, joint full-state covariance floor after 16 s, nuisance
covariance upper comparison after 17 s, and stationary A21 detectability.
Next prove the full transported 21-state loss with nuisance cancellation,
then actual capture/release into a retained nonlinear region and arithmetic
transfer. General capture, AG covariance upper bounds and strict dissipativity
remain open; no proof-completion percentage or full theorem claim is justified.

Reproduce with:

```
python3 -m unittest discover -s tests/validation -p 'test_ou3_*.py'
python3 tools/stability/ou3_theorem/build_evidence.py
```
