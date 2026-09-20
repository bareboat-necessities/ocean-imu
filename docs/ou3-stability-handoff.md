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
Read `ou3-corrected-word-proof.md` next. The nuisance floor is now transported
to roots immediately before prediction, without an injection or nominal AW
bound. A Schur-complement argument proves that a uniform SPD bound J on the
six AG columns of the actual complete loss implies a full covariance upper
bound; the first prediction then supplies strict 21-state contraction. J is
still open. This is not the invalid lifting of two-column magnetic service.

The same document proves actual-gain finite-error word composition and an
explicit finite-angle reset bound, including the real small-angle quaternion
polynomial defect. The reset remainder has an injection-squared times error
term; do not call it purely quadratic in error at a nonzero injection. No
second filter with identical gains or event decisions is assumed.

Next discharge the actual six-column J premise, then combine these supplies,
projection sector and physical mismatch with capture/release and every-prefix
retention. General capture, full AG upper covariance, strict uniform loss and
arithmetic remain open. No completion percentage or full theorem is justified.

Reproduce with:

```
python3 -m unittest discover -s tests/validation -p 'test_ou3_*.py'
python3 tools/stability/ou3_theorem/build_evidence.py
```
