# OU-III stability handoff

This continuation starts from main `257dce56d195276bd6d237c0f8cd653c016eec24`
and stays in one branch/PR. Continue from its current head until merged;
afterward start from latest main. Do not restart an older proof branch.
The theorem is **not closed**; merging this work does not certify stability.
Keep the single same-history construction -> capture -> H18 -> release -> A21
architecture and the unchanged shipping estimator.

Read in order, with paths relative to the repository root:

1. `AGENTS.md`.
2. `docs/ou3-proof-research-state.md` and
   `docs/ou3-controlling-proof-obligations.md`.
3. `reports/results/ou3_stability/theorem-status.json` for the machine-checked
   distinction between established implications and undischarged premises.
4. `docs/ou3-corrected-word-proof.md`, `docs/ou3-sampling-fidelity.md`, and
   `docs/ou3-nuisance-upper-proof.md`, alongside their modules under
   `tools/stability/ou3_theorem/` and the shipping source operations they cite.
5. `docs/ou3-ag-readout-proof.md` and `ag_readout.py`: the historical
   six-column action and the exact forward-prior obstruction.

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
polynomial defect. The exact covariance recursion also bounds all prediction
and measurement inputs jointly by the square root of their summed actions;
reset/projection defects retain separate norm supplies. This avoids an
unnecessary long-word factor from summing every individual input norm.
The reset remainder has an injection-squared times error
term; do not call it purely quadratic in error at a nonzero injection. No
second filter with identical gains or event decisions is assumed.

Next discharge the actual six-column J premise, then combine these supplies,
projection sector and physical mismatch with capture/release and every-prefix
retention. General capture, full AG upper covariance, strict uniform loss and
arithmetic remain open. No completion percentage or full theorem is justified.

The latest addition proves a conditional historical AG covariance comparison.
Backward factor action cancels an arbitrary AG root exactly, preserving all
nuisance/process correlations and applied corrections/resets. A uniform
action ceiling would supply J through a matrix first-prediction comparison.
The implementation rejects even an arbitrarily small uncancelled AG root;
one fixed numerical reader cannot be reused over varying coefficients.
Uniform source rank/action and inherited entry linkage remain open. The
supplied exact/80-digit examples are not shipping reachability evidence.

The exact failed relaxation is `D_AG,AG >= 10^-6 I6` for unbounded priors:
at `diag(10^12 I6,I15)` the Rayleigh margin is at most `-9.99999e-7`.
No reachability or magnetic-service membership is claimed for that prior.
Next test coefficient-dependent historical readers against actual varying
source windows and enclose their 6x6 action; do not refine the invalid
forward-only unrestricted-prior tactic or promote a finite replay.

Specifically, at pre-prediction roots, establish a single positive-definite
matrix J with `D_word[AG,AG] >= J` uniformly over admitted histories, where AG
contains all three attitude-error and three gyro-bias coordinates. This is a
principal block of the **actual complete corrected-word loss**, including the
realized gain and reset transports. A physical Gramian, stationary-system
detectability, or a positive diagnostic on one supplied word does not discharge
this premise. Preserve loss factors and nuisance coupling when seeking the
bound. Then apply the existing Schur/first-prediction implication, quantify its
margin, and charge the finite-error supplies in the same storage before claiming
nonlinear capture, release, or retention. Capture may have a history-dependent
finite time; no such time has yet been certified.

The main article is `doc/kalman_ou_iii/kalman_ou-w3d.tex`; the stability study is
a separate document. PR #560 moved the dataset description into evaluation
methodology, added engine-noise prefiltering to the deployed-filter diagram,
made the four requested lever-arm figures single-column, completed bias and
comparative-results prose, and removed the identified implementation names.
The revised main article compiled to 32 pages and its affected pages were
visually checked. Two existing ensemble-table overflows and one existing
full-state-equation overflow remain. These editorial changes are not proof
milestones.

Reproduce with:

```
python3 -m unittest discover -s tests/validation -p 'test_ou3_*.py'
python3 tools/stability/ou3_theorem/build_evidence.py
```

Use a Python environment with NumPy and SciPy; diagnostics may also need
mpmath. The broader publication/evidence suite uses Matplotlib and pandas:

```
make -C tests/validation evidence-test
```

At the article handover, that suite passed 463 tests with one existing
data-dependent skip. Recheck CI at the final PR/merge commit rather than
assuming an earlier run covers later edits. Evidence validation must continue
to report `theorem_closed: false` until the remaining obligations are proved.
