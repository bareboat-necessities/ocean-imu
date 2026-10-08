# OU-III stability handoff

Continue from the latest head of the requested active PR branch and retain
one proof path. A merged historical proof commit is not a stability certificate.
The theorem is **not closed**; merging this work does not certify stability.
Keep the single same-history construction -> capture -> H18 -> release -> A21
architecture and the shipping estimator, including its residual gyro-bias protection.

Read in order, with paths relative to the repository root:

1. `AGENTS.md`.
2. `docs/ou3-proof-research-state.md` and
   `docs/ou3-controlling-proof-obligations.md`.
3. `reports/results/ou3_stability/theorem-status.json` for the machine-checked
   distinction between established implications and undischarged premises.
4. `docs/ou3-corrected-word-proof.md`, `docs/ou3-sampling-fidelity.md`, and
   `docs/ou3-nuisance-upper-proof.md`, alongside their modules under
   `tools/stability/ou3_theorem/` and the shipping source operations they cite.
   Its section 7 (`word_diameter.py`) reduces the word contraction to the
   kernel-bounded Riccati diameter; that is the current controlling target.
5. `docs/ou3-ag-readout-proof.md` and `ag_readout.py`: the historical
   six-column action and the exact forward-prior obstruction.
6. `docs/ou3-regime-design.md`, `docs/ou3-moving-six-pivots.md`, and the
   retained `docs/ou3-signed-temporal-continuation.md`.
7. `docs/ou3-world-frame-rows.md`: attitude-free historical rows, the literal
   injection budget, the same-cell collinear witness and the AW transfer.

MARINE MOTION retains locally absolutely continuous physical acceleration
with jerk <=100 m/s^3. Excitation is required only on every complete T_E window contained in one
maximal physical moving episode between nondegenerate rest intervals. Isolated
zero-rate instants do not restart an episode. T_E/theta_E remain symbolic.
Windows crossing rest boundaries carry finite-transition obligations; no
arbitrarily short departure window owes a full positive span. Read
`ou3-regime-design.md` for the exact quantifiers, stationary observability,
indistinguishable rest/motion histories and detector requirements. Quiet
packets cannot certify physical STILL with both finite entry and guaranteed
finite exit under the existing bias bounds. The shipping estimator is unchanged;
stationary practical robustness and certified transition retention remain OPEN.

The precise zero-residual quiet nominal subcase now has a root-independent
historical action, every-operation full covariance ceiling and qualitative
homogeneous linear loss; see `ou3-stationary-detectability.md`. These are not
physical compatible-class or nonlinear results. The moving reader uses
C=A^-1 B to retain resets without erasing accumulated gyro action, and
actual acc/mag groups within one prediction cell factor exactly with E=0.
Uniform geometry, injection/action bounds and nonlinear retention remain OPEN.

The moving six-column premise is now an aggregate world-frame statement.
Every historical row is `-R_k[f_k]x[A~_k R_0',B~_k]`, so the attitude estimate
enters only through world injections and the nominal rotation integral. A
1-Hz collinear motion/bias history degenerates same-cell groups, although its
cadence fails MAGNETIC SERVICE; a jerk lemma ties same-cell geometry to the
applied magnetic cadence. Corollary A transfers the physical transverse force
to the attitude columns at the cost of the AW tracking error
(<1.12383 m/s^2 on 16-s windows). The AW covariance ceiling `(1+eps)16` is
tight at every sync and covariance-normalized storage fails on the carried
collinear motion, so bound that error physically from the literal
correction loop, then the aggregate gyro columns; charge injections by their
signed world sum.

Read `ou3-sampling-fidelity.md`: it proves the sharp nonuniform trapezoidal
mean bound, exclusion of every constant-attitude stationary-sample alias of
at least six degrees over 32 s, a positive joint 3-D measured-vector Gramian
after 64 s (both magnetic residual envelopes charged), finite-rotation
coercivity, and physical LIN prediction supply in the full covariance metric.
The exact constants are reproduced by `sampling_fidelity.py`. The physical
Gramian uses true world transport; it is not the actual corrected full loss.

The 200-Hz counterexample and native `sampled_capture-test` now have the
explicit scope of the previous domain without jerk or tilt-span requirements.
They remain useful historical regression checks. They do
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

Next bound the kernel-bounded word diameter of `ou3-corrected-word-proof.md`
section 7 source-uniformly, then combine these supplies, projection sector
and physical mismatch with capture/release and every-prefix retention.
General capture, full AG upper covariance, strict uniform loss and
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
Carried quiet and moving windows now have coefficient-dependent factor readers;
both exported coefficient sequences have exact full-matrix action enclosures.
These checks are not source-uniform. Independent nominal coefficient ranges
fail the all-row Schur test at rank four, with exact Rayleigh margin -1
against I6. See the source-shaped annihilator in `ou3-ag-readout-proof.md`.
The missing argument must use the carried nominal mean/innovation/reset
recursion to exclude sustained near-null behavior; do not refine the invalid
forward-only unrestricted-prior tactic or promote a finite replay.

The inherited BA covariance bound independently gives a projection guard:
pre-projection sqrt(V)<=6 makes the literal projection inactive, with
.02483339501604595 m/s^2 interior slack. Every-prefix retention and entry
into that guard remain open. The finite-angle reset defect is still present.

Reproduce the optional source experiment with
`python3 -m tools.stability.ou3_theorem.ag_readout_source_diagnostic --output /tmp/ou3-source-action.json`
(`--eigen` selects an alternate Eigen path). It compiles an observer and an
untapped control, checks exact terminal parity, runs 80-digit diagnostics and
a rational exported-word action check. This is not float32 totality.

Specifically, bound `kappa_nu=lambda_max(Pi^-1 P_nu)` uniformly over admitted
words. Here `Pi` is the known-root terminal covariance and `P_nu` the terminal
covariance started with information only along the physical tilt/BA kernel
`nu`. Bound it with the realized gains and resets, and keep the scalar kernel
ceiling invariant. The first-prediction route is capped (DEAD_END 25) and a
full root covariance ceiling is not needed for `rho_0` (DEAD_END 27). A physical
Gramian, stationary-system detectability, or a positive diagnostic on one
supplied word does not discharge this premise. Then quantify the margin, and
charge the finite-error supplies in the same storage before claiming
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

The complete moving sync/symmetry operation now has an exact signed-factor
upper envelope, including its -2^-44 Rayleigh defect and all 21 coordinates.
Quiet and moving exported words both pass exact action enclosure. Prediction,
solve, Joseph, reset and state arithmetic are not thereby certified.
The nominal-history attack also tests a zero-innovation full-turn gyro alias:
h=.005, b_hat_g=-400 pi e_z, quiet inputs. It satisfies the regular nominal
mean recursion in real arithmetic, but its AG array has rank four. Neither
construction reachability nor magnetic service is certified. Next derive a
construction-linked estimate gyro bound; a physical gyro bound or innovation
bound alone does not supply it. All uniform and end-to-end claims stay false.

The latest construction attempt is reproducible with
`python -m tools.stability.ou3_theorem.construction_history_diagnostic --output /tmp/ou3-construction.json`.
Its R=I translation fails the current tilt-span assumption; it is retained
only as an old-domain diagnostic. It obeys the other motion/bias envelopes and
passes through the unchanged wrapper from begin to 600 s. The 400--600 s
native tilt is about 8.21 degrees, so a six-degree entry set is not inherited
merely from refinement/release. This finite observation is not an all-time
magnetic-service certificate or a refutation of eventual capture. The nominal
acceleration reaches 9.776391 despite physical norm <=sqrt(72)<8.8. Use the
actual construction for the joint mean/covariance enclosure; no assumed
nominal-state box or free-root family can replace it.

The recorded terminal full covariance is also verified SPD exactly. The known
zero true BA and its correlated 3x3 marginal give a full-storage lower bound
>11492.6752, hence the candidate V<=36 guard has margin <-11456.6752 at that
endpoint. `build_evidence.py` reproduces this exact certificate and the native
driver fingerprint. This finite endpoint exclusion is not eventual-capture
refutation and does not establish all-time magnetic service.

The construction mean-action attempt is in `ou3-construction-mean-action.md`
and `construction_mean_action.py`. It carries 89998 predictions and 122251
accepted rank-three corrections from zero means, retains the full BG/AW
matrix, and encloses actual mean rounding. Rational interval factors prove
||b_hat_g||<1 at every recorded prefix; exact checks prove force/field sine
>2/5 at each of the 40000 tail accelerometer observations. These finite
closures are not a uniform exclusion or readout action ceiling.

The energy-only enclosure itself fails `E_col-E>0` with exact margin in
[-817884.035305,-817884.035304], even though the actual finite tail remains
separated. No new physical counterexample follows. Retain the signed
innovation relations and physical/nominal integral dynamics in the next
whole-window argument; do not tighten this failed pointwise ellipsoid or
substitute a nominal-state cap. The committed mean-action summary is checked
with rational interval LDL, and CI reproduces the full binary transcript,
80-digit diagnostic and outward enclosure. Uniform B_*, J and rho0 remain
false. The newer MARINE MOTION tilt-span
restriction excludes this particular moving diagnostic.

Read `ou3-signed-temporal-continuation.md` before continuing. The exact carried
adjoint criterion and signed physical bias/velocity identities are established.
Zero-mean weights do not preserve adjoint compatibility. Positive temporal
margins, their quantitative six-pivot implication and coefficient compactness
remain open. The implemented .5 rad/s gyro ball excludes complete-turn bias
aliases and gives a one-step transport floor on qualified 4--6 ms predictions;
it does not close the signed temporal gyro margin. Preserve both
compatibility residual sums and the quiet-bias stillness ambiguity.

The chronological two-group lemma gives a sufficient actual-row defect budget
for all six historical pivots. Its uniform reset/rotation/reference/asynchrony
premises remain open. The excited moving diagnostic already fails the recorded
signed norm relaxations; deleting indefinite stillness cannot make those bounds
positive. Prove the coupled budget, then finite B_*, J_AG and strict full loss.
Do not enable a quiet-packet bias hold or zero-motion observation before a
stationary theorem covers the entire measurement-compatible class. Repeated
finite bridges also require a proved cumulative storage budget.


## PR #643 handoff — current controlling state (2026-10-02)

PR #643 on branch proof/ou3-linked-local-defects is the active continuation. The single path remains construction -> capture -> H18 -> release -> A21 -> regional practical stability. Do not revive historical architectures marked dead in the research ledger.

The decisive executable path is the source-uniform physical HistoryCell release enclosure. The literal leaf carries AccelVibrationGuard, Mahony, adaptive band/variance/frequency/tuner, same-history Racc, interval mean/covariance/gain, K r, reset and S scheduler. Coupled raw-force and quaternion norm invariants prevent Cartesian dependency artifacts. FAST guard forcing is charged through a complete LP-I signed-primitive Abel state instead of 2 B_a,f at every sample. Integrated guard-reader FAST forcing is paired samplewise; the terminal AW reader is generated by backward interval adjoint through literal events.

Enclosure failures carry actual time. Subdivision may split only physical/sensor root histories (p, v, slow accel/gyro, FAST primitives, gravity direction, magnetic field), never generated w, Mahony state, tau, sigma_aw, R_S, T_S, covariance, gain or scheduler phase. A converged admitted physical leaf that literally violates a required condition is a counterexample, not indefinite class D.

The stability-study manuscript is doc/kalman_ou_iii/kalman_ou-w3d-stability-study.tex. Detailed current proofs are in doc/kalman_ou_iii/kalman_ou-w3d-stability-proofs.tex-part. Historical witnesses, dead architectures and retracted lemmas were removed. The live lemma chain is labeled PROVED/CONDITIONAL/OPEN; conditional implications are proved with their unresolved premises named.

Machine mapping is reports/results/ou3_stability/lemma-manifest.json. It maps every current LaTeX lemma/theorem to stable ID, status, assumptions, proof modules, certificates and open dependencies. tests/validation/test_ou3_lemma_manifest.py enforces exact manuscript coverage and referenced paths. Keep it synchronized on any lemma change or promotion.

The regional theorem remains OPEN. Do not promote continuous MAGNETIC SERVICE transport-stratum coverage; uniform A* physical-to-nominal margin over all HistoryCells/prefixes; finite source-uniform H18 release completion; complete physical HistoryCell release/A21 enclosure; recurring finite-error A21 contraction/prefix retention; linked sup_h inf_gamma chi_gamma/gamma < .0225; nonlinear remainder; float32 transfer; or regime composition until their gates pass.

A separate existing theorem-workflow failure remains in ag_readout_source_diagnostic.py: event 4 lacks physical lift fields. It is not the causal HistoryCell enclosure obstruction.

Before continuing, read docs/ou3-proof-research-state.md, docs/ou3-controlling-proof-obligations.md, this handoff, both stability-study TeX files, lemma-manifest.json, theorem-status.json and literal shipping OU-III source. Check CI at the actual continuation head.
