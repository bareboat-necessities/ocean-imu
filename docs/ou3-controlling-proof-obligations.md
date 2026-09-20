# OU-III controlling proof obligations — PR #560

The exact status is reproduced by `theorem_status.py` and committed in
`reports/results/ou3_stability/theorem-status.json`. Algebraic infrastructure
is distinct from a shipping source-uniform theorem certificate.

| Obligation | Current state | Required certificate |
|---|---|---|
| Joint recurring lower covariance | CLOSED in real arithmetic at regular A21 post-prediction roots after a 16-s window | `root_covariance_certificate.py`: convex combination of fresh AG/BA injection and corrected LIN matrix action, with all cross covariance retained |
| Full A21 information/loss and rho0 | OPEN; controlling gap | Bound the complete transported loss `D_word >= delta P_root^-1`, delta>0, in all 21 coordinates; restricted magnetic service cannot be lifted to independent heading information |
| Uniform covariance upper bound | OPEN | Joint source-uniform detectability/path-action comparison including literal corrections, resets and PSD sync |
| Explicit nonlinear retained radius | OPEN | Bound the complete nonlinear remainder, including projection/reset/tuner behavior, against the verified strict linear margin |
| Whole-word float32 supply | Composition only | Literal operation counts, magnitude envelopes and certified prefix gains; real-arithmetic covariance positivity is not float32 totality |
| Finite startup/capture | OPEN | History-dependent finite capture from construction into the retained domain; no common deadline is assumed |
| Finite H18 retention | Composition only | Actual entry set, history-dependent bridge duration and supply/gain bounds |
| Reference refinement and bias release | Conditional captured-domain completion only | Retain the actual refinement state machine, accepted-update count, one-second guard and covariance release |
| Release into the A21 region | OPEN | Compare the certified release set with the nonlinear retained region/projection-sector bound |
| Every-prefix tail retention | Composition only | Uniform bounds at every intermediate operation, including between recurring prediction roots |
| Recurring magnetic service | Continuation schema only | One-history every-window certificate from actually applied informative corrections and actual innovation covariance |
| Physical/sensor/bias qualification | Composition only | Simultaneous MARINE MOTION, IMU BIAS and MAGNETIC SERVICE continuations plus assembled sensor/mount/calibration qualification |
| Implementation/arithmetic totality | OPEN | Every finite branch, innovation solve, projection/reset, scheduler and arithmetic enclosure |

Do not revive the retired entrywise Riccati subdivision, endpoint-batch floor,
restricted-information lifting, or scalar normalization of an unproved full
Gramian. Local interval kernels remain available as conditional algebra; no
recurring interval box is currently certified. The active construction retains
matrix factors and the full covariance-energy loss identity.

The complete stability theorem remains open. Consult `ou3-proof-research-state.md`
for the current failure classification and next falsifiable experiment.

The supporting [literature map](ou3-literature-applicability.md) identifies
exact imported-result hypotheses and their remaining shipping gaps. The
[factor construction](ou3-factor-loss.md) preserves three-row losses,
process blocks and singular nuisance ranges. Its carried quiet-water full21
rho≈.833 is diagnostic only; the full-state loss row above remains OPEN.
