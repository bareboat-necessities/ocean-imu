# OU-III controlling proof obligations — PR #560

The exact status is reproduced by `theorem_status.py` and committed in
`reports/results/ou3_stability/theorem-status.json`. Algebraic infrastructure
is distinct from a shipping source-uniform theorem certificate.

The amended MARINE MOTION contract adopts locally absolutely continuous
acceleration with jerk <=100 m/s^3. This excludes the old 200-Hz witness.
`ou3-sampling-fidelity.md` proves the resulting sharp sampled-mean bound,
32-s exclusion of fixed-attitude stationary-sample aliases, joint 3-D vector
information and physical LIN prediction supply in the full covariance metric.
The full corrected state loss and general capture remain open.

| Obligation | Current state | Required certificate |
|---|---|---|
| Joint recurring lower covariance | CLOSED in real arithmetic at regular A21 post-prediction roots after a 16-s window | `root_covariance_certificate.py`: convex combination of fresh AG/BA injection and corrected LIN matrix action, with all cross covariance retained |
| Full A21 information/loss and rho0 | OPEN for the conditional tail | Bound the complete transported loss `D_word >= delta P_root^-1`, delta>0, in all 21 coordinates; restricted magnetic service cannot be lifted to independent heading information |
| LIN/BA nuisance covariance upper bound | CLOSED for the regular default A21 profile after 17 s | `nuisance_upper_certificate.py` and `ou3-nuisance-upper-proof.md`: cancel the neutral root with three actual S observations; bound OU forcing, source Q defects and actual PSD sync; retain all nuisance cross covariance |
| Full covariance upper bound | OPEN: remaining AG6 block | Bound attitude/gyro covariance under varying realized coefficients and actual corrections/resets |
| Explicit nonlinear retained radius | OPEN | Bound the complete nonlinear remainder, including projection/reset/tuner behavior, against the verified strict linear margin |
| Whole-word float32 supply | Composition only | Literal operation counts, magnitude envelopes and certified prefix gains; real-arithmetic covariance positivity is not float32 totality |
| Finite startup/capture | OPEN under the jerk-bounded domain | Fixed-attitude stationary-looking aliases >=6 degrees are excluded; prove actual moving-attitude capture and retention through the shipping proxy, reference refinement and release |
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

The complete stability claim remains unproved under the revised domain. Consult `ou3-proof-research-state.md`
for the current failure classification and next falsifiable experiment.

The stationary A21 detectability proof in `ou3-stationary-detectability.md`
finds no nondecaying unobservable mode at rest with nonparallel gravity and
magnetic field. It is not a uniform varying-history theorem. Quiet water is
not excluded. The old sampling obstruction uses nonzero motion and is excluded by the
jerk condition, not by deleting quiet water.
