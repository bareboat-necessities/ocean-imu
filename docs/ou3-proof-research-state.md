# OU-III proof research state

## Current hypothesis

One shipping-faithful path remains:
`construction -> finite capture -> finite H18 bridge -> finite reference
refinement/release -> recurring magnetically informed A21 -> regional practical
stability`. All stages inherit one physical execution. The estimator, physical
assumptions, and quality thresholds are unchanged.

The controlling linear inequality is now stated without an information-lifting
shortcut. Let M_k transport the root error through every preceding linear
prediction, applied correction, and covariance reset. For V=e'P^-1 e, exact
Kalman/Joseph algebra gives

`M_end' P_end^-1 M_end + D_word = P_root^-1`.

D_word includes every transported measurement loss H'S^-1H and every prediction
loss P^-1-F'P_pred^-1F. The required source-uniform inequality is
`D_word >= delta P_root^-1`, delta>0. It yields rho0<=1-delta. The nonlinear
small-gain, physical/arithmetic supply and every-prefix retention obligations
then enter the same finite-error theorem. A finite supplied word is algebraic
verification, not source-uniform evidence.

## Evidence

- `word_energy.py` proves the telescoping identity by exact rational matrix
  operations. Regression words include interleaved process noise, coupled
  measurements, a nonorthogonal reset, and singular PSD covariance inflation.
- The prior principal-block lifting implication is false. A positive-noise
  Kalman correction has full loss
  `D=[[1,0,1],[0,1,0],[1,0,1]]`, with restricted heading/bias loss I_2.
  Nevertheless e=(1,0,-1) has zero loss and preserves energy 20 exactly.
  This is an algebraic counterexample, not a shipping marine-history witness.
- The LIN endpoint path-action construction survives this distinction. For
  endpoint y=E_LIN x, a zero-root degree-seven Hermite path bounds
  `y'P_end^-1 y <= x'A x`. Therefore
  `P_end >= E_LIN A^-1 E_LIN'`, including carried cross covariance. This is a
  singular full-state lower bound; it cannot be added to unrelated marginal
  block floors without a joint comparison.
- `lin_matrix_certificate.py` exports the full rational 4x4 A and A^-1, with
  positive exact LDL pivots, rather than using the smallest pivot as an
  eigenvalue. It includes arbitrary piecewise-constant tau in [0.02,12],
  sample intervals [0.004,0.006], all regular interleaved accelerometer and
  integral corrections, the small-x source defects, and default PSD sync.
- The bound uses the existing default SpectralMSE configuration: sigma_aw>=.05,
  S_factor=1, nominal acceleration standard deviation>=.05, and actual S-noise
  standard deviation>=.15*.5=.075. It does not cover arbitrary user setters or
  alternate noise/cadence policies. Float32 implementation transfer is open.
- With fixed scales (2.4,18,132,4), the matrix comparison with raw neutral S
  observations gives mu_translation>=3.5959862602014e-11. Event times cover
  [0,.156], [8,8.156], [16,16.156] analytically. This is not a closed-loop loss
  bound or a full A21 rho. A 70-digit feasibility calculation preceded the
  rational enclosure; nominal event times gave approximately 4.41917e-11.

## Current limiter

The former mu_N>=2.04e-3 and matched-coordinate rho0 promotion relied on lifting
`E_hb' J E_hb >= I_2` to an independent full-state heading contribution. The
principal-block inequality does not control nuisance cross terms. Its use as
an unconditional full-state theorem premise is withdrawn. The geometric and
raw translation inequalities remain conditional algebra.

The carried, closed-loop M_k also includes gains/resets. It is not the bare
kinematic theta_0+C_k b_g used by the former long-word accumulation argument.
Positive information on a neutral restriction does not by itself supply the
full 21-state finite-word margin for the stable a_w/BA modes.

Recurring AG/BA factors through applied corrections and resets are not
certified by one-second process probes or a scalar bias iteration. In
particular, an LDL pivot is not a singular/eigenvalue floor. The additive
single-prediction gamma=1 lemma remains valid for factors of that actual Q;
it cannot substitute a corrected multi-step marginal factor for Q.

Consequently full covariance coercivity, full transported loss, explicit
nonlinear radius, float32 practical radius, finite capture, H18/release
retention and every-prefix tail retention remain open. The theorem is not
claimed. Physical qualification and recurring applied service remain required.

## Failed approaches / DEAD_ENDS

- Entrywise midpoint-radius Riccati refinement failed twice by dependency;
  prediction enclosure [-18.7907040,2502.41442] and failed acc/mag inverse
  boxes. Do not resume subdivision of that mechanism.
- Absolute cross-ceiling/Gershgorin reduction gave gamma=-5.1223e8. Preserve
  factors and joint action instead of dividing unrelated covariance ceilings
  by small scalar floors.
- Accumulated endpoint process noise plus per-event information ceilings does
  not survive arbitrary interleaved corrections; the exact 2-D counterexample
  is retained in the tests. Endpoint path action accounts for the sequence.
- Restricted magnetic information cannot be embedded as independent full-state
  information. Tighter scalar factors or longer words cannot fix that logical
  implication. Review full-loss/Schur or coupled path-action constructions.
- Inherited head 896d6912 failed Python compilation due to literal backslash-n
  separators; tests, YAML, JSON, and the scale-argument call were also broken.
  These serialization/API failures are repaired. The new matrix test initially
  passed string numerator/denominator values to Fraction; explicit integer
  parsing repaired that test harness failure. Neither changes mathematics.
- Local article rendering lacked IEEEtran.cls/luaotfload. The standard package
  installer failed with setgroups/seteuid permission errors. This is an
  environment dependency failure; the unchanged CI renderer supplies the
  required TeX distribution. No estimator or proof premise is affected.
  The workflow-contract command initially used the repository root and could
  not import its test modules; rerunning from tests/validation fixes that
  invocation error.
- Git transport can fetch but cannot push here (no HTTPS username credential).
  Publish the identical verified tree through the authorized GitHub connector
  with a non-forced update to the existing PR branch.

## Retained facts

The physical contracts, same-history bias mismatch, zero-DC displacement
condition, accepted-event service definition, Euclidean projection sector,
finite-bridge algebra, positive process-noise mechanisms, single-event Joseph
identity, and exact LIN path action remain. No shipping instability is inferred
from failure of a sufficient certificate.

## Alternatives

Within the same theorem, either bound the full D_word using block Schur
complements and the actual carried transport, or construct a coupled full-state
path-action/storage bound including AG/BA/LIN and the projection sector. A
principal-block restriction is insufficient in either construction.

## Next falsifiable experiment

Construct the full same-history A21 loss with carried covariance and actual
correction/reset transport, retaining all cross blocks and stable modes. Before
rigorous enclosure, compute a non-promoting high-precision feasibility ratio
for that construction. Show a positive normalized full-loss margin and a
quantified path to source-uniform enclosure; otherwise analyze its failing
mode. Do not report a new rho from the LIN matrix or the conditional raw S
normalization alone. Then certify the finite nonlinear region and additive
supply before attempting capture/release-to-tail composition.
