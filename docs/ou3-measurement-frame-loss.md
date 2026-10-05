# Moving-frame magnetic loss and the remaining word obligation

This continues the [governing strategy](ou3-information-shear-strategy.md).
The proof is in the authoritative appendix, `app:measurement-frame-loss`;
`measurement_frame.py` reproduces the rational inequalities and identities.

## Lemma and finite-error role

**PROVED — analytical, qualified implication.** For every regular real-operation
planar word on the specified wrapper profile, with actual SPD covariance,
qualified reference cone and applied magnetic service, the pitch/BG covariance
reader below is uniform over the inherited root covariance. At each covered
magnetic substep, additionally suppose the retained nominal/central-physical
pitch chart is at most six degrees and the inherited comparison energy
`e'P^-1 e <= 100`. The complete linked mean/covariance *magnetic substep* then
satisfies

\[
 W^+-W\le-\frac9{10}(L_\eta+L_P),\qquad
 W=\eta^TP^{-1}\eta+\operatorname{tr}(P^{-1}dP P^{-1}dP).
\]

These are activation-domain implications, not new physical assumptions and
not forward-invariance assertions. The cap 100 is a sufficient comparison
energy for this inequality, not a claimed reached bound. Reference acquisition,
chart/energy retention and full binary32 transfer remain OPEN. The lemma does
not assume or prove a 20-second periodic return. Existing service is assumed
for this result; it does not admit the special planar witness.

The literal structure used is the same P/K/Joseph correction in the moving
nominal-world frame `T=diag(R,R,I,I,I,I,R)`. The transformed rows are
`Hacc=[-[aw-g]x,0,0,0,0,I,I]` and `Hmag=[-[Bref]x,0,...,0]`.
Attitude-dependent row derivatives cancel, while frame connection, AW,
reference/noise changes and injection/reset transport remain. This is a
coordinate change of the shipping estimator, not a replacement observer.

The magnetic covariance-induced residual term is absorbed with its square
charge into the loss of the same correction. No independent K/r/P extrema,
finite replay factor, fixed-root covariance ball or per-AW norm product is used.

## Actual calculation

For the planar fixed-reference/noise substep, let `h` be its scalar even-sector
pitch row, `R_m=RN32(.8)^2`, `s=hPh'+R_m`, and

\[
 t=hPh^T/R_m,\quad \kappa=B^Tr/|B|^2,\quad
 x=h\eta,\quad z=h\,dP P^{-1}e,\quad E=e^TP^{-1}e.
\]

The frame connection reduces **exactly** to
`eta+ = A eta - kappa K h (eta+dP P^-1 e)`.
Only the radial innovation enters this endogenous magnetic term. The exact
rank-one Fisher loss gives `z^2/s <= E*(1+t)/(2+t)*L_P`. Therefore

\[
 W^+-W=-aL_\eta-L_P
 -2\kappa(1-\kappa t)xz/s+\kappa^2t z^2/s,
 \qquad a=1+2\kappa-\kappa^2t.
\]

The last positive square is retained. On `|kappa|<=3/500`, `t<=6`, `E<=100`,
the two-by-two loss matrix minus `(9/10)I` has positive diagonal entries and
determinant at least `4673/1250000`. This is an exact rational certificate.
It proves the displayed loss inequality for all operands in that stated domain.

The required pitch covariance cap is now **derived**, not assumed. In the
planar pitch/BG marginal the literal transition is `[[1,h],[0,1]]`, reset is
identity, and AW leaves the marginal unchanged. Other optimal corrections
decrease it. Reached held-BA corrections have the proved active optimal
reduction. A two-observation unbiased reader eliminates the inherited root.
With IMU steps in `[.004,.006]`, choose prior service windows at the last IMU
roots no later than `t-3-hmax` and `t-1-hmax`. Both windows finish strictly
before the correction at t. Their event separation is in `[497/500,1503/500]`
and last-event age is at most `253/250`. After `753/250` seconds of qualified
service, exact reader variances give

\[
 P_{\theta_y\theta_y}<3/5000,\qquad
 P_{b_{g,y}b_{g,y}}<1/4000.
\]

This is a covariance-order reader bound on the actual marginal; the shadow
reader is an upper bound, not substituted estimator dynamics. No physical
noise process is required to follow the estimator's process prior. The wrapper
passes `b0=1e-10f` through `initialize_ext`; the standalone core default differs.
The proof is a real-operation lift with represented configuration constants;
the optional PSD-repair outer charge is not a floating-point error budget.

The qualified reference cone and six-degree chart imply `|kappa|<3/500`.
The covariance cap and reference norm imply `t<6`. They further give the
uniform **pitch-row** bounds

\[
 L_\eta\ge1250\eta_{\theta_y}^2,\qquad
 L_P\ge1400(dP P^{-1}dP)_{\theta_y\theta_y}.
\]

Thus the coupled correction loses at least
`1125*eta_pitch^2+1260*(dP P^-1 dP)_pitch,pitch`. This is genuine positive
row coercivity, not a positive full transverse word gap. Odd covariance Fisher
loss remains present; the nominal tangent here is restricted to the proved
planar invariant stratum, not every 21-state perturbation of its orbit.

## Linked word, kernel and minimal remaining dependencies

At S corrections the world S row is constant, `H_S Omega=0`, and default
one-way tuner dependence gives `dR_S=0` for the MEKF-root homogeneous variation.
Hence their pre-injection coupled loss is exact without an NIS cap. Source
variation of R_S and held-bias ports remain in the forced word.

Split literal operations at additive correction/injection/frame transport.
Use the exact linked word identity with the same moving-frame storage. On the
covered homogeneous planar domain it now gives

\[
 W_N-W_0\le-\widetilde{\mathcal A}_W+\mathcal F_{\rm remaining},
\quad
 \widetilde{\mathcal A}_W=
 \tfrac9{10}\sum_{\rm mag}(L_\eta+L_P)
 +\sum_S(L_\eta+L_P)+\sum_{\rm other}\ell_i.
\]

`F_remaining` is the exact signed cross-plus-square work of accelerometer
AW-row feedback, prediction/source/dF/dQ, reset and frame changes, projections,
AW faces/targets, reference/noise variations and moving quotient conversion.
No endogenous term is declared an independent disturbance. Gauge/source
forcing and fibre curvature are added through their actual linked ports.

The magnetic coupled loss has zero action exactly when `h eta=0` and `h dP=0`
in the active even sector (plus the ordinary odd covariance loss kernel).
The radial feedback also vanishes there. The positive weights leave the base
word action kernel unchanged. They do **not** identify it with physical
compatibility: same-record physical variation leaves nominal R/P unchanged,
and the physical BA variation carries OU mismatch. The physical roll/BA gauge
is not a planar nominal-pitch perturbation. No gauge direction is deleted here.

For a justified transverse root lift B, let `G_Atilde` be the actual-prefix
Gram of this improved action, `G_remaining` the signed work matrix, and
`Gamma_i` the endpoint gauge-energy matrices. The precise remaining gap is

\[
 B^T(G_{A\!\sim}-G_{\rm remaining}-\Gamma_0+\Gamma_N)B
 \succeq c B^TJ_{\perp,0}B,\qquad c>0.
\]

| Activation fact | Role | Status |
|---|---|---|
| Planar parity, CoG profile, reached held-BA manifold | Exact rows, pitch marginal and masks | PROVED analytical on stated branches |
| Applied service and qualified reference cone | Two prior observations and their effective noise | Service is an existing premise; planar admission and float reference acquisition OPEN |
| Pitch/BG marginal covariance ceiling | `t<6`; row coercivity | PROVED analytical implication after the prior-window horizon |
| Six-degree chart and comparison energy `E<=100` | Radial fraction and absorption | Uniform implication proved; inherited every-prefix activation OPEN |
| Same-history reference/noise and event branches | Decide which ports vanish | Fixed during covered magnetic/S substep; variations retained elsewhere |
| Remaining signed work and complementary transverse action | Uniform complete-word `c>0` | OPEN; the pitch-row estimate alone does not control all states |
| Gauge/source, chart conversion, nonlinear/arithmetic charges | Finite-error quotient/tube, every-prefix retention | OPEN beyond retained local formulas |

No absolute p/v box is introduced. No radius or all-future magnetic floor is
solved before the complete-word gap and forcing charges close.

**Structures preserved:** all 21 states/P/K/Joseph, held masks, complete
frame derivative, physical bias mismatch, reference/gates and inherited clocks.
**Relaxations introduced:** real-operation regular-branch scope; qualified
planar domain, scalar unbiased-reader upper bound, rational Schur sufficient
bound. These enlarge only the stated inequalities, not the physical premises.
**Failure record:** E for the corrected pre-correction window endpoint
bookkeeping; the current implementation uses strictly prior windows. D remains
for the unclosed complete-word sufficient gap. No admitted counterexample.
**Next calculation:** absorb the accelerometer AW-row port and prediction/reset
work against the linked process/S and retained magnetic loss, keeping physical
gauge/source forcing. Prove precisely the activation facts that estimate uses.
