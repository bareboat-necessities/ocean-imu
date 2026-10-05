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
**Word continuation:** [linked suffix scores](ou3-linked-word-score.md) regroup
this same exact word work, retaining covariance variations created inside it.
They do not add a second magnetic charge.

**Next calculation:** absorb the accelerometer AW-row port and prediction/reset
work against the linked process/S and retained magnetic loss, keeping physical
gauge/source forcing. Prove precisely the activation facts that estimate uses.

## Exact AW shear removes the remaining nominal accelerometer row derivative

**PROVED — analytical**, on each regular CoG real-operation branch. The
continuation in `app:aw-shear-loss` uses the invertible full-state frame

\[
 N(a)=-E_{aw}[a]_\times E_\theta^T,\quad L(a)=I+N(a),\quad
 L(a)^{-1}=I-N(a).
\]

`N(a)N(b)=0` for every a,b. At the actual nominal AW state,
`Hacc L^-1=[+[g]x,0,0,0,0,I,I]`: the row is independent of AW. Magnetic and
S rows are unchanged. Held-BA masks commute with this transformation, so the
literal innovation/gain/Joseph operands are preserved, including arbitrary
held cross blocks. Optimal loss formulas still require the reached active
reduction and retain `dR_eff=dRacc+dPBA`.

This is a change of proof coordinates, not an invariant-EKF replacement.
No bound on nominal AW is needed for invertibility. For the comparison
orientation `e_aw=aw_hat-a_physical`, its AW component is
`e_tilde_aw=(I+[e_theta]x)e_aw-[a_physical]x e_theta`.
The inverse of `I+[e_theta]x` has norm at most one. Thus no tangent direction
is deleted. Physical/source variation retains the derivative of physical AW.
Point physical energy is unchanged: `e_tilde' P_tilde^-1 e_tilde=e'P^-1e`.
The same-record physical fibre has unchanged nominal L and transforms with
its precision Gram intact; a constant gravity row is not bias identification.

With `Gamma=dL L^-1=N(daw)`, the transformed differential is

\[
 d\widetilde P=LdPL^T+\Gamma\widetilde P+\widetilde P\Gamma^T,
 \qquad
 \widetilde\eta=L\eta-\widetilde P\Gamma^T\widetilde J\widetilde e.
\]

The actual accelerometer covariance row port is exactly the boundary difference

\[
 -K\,dH\,C-C\,dH^TK^T
 =A(\Gamma P+P\Gamma^T)A^T-(\Gamma C+C\Gamma^T).
\]

It cancels in the transformed correction. On the fixed-noise substep,
`dC_tilde=A_tilde dP_tilde A_tilde'` and
`eta_tilde+=A_tilde eta_tilde+K_tilde m`, with the complete linked mismatch
`m=dr+H_tilde de_tilde`. The exact joint balance is

\[
 W^+-W=-\|\widetilde H\widetilde\eta-m\|_{S^{-1}}^2
       -\lambda L_P+\|m\|_{R^{-1}}^2.
\]

Thus the full covariance Fisher loss is retained, without independently
charging the AW-row norm. The remaining mean mismatch is endogenous; its
positive square cannot be treated as independent noise. In the preceding
rotating frame it is
`m_acc=m_extra-[omega]x r-Hacc Omega e-[daw]x e_theta`.
All physical/noise/reference/chart contributions remain in their literal ports.

The cost is explicit and must not be omitted. A correction changes AW by
`Delta a=K_aw r`; its next-frame map and derivative are

\[
 B_L=I+N(K_{aw}r),\qquad dB_L=N(dK_{aw}r+K_{aw}dr).
\]

These produce covariance and mean connection work in the *same* word, including
the suffix-score contribution of generated covariance. The literal OU
prediction gives
`F_tilde_aw,theta=phi[aw]x(I-W_theta)` and
`F_tilde_aw,BG=-phi[aw]x B_theta`. The integrated v/p/S rows retain
`F_y,aw[aw]x`. Do not set the literal attitude-frame discrepancy to zero or
delete BG/LIN coupling. Process covariance and reset retain both endpoint
frames and their complete derivatives.

The actual AW PSD increment is unchanged by L and has `Gamma Delta_aw=0`.
Its active face must still be evaluated using the original AW marginal
recovered from `L^-1 P_tilde L^-T`; it is not the transformed marginal's floor.
Pending target, same-history face derivative and clock remain inherited.

On the planar stratum L preserves parity, the attitude/BG covariance marginal,
physical comparison energy and magnetic innovation. The same qualified 9/10
magnetic-loss proof therefore applies to this new joint storage, as do its
pitch/BG covariance ceilings. This is a reapplication in one consistent
storage, not adding decrements from two different storage functions. Fixed
S loss is retained; noise/reference variations and held ports remain explicit.

**Minimal remaining obligation:** combine the linked mean mismatch with the
actual correction-frame jump, literal OU/reset/BG/LIN, noise/reference and AW
face ports in the suffix-score balance. Use the inherited comparison supply
budget. Prove a positive uniform transverse gap on its necessary activation
domain. That absorption, all-time domain, precision transfer and special planar
admission remain **OPEN**. No radius is solved. Structures preserved: complete
shipping state/covariance/history and all frame derivatives. Relaxation: regular
real-operation CoG scope; no new physical assumptions or runtime changes.

### Physical substitution narrows the remaining mean-port domain

**PROVED — analytical** for the fixed-input central planar physical comparison,
not the roll/BA fibre. Put `e_theta=t e_y`, nominal pitch variation `omega`,
`J_y=[e_y]x` and `R_nom'R_phys=R_y(-t)`. The sensor equation gives
`r+e_BA=-e_aw+(R_y(-t)-I)(a_physical-g)+nu`. Thus the BA estimate cancels,
and the linked mismatch becomes

\[
 m_a=J_y[t\,daw+\omega e_{aw}-\omega tD(t)f_{phys}-\omega\nu],
 \qquad D(t)=(R_y(-t)-I)/t.
\]

The continuous extension is `D(0)=-J_y`, and `||D(t)||<=1` analytically.
The source nu retains existing FAST and explicit model/arithmetic defects;
it is not a new unrestricted residual channel. Let
`k=daw-omega D(t)f_physical` and let `B_v` have pitch column k and AW block
`omega I`. For `m0=J_y B_v e`, the same-history bound is

\[
 \|m_0\|_{R^{-1}}^2
 \le V\operatorname{tr}(R^{-1}J_y B_v P B_v^T J_y^T).
\]

This follows from `e=P^(1/2)z`, `||z||^2=V`. Writing
`U=J_y'R^-1 J_y`, the coefficient is exactly
`P_pitch,pitch*k'Uk + 2*omega*k'U*P_aw,pitch + omega^2*tr(U P_aw,aw)`.
The signed pitch/AW cross covariance is retained. The source cross-plus-square
is also retained: `-2*omega*m0'R^-1 J_y nu + omega^2*nu' U nu`.
No nominal AW or BA box appears. Only comparison energy, the joint pitch/AW
covariance block, actual physical force/rotation and existing source bounds
are needed for this port. Uniform activation and absorption together with the
applied-increment frame jump remain OPEN. This is a derived bound, not a
sampled covariance enclosure or a complete-word contraction assertion.
