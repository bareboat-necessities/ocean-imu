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

## Complete-word frame cancellation and directional endpoint absorption

**PROVED — analytical**, `app:aw-frame-word`. For every regular finite word,
with SPD actual covariances and differentiable linked frame/history, put
`C_i=L_i^-1 dL_i`, `Z_i=C_i P_i+P_i C_i'`. In the suffix-score normal form,

```
Btilde_i = L_(i+1) B_i L_i^-1
Utilde_i  = L_(i+1) (U_i+Z_(i+1)-B_i Z_i B_i') L_(i+1)'
vtilde_i = L_(i+1) (v_i+C_(i+1)e_(i+1)-B_i C_i e_i-B_i Z_i d_i)
dtilde_i = L_i^-T d_i.
```

The last term in `vtilde` is compulsory. Substitution into the **complete**
normal form, including `U_i q_(i+1,N)`, cancels all internal connections.
The terminal differential is exactly the derivative of the terminal frame.
Thus the net storage change is

\[
\widetilde W_N-\widetilde W_0
=-\mathcal A_W+\mathcal F_W+\Phi_N-\Phi_0.
\]

This equality refers to the input-frame action and work. It does **not** permit
keeping transformed improved action while deleting the associated connection
work. Cancellation is in the complete action-minus-work expression, not in
isolated nonnegative port charges. Physical mismatch, OU/BG/LIN, reference,
noise, reset/projection, AW-face and source/gauge work all remain literal.

For the AW shear alone, write `a=daw`, `C_j=N(e_j)`, `D=dP`, `J=P^-1`:

\[
\Phi=2h^Ta+a^TQa,
\quad h_j=-e^TJC_j\eta+2\lambda\operatorname{tr}(DJC_j),
\]
\[
Q_{jk}=e^TJC_jPC_k^TJe+2\lambda\operatorname{tr}(JC_jPC_k^T).
\]

`Q` is positive definite because `N(a)` is nonzero for every nonzero a and
its Fisher square is strictly positive. Exact completion gives

\[
\Phi=(a+Q^{-1}h)^TQ(a+Q^{-1}h)-h^TQ^{-1}h
\ge-h^TQ^{-1}h\ge-W.
\]

This is a lower bound; it is not an affordable upper bound on terminal work.
The actual a is inherited and cannot be replaced by the formal minimizer.
Writing `j_a=(Je)_aw`, the coefficient simplifies to

\[
Q_{jk}=\operatorname{tr}[(j_a j_a^T+2\lambda J_{aw,aw})
[e_j]_\times P_{\theta\theta}[e_k]_\times^T].
\]

This removes any separate domain obligation for every intermediate frame
increment. The remaining endpoint calculation uses **conditional AW precision**
`J_aw,aw`, the actual attitude covariance block and linked tangent rows. An AW
marginal floor does not prove a conditional-covariance floor; the existing
pitch ceiling does not bound the other attitude directions. Both parities
remain, including on the planar nominal stratum.

On a justified stored root lift let `h_i=H_i^L w`, `a_i=A_i^L w`. The endpoint
work matrices are `K_i=(A_i^L)'Q_i A_i^L+(H_i^L)'A_i^L+(A_i^L)'H_i^L`.
Let G be the exact input-frame word gap, Gamma the **transformed** gauge-energy
matrices, and Jperp the transformed transverse root metric. For a proposed c,
set

```
Z   = A_N^L + Q_N^-1 H_N^L
R_c = G + K_0 + (H_N^L)'Q_N^-1 H_N^L - Gamma_0 + Gamma_N - c Jperp.
```

If `R_c>0`, the complete remaining endpoint test is exactly

\[
Q_N^{-1}-Z R_c^{-1}Z^T\succeq0.
\]

This is a three-row Schur test, with threshold one for the squared singular
norm of `Q_N^(1/2) Z R_c^(-1/2)`. If `R_c` is not positive, retain the full
signed matrix; a pseudoinverse must not remove kernel forcing. No uniform c,
positive complementary block or threshold margin is established here.

**Dependency discharged:** separate accumulation/absorption of internal
applied-increment frame derivatives. **Remaining OPEN:** absorb the explicit
endpoint quadratic together with actual physical word work; prove the required
conditional precision, reference/chart/source activation and service-to-coupled
comparison on the same history. This is not a new source assumption, an extra
gauge quotient, a radius proof, or an all-time planar admission certificate.

**Structures preserved:** complete state/P/K/Joseph/masks, physical mismatch,
all generated covariance scores, both parity blocks and endpoint frames.
**Relaxations:** regular real-operation derivatives; the formal completion
supplies only its stated algebraic lower bound. **Failure:** D for unclosed
uniform absorption, not an admitted counterexample. No numerical campaign.

**PROVED — analytical AW boundary monotonicity.** For the actual PSD increment
Delta, let `V_a=P_aa-P_ao P_oo^-1 P_oa` and `v_a=e_a-P_ao P_oo^-1 e_o`.
The shipping AW synchronization gives exactly `V_a+=V_a+Delta`, leaving v_a
and the attitude covariance unchanged. Hence `J_aa+=(V_a+Delta)^-1<=J_aa`
and `V_cond=v_a'J_aa v_a` do not increase. The endpoint Fisher coefficient
`K(a)=tr(J_aa [a]x P_theta,theta [a]x')` also does not increase. Moreover
`a'Q a <= (V_cond+2 lambda) K(a)` by same-history Cauchy--Schwarz.
This is a monotone coefficient bound, not an extra storage decrement or an
upper bound on the signed h term. The floor is still evaluated on P_aa;
dDelta retains the original target/face derivative. Prediction/reset and
corrections still require inherited conditional-precision control.

## Conditional AW action and the endpoint Schur blocker

**PROVED — analytical, qualified:** Appendix `app:aw-conditional-loss`
splits the actual joint storage without dropping cross covariance. In the
partition `o=non-AW, a=AW`, let `B=P_oo`, `T=P_ao B^-1`,
`C=P_aa-T P_oa`, `s=J eta`, and `ell=s_o+T' s_a`. Then

```
W_a = s_a' C s_a + lambda tr(C^-1 dC C^-1 dC)
W_o = ell' B ell + lambda tr(B^-1 dB B^-1 dB)
      + 2 lambda tr(C^-1 dT B dT')
W = W_a + W_o.
```

Every optimal additive mag/S correction has `H_aw=0`, so `C,T,j_aw`
and their linked variations are invariant. Thus `W_a` is invariant, including
non-AW row/noise/reference/residual variations. The held-BA branch requires
its reached active-block reduction. Subsequent reset/frame/face operations
remain explicit. The endpoint coefficient decreases exactly by

`a'(Q-Q+)a = tr((j_aw j_aw'+2 lambda J_aw,aw) [a] K_theta S K_theta' [a]') >= 0`.

This is coefficient bookkeeping, not another loss to add to the 9/10 result.
For the actual `b=d(K_aw r)`, the endpoint work still includes
`2a'(h+-h)+2b'(Q+ a+h+)+b'Q+ b` with its signs intact.

For the homogeneous base action, conditional directions
`eta=E_aw C s_aw`, `dP=E_aw dC E_aw'` are missed by all such mag/S rows.
They are three mean and six symmetric-covariance directions, **not physical
gauges**. Literal integrated OU noise controls them: with
`Q_process >= F E_aw q_* E_aw' F'` and the existing `C<=m I`,

`A_process,conditional >= c_aw W_a`,
`c_aw=q_*/(m+q_*) > 1/100000000000000`.

The exact rational coefficient is
`248188600375173268173354737 / 2611278820575057006547769109785667114737`.
This is a uniform **base-action block** bound on the existing qualified real
isotropic regular profile, not a complete process-tangent or word contraction.
The proof shorts the full integrated OU covariance against the actual
AW-to-v/p/S column, using the existing small-argument source defect and AW
ceiling. It retains the covariance-generated prediction score; that score
belongs to the remaining signed word work. No sampled enclosure is used.

The exact conditional precision recurrence further preserves

```
phi^2 Jnext_aa = J_aa - L_aa - C_LIN
L = J - F' Jnext F
C_LIN = phi (Jnext_ao B_F+B_F' Jnext_oa)+B_F' Jnext_oo B_F
B_F = E_o' F E_aw.
```

The resulting connection-square balance contains a positive process loss,
signed integrated-chain work, and signed attitude/BG covariance work (CA7).
It does not discard the latter two terms. Source variation of phi retains
`aw*dphi` and its cross-plus-square charge.

The existing LIN path matrix also gives
`J_aw,aw <= (A_44/16) I < 15942618 I` after its qualified regular 16-second
A21 window. This includes cross covariance; a conditional precision ceiling
is now available on that domain. It does not automatically cover H18,
startup, another profile or float32.

| Required fact | Role | Current status |
|---|---|---|
| Actual P SPD, regular operation branch | Define precision and tangent | Qualified premise; global arithmetic/branch retention OPEN |
| Existing isotropic AW ceiling and LIN process range | Positive conditional base loss | PROVED analytical implication on the bound profile |
| Qualified 16-second LIN path activation | Finite conditional precision | PROVED analytical implication; chronology/other-profile transfer separate |
| Applied mag/S rows and reached masks | Complementary measurement action | Exact identities; existing service assumed for general theorem |
| Actual net mixed blocks and endpoint tangent maps | Positive complementary Schur block | OPEN |
| Reference/noise, physical/source/gauge and generated covariance work | Absorption in the same storage | OPEN |
| Finite-error chart and target arithmetic | Nonlinear every-prefix theorem | OPEN |

**OPEN — endpoint Schur margin:** the preceding facts remove a zero-action
block and establish conditional precision finiteness, but do not prove
`R_c>0` or `Q_N^-1-Z R_c^-1 Z'>=0`. In particular `G` in that expression is
already the full signed gap, not a known positive input. For net block matrix
`[[D_a,X],[X',D_o]]` in root-metric orthonormal coordinates, the exact
threshold for a proposed c is `D_a-c I>0` and
`D_o-c I-X' (D_a-c I)^-1 X >= 0`. Both are **net signed** blocks. The base loss above
cannot be substituted for that net block while discarding its work or X.
No constant-tightening campaign or radius calculation is justified by this
small positive block coefficient.

Structures preserved: all conditional covariance correlations, actual
integrated OU noise, actual mag/S gains, held reduction, AW floor, frame and
suffix chronology. Relaxations: existing regular real/profile qualifications
and a conservative lower bound for one base-action block. The incomplete
Schur domination is D, not a shipping counterexample. All-time planar
admission remains separate OPEN.


The signed mixed-block continuation is now derived in
[the complete suffix-score note](ou3-linked-word-score.md#actual-conditional-mixed-block-bound)
and `app:conditional-mixed-word`. It uses the actual conditional C, transported
precision T and linked score to give an explicit three-row root-score test,
then bounds aggregate generated-port work while retaining half the net root
loss. The endpoint adjustment remains explicit. The process lower bound
justifies the base inverses; uniform signed score/packet absorption remains
OPEN and no radius is solved.

## Complete odd-covariance gap: retain the active-face derivative

The common-tangent substitution in app:correlated-complete-gap is the
controlling calculation. Let p=P_aa, alpha_f be the actual queued target,
mu=D_aa-dalpha_f, C_plus=C+Delta, q=dC/C and a=C/C_plus.
The SAME Fisher face decrement is kept unsplit:

\[
 L_{\rm face}=q^2-(a q-\mu/C_+)^2
 +2(1/C-1/C_+)dT B dT^T.
\]

There is no division by Delta or independent minimization over q and mu.
At zero gap use the actual directional map, without asserting a smooth
crossing. Before polarization, use
the actual AA receipt (CR5). On the qualified fixed-input covariance
fibre after an actual active floor, mu is exactly the weighted sum
of prior negative information-decrement derivatives. The sharp bound
(dI_cov)^2 <= I_a*(2P_aa-I_a)*L_P uses that same actual acc/S/magnetic
correction and its directional Fisher loss. Row/noise and target/tuner
ports are retained by causal origin, including the queued target lag
and literal polynomial/PSD-repaired process Q_aa.

The S/magnetic regression identities (CR7) link dI and dC to the same
dB,dT and actual cross covariance. Releasing that link yields the
D-class failed Schur payment (CR9)--(CR10); no scalar tightening or
shipping instability inference follows. The actual linked image must
pass (OF5) and the full homogeneous cross Schur test (CR3).
Neither uniform margin is proved.

CR12 constructs all prefixes from one inherited tangent on the actual
history-image intersection; its rank and uniform parameterization are not
proved by an ambient covariance basis. CR13a pulls this unsplit formula
back to
`G_OO/lambda=H_mu+Q_mu' L_mu+L_mu' Q_mu-L_mu' L_mu`.
Its whole-word signed reserve exceeds `10^-37` on the nuisance-supported
zero-receipt intersection after 17 s regular default A21, including the first
prediction. Zero-gap faces also require zero target-relative AA tangent;
that does not follow from the strictly active receipt rows. All partial
matrices use the same zero-gap-compatible image; nonzero crossing directions
remain unresolved. Strictly inactive and absent floors add no constraint.
This is a conditional subspace result,
not uniform deficit absorption or a nonempty-image certificate.

CR17 retains the exact full cross cost `X_c' N_c^-1 X_c` and root-metric
cross entries. CR17a's lossless receipt matrix is
`Z^-1+QY+Y'Q'-E'E-Q K_sharp Q'`; it keeps the signed q/mu link and
every dependent face row. Its inverse qualification and nonnegative
uniform sign remain unproved. The complete odd blocker and homogeneous
gap are OPEN.

CR18--CR24 now pair an active face with the actual same-cycle applied acc
correction, including the intervening S/reset transport. The acc row's AW
coefficient is one. Its conditional Fisher completion pays a fraction
`eta=2*Cplus*(2*s-r-Cplus)/(s*(2*s-r))>Cplus/s>0`, with
`r=(h_o+T_acc) B_acc (h_o+T_acc)'`. The exact shared constraint is
`dC_acc=-d(beta_face)`. Combined work retains the positive completed acc
square and `eta/omega*t^2`, subtracting
`omega*(b+eta/omega*t)^2`, where omega=1-eta, b=d(beta_face)/Cplus and
t=dT_acc B_acc(h_o+T_acc)'/Cplus. This is the original Fisher storage.

The complete weighted reader threshold (CR22) and full cross remainder
remain unproved. Actual process/S/magnetic losses stay in H_acc and pay
their transported covariance columns; the acc loss is counted once.
At h_o+T_acc=0 a weighted T dB T' reader survives. CR24's actual integrated
cross-covariance generator, not a generic positive-noise model, must link
that reader to service loss. Neither the pointwise fraction nor CR23's
possibly empty constrained-kernel reserve closes the AW blocker. No uniform
margin, inherited origin/rank, source, activation or arithmetic promotion.

Appendix `app:odd-covariance-signed-gap` specializes the **same** fixed-weight
complete signed gap to a regular central-planar A21 pure odd-covariance fibre,
with inherited delivered/private state fixed. Literal parity and zero odd
residual give zero mean/reference/score and non-AW generated self-ports.
This is a partial derivative on the justified causal image, not a substitute
for the full correlated root/source lift or for arbitrary odd mean variations.

For scalar AW_y, the existing conditional decomposition gives
`beta=T B T'` and `C_plus=C+Delta`. The actual active floor retains B and T
and gives `dC_plus=-d(beta)`. Its signed Fisher decrement is

```
(W-W_plus)/lambda = (dC/C)^2
  + 2*(1/C-1/C_plus)*dT B dT' - (d(beta)/C_plus)^2.
```

After all literal covariance prefixes, including every active AA deletion,
the complete odd partial matrix is `lambda*(H_W-R_W'R_W)`. H_W contains
actual process, acc, recurring S and applied magnetic Fisher losses plus the
positive face terms; R_W contains only the remaining regression rows.
For a scalar applied mag/S row h the directional Fisher action is exactly
`(D h')'*(2 J/S-h'h/S^2)*(D h')`, with the same actual covariance/noise and
chronology as the innovation-whitened service probes. This retains the
transported weak directions and does not promote a two-column physical
service floor to a full covariance floor.

If `H_W-c H_0>0`, the exact uniform-margin candidate requires
`I-R_W*(H_W-c H_0)^-1*R_W' >=0`; its threshold is one. A formal anisotropic
correction--prediction--floor scope test exceeds one even with fixed root
AA tangent and positive process. Classification **D**, not a reached
shipping counterexample. Conditional precision nonincrease cannot pay the
face derivative automatically. The actual reader comparison, full even/source
cross Schur work, inherited activation and uniform complete margin remain
**OPEN**. Existing frame cancellation, magnetic loss and conditional process
results are unchanged. No new storage, physical assumption or radius.
