# Linked suffix score: cancellation before the complete-word bound

This is a continuation of the [governing strategy](ou3-information-shear-strategy.md),
not a new contraction search. The authoritative proof is
`app:linked-word-score`; `information_shear_word.py` checks its algebra with
exact rationals. It extends the existing `planar_cell_transport` entry without
promoting its open uniform premises.

## Precise lemma and role

**PROVED — analytical.** For every finite regular real-operation word on one
inherited shipping history, split each operation into its dissipative base and
its **exact linked** differential ports as below. The suffix-score formula is
an identity for its complete tangent, including covariance variations created
inside the word. No bound on nominal position/velocity, residual/NIS or initial
covariance is needed for the identity. SPD boundary covariances and valid local
charts/branches are needed. Their all-time activation is not asserted here.

The new cancellation eliminates the explicit covariance-dependent gain times
residual at every optimal additive correction from the *local score*. Its
influence through subsequent comparison states and generated ports remains.
Thus one need not charge a new independent NIS maximum at every correction.
This does not eliminate accelerometer row feedback, reference/noise variation,
physical OU mismatch, reset or AW feedback.

## Exact normal form

Write `D_i=dP_i`, `J_i=P_i^-1`. For a literal operation choose base `B_i`, and
retain its differential exactly as

\[
 D_{i+1}=B_iD_iB_i^T+U_i,\qquad
 de_{i+1}=B_i de_i+B_iD_i d_i+v_i.                 \tag{1}
\]

`U_i,v_i` are functions of the actual complete inherited variation, including
auxiliary state and source history. They are not independently variable noise.
For an optimal additive correction `B=A=I-KH`, `C=AP`,

\[
 d=H^TS^{-1}r,\quad U=-K\,dH\,C-C\,dH^TK^T+K\,dR\,K^T,
\]
\[
 v=dK_{H,R}r+K(dr+Hde),\quad
 dK_{H,R}=[P\,dH^T-K(dH\,PH^T+HP\,dH^T+dR)]S^{-1}.
\]

This retains the actual innovation derivative and the gain/residual square
charge when storage is formed. Reached held-BA corrections use the active
18-state reduction, `R_eff=R_acc+P_BA`, with its derivative and held-bias ports.
Arbitrary masked updates do not satisfy this optimal formula. The full held
coordinates remain in the word; no strict held-BA contraction is asserted.

Prediction has `B=F`, `e+=Fe+s`,
`U=dF P F'+F P dF'+dQ`, `v=dF e+ds`. Reset has `B=G`,
`P+=GPG'`, `e+=Ge+s`, `U=dG P G'+G P dG'`, `v=dG e+ds`.
Default PSD AW addition has `B=I`, `d=0`, `U=d(increment)`, `v=0`;
projection has `B=I`, `U=0`, `v=d(shift)`. Nonlinear attitude/chart mismatch
is included in `s,ds`, not equated to the covariance reset. Moving-frame
connection and branchwise solver/noise derivatives belong to the same ports.
A non-additive correction comparison carries its extra `s,ds` too.

Define only a **base** transport `M_{b,a}=B_{b-1}...B_a`, never an approximation
to the full derivative. Set

\[
 a_i=J_i e_i+d_i-B_i^TJ_{i+1}e_{i+1},\qquad
 q_{a,N}=\sum_{i=a}^{N-1}M_{i,a}^Ta_i.
\]

The score telescopes:

\[
 q_{a,N}=J_a e_a+\sum_{i=a}^{N-1}M_{i,a}^Td_i
                  -M_{N,a}^TJ_N e_N.                       \tag{2}
\]

At an optimal **additive** correction, `A'C^-1=J` and
`JK=H'S^-1`, so `a_i=0` for every residual, without assuming `dH=dR=0`.
Those derivatives remain in `U,v`. At prediction and AW respectively,

\[
 a_i=(J_i-F^TJ_{i+1}F)e_i-F^TJ_{i+1}s_i,
 \qquad a_i=(J_i-J_{i+1})e_i.
\]

Reset contributes `-G'Jnext s`; projection contributes `-J s`.
A correction with extra comparison mismatch contributes `-A'C^-1 s`.
Thus physical forcing is not silently removed by the cancellation.

For `eta_i=de_i-D_i J_i e_i`, variation of constants in (1) gives

\[
 D_N=M_{N,0}D_0M_{N,0}^T+
       \sum_i M_{N,i+1}U_iM_{N,i+1}^T,                       \tag{3}
\]
\[
 \boxed{\eta_N=M_{N,0}\eta_0+M_{N,0}D_0q_{0,N}
 +\sum_i M_{N,i+1}[v_i-U_iJ_{i+1}e_{i+1}+U_iq_{i+1,N}].} \tag{4}
\]

The last `U_i q_(i+1,N)` is essential: covariance variation created by a
changing accelerometer row/noise, reset or AW face is transported through
later gains. Omitting it is not justified by cancellation of the local score.
Equations (3)-(4), inserted together into the joint storage, give the existing
exact signed cross-plus-square word work. Do not take separate global maxima
of these terms or combine this accounting with a second charge for the same
magnetic feedback already absorbed by the moving-frame lemma.

## Linked Fisher-loss charge, including singular action

**PROVED — analytical.** For every `P,C>0`, base transport `M` with
`G=P^-1-M'C^-1M>=0`, symmetric `X`, and score `q` in `range(G)`, choose any
solution `Gz=q` and put `chi=q'z`. Then

\[
 \boxed{\|MXq\|_{C^{-1}}^2\le\frac{\chi}{2}L_P(X)},\qquad
 L_P(X)=\operatorname{tr}(JXJX)
 -\operatorname{tr}(C^{-1}MXM^TC^{-1}MXM^T).              \tag{5}
\]

Proof: set `D=P^(1/2)M'C^-1MP^(1/2)`, `E=I-D`,
`Z=P^-1/2 X P^-1/2`, `qbar=P^(1/2)q`. The range hypothesis is exactly
`qbar in range(E)`, and `chi=qbar'E^dagger qbar`. Since `0<=D<=I`,

\[
 L_P=2\operatorname{tr}(ZDZE)+\operatorname{tr}(ZEZE)
 \ge2\|D^{1/2}ZE^{1/2}\|_F^2.
\]

Write `qbar=E^(1/2)w` with minimum `|w|^2=chi` and apply Cauchy--Schwarz.
This proves (5), including noncommuting covariance variations. If `q` is not
in `range(G)`, a pseudoinverse would discard forcing in a zero-action direction.
The exact-rational helper instead rejects the range failure. This is not a
shipping counterexample and does not authorize deleting a physical gauge.

For the *loss-generated* score, let
`Q_i=J_i-B_i'J_(i+1)B_i>=0` on prediction/AW operations and let

\[
 q_{loss}=\sum_{i\in\mathcal L}M_{i,0}^TQ_i e_i,
 \qquad E_{loss}=\sum_{i\in\mathcal L}e_i^TQ_i e_i.
\]

Stack `A_i=Q_i^(1/2)M_(i,0)` and `w_i=Q_i^(1/2)e_i`.
Then `q_loss=A'w`, `A'A<=G`, hence `q_loss in range(G)` and

\[
 \boxed{\chi(q_{loss})\le E_{loss}.}                       \tag{6}
\]

For example, using `Gz=q_loss`,
`chi=(Az)'w<=sqrt(chi)*sqrt(E_loss)` proves the bound.
This is a same-history nonmeasurement comparison-loss budget; it is not a
new box on every nominal coordinate. The remaining source/chart score is
`-sum M_(i,0)'B_i'J_(i+1)s_i`, including physical BA OU mismatch. Retain its
sign, range and kernel components. No uniform bound on it is asserted.

## Inherited comparison energy closes the loss-score bookkeeping

**PROVED — analytical.** For every covered finite word and each of its
prefixes, let `V_i=e_i'J_i e_i`. At an optimal additive correction define
`delta_i=r_i+H_i e_i`. Since `J+=J+H'R^-1H`, `J+ K=H'R^-1` and
`K'J+ K=R^-1-S^-1`, expansion gives

\[
 V_{i+1}-V_i=\delta_i^TR_i^{-1}\delta_i-\mathrm{NIS}_i.
\]

At a noncorrection step, `e+=B e+s` gives

\[
 V_{i+1}-V_i=-e_i^TQ_i e_i+
 2(B_i e_i)^TJ_{i+1}s_i+s_i^TJ_{i+1}s_i.
\]

Let `Supply_W` sum the actual correction defect energies and these **signed**
noncorrection cross-plus-square terms. Reset/projection have zero base loss;
include them in the same sum. Telescoping proves

\[
 \boxed{\chi(q_{loss})\le E_{loss}
 =V_0-V_N+\mathrm{Supply}_W-\sum_i\mathrm{NIS}_i
 \le V_0+\mathrm{Supply}_W-\sum_i\mathrm{NIS}_i.}       \tag{8}
\]

Thus the loss-score budget needs root comparison energy and actual signed
supply, not independent bounds on every intermediate state. The same identity
controls terminal comparison energy after retaining nonnegative `E_loss`.
Realized NIS is a dissipation term in this identity; it is still **not** an
observability lower bound on arbitrary perturbations.

At the shipping `S=0` update, `r=-S_hat` and `e_S=S_hat-S_physical`, hence
`delta=-S_physical`. Its supply is not zero and physical S is not restarted.
Keep correction defect energy and NIS linked. Accelerometer/magnetic defects
retain physical/model/reference mismatch. Reached held BA uses effective noise
and the constant held comparison metric at correction; its physical evolution
remains in prediction forcing.

The first-Live bound is usable only at that actual central-fibre boundary
under its existing hypotheses. Applying it anew at each future word would be
an unjustified restart. Uniform signed source supply, the remaining source
score and generated ports remain open. Equation (8) proves the reduction of
the activation obligation, not its all-time discharge.

## Precisely isolated coercivity obligation

For the root terms of (3)-(4), suppose a justified space has
`G>=cJ`, `0<c<1`, and the **combined** score has `chi<=chi_*`.
These are currently open uniform premises, not assumed new physical laws.
The root joint storage loss is bounded below by the quadratic form in
`(sqrt(eta'Jeta),sqrt(L_P))` with matrix

\[
 \begin{pmatrix}
 c&-\sqrt{(1-c)\chi_*/2}\\
 -\sqrt{(1-c)\chi_*/2}&\lambda-\chi_*/2
 \end{pmatrix},\quad
 \det=c\lambda-\chi_*/2.                                  \tag{7}
\]

Thus `chi_*<2c lambda` gives a positive root-block loss. This is a
**CONDITIONAL implication**, not the full word gap: covariance loss must also
control the retained covariance directions, and all generated ports in (4)
and (3) must be included before applying it to shipping. On a physical
quotient, moving projector work, endpoint gauge energy and `C_Q alpha` remain.
There is no claim that the physical fibre is the kernel of `G` or the coupled
action. Frozen held directions need tube treatment, not strict contraction.

| Minimal obligation | Exact role | Status |
|---|---|---|
| Regular SPD operation and reached held manifold | Identities (1)-(4), loss PSD | Analytical branch identities proved; inherited all-time branch/precision activation OPEN |
| Comparison loss on prediction/AW and physical mismatch score | Bounds (5)-(6); possible range/kernel forcing | Root-energy plus signed-supply budget proved; uniform supply and source conversion OPEN |
| Actual acc/reference/noise/reset/AW coefficient derivatives | Generated ports `U,v` and their suffix scores | Exact derivatives retained; uniform signed absorption OPEN |
| Existing magnetic service plus complementary process/S action | Quantitative `c`, covariance action and justified quotient | Magnetic pitch-row loss qualified; full transverse comparison OPEN |
| Reference/chart/energy and event margins | Activate the retained magnetic lemma and bound only the ports used | Qualified local implications; every-prefix activation OPEN |
| Physical fibre, moving projector and finite precision | Gauge/tube and finite-error transfer | Physical compatibility proved on its domain; general strata and full transfer OPEN |

No all-state reachable box, fixed-root covariance ball, per-AW norm product,
finite replay constant, independent K/r/P bound or private Mahony contraction
is used. No invariant radius or all-future service floor is solved prematurely.
The planar wrapper and AtomS3R remain distinct profiles. General dissipativity
assumes existing service; this result does not admit the planar witness.

**Structures preserved:** complete state/P/K/Joseph, reached held masks, actual
row/noise/innovation derivatives, covariance created inside the word, reset,
projection, AW faces/targets and inherited physical/source/reference clocks.
**Relaxations introduced:** regular real-operation identities, followed only
by the explicit Cauchy--Schwarz sufficient bound (5) and conditional bound (7).
No shipping equations or assumptions are changed. **Failure classification:**
D for the still-unclosed uniform complete-word sufficient gap; no new failed
shipping execution and no admitted counterexample. **Next calculation:** use
(3)-(4) to form the *combined* acc/AW-row and prediction/reset coefficient work,
then bound its signed action and comparison-loss/source score on the minimal
inherited activation domain. The missing uniform margin is not a NIS maximum.
