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

The frame continuation `app:aw-frame-word` applies (3)-(4) to both endpoint
frames of every operation. It proves exact cancellation of internal connection
terms, including generated covariance suffix scores, and derives the signed
three-dimensional AW endpoint quadratic. Its absorption remains open; it does
not delete the physical/reference/noise/reset/AW-face ports in (4).

## Actual conditional mixed-block bound

**PROVED — analytical identity and qualified inequality.** The complete
argument is `app:conditional-mixed-word`. For every regular causal word
beginning with a qualified prediction, let `E=E_aw`, `M=M_N,0`, and `q=q_0,N`
from the complete suffix formula. Set

```
C = (E' J_0 E)^-1       L = M E
T = L' J_N L           q_a = E' q
A = C^-1-T             K = T+T A^-1 T.
```

The existing conditional process loss proves `A>=c_aw C^-1>0`. Subsequent
actual dissipative base maps, including S and magnetic corrections, preserve
that inequality. This uses their actual gains and whitening; M is still only
a base transport. No complete derivative has been replaced by M.

For an unscaled symmetric basis `E_j`, write `Y=sum y_j E_j` and define

```
(G_F)_ij = tr(C^-1 E_i C^-1 E_j-T E_i T E_j)
E_q y = Y q_a.
```

The Fisher Gram is positive by the same process loss. The **actual conditional
root** mean/covariance signed block is

```
D = [[A,        -T E_q],
     [-E_q' T,  lambda G_F-E_q' T E_q]].
```

Its Schur complement is `S_lambda=lambda G_F-E_q' K E_q`. The exact test is

`D>0 <=> S_lambda>0 <=> lambda K^-1-E_q G_F^-1 E_q'>0`.

This explicit three-row test uses C,T,q rather than the unknown full signed
word gap. It resolves the root covariance-score coupling without independently
bounding q and the weak process loss. It **does not certify** a uniform margin
on shipping histories. Singular margins fail closed; no pseudoinverse deletes
forcing or an unobserved direction.

For clarity, after whitening C and diagonalizing T with eigenvalues r_i, the
reader has entries

`B_q[i,l] = delta[i,l]/2 sum_j q_j^2/(1-r_i r_j) + q_i q_l/(2(1-r_i r_l))`.

In the precisely restricted case `T=r I`, `q=(1-r)e`, the threshold reduces
to `lambda > r/(1+r) |e|^2`: the small loss cancels exactly. This checks the
same-history mechanism; isotropic OU noise does **not** make the reached
conditional covariance isotropic, so this simplification is not imposed on
the general shipping word.

### Generated ports remain in a single linked quadratic

Decompose the actual root tangent with the conditional metric:

```
u = C E' J_0 eta_0        Y = C E' J_0 dP_0 J_0 E C
eta_o = eta_0-E u         D_o = dP_0-E Y E'.
```

Let `W_o` be the complementary root storage from CA1. Define from the **full**
terminal tangent

```
t = eta_N-L(u+Y q_a)
  = M(eta_o+D_o q)+sum M_N,i+1 [v_i-U_i J_i+1 e_i+1+U_i q_i+1,N]
Z = dP_N-L Y L'
  = M D_o M'+sum M_N,i+1 U_i M_N,i+1'.
```

These packets can depend on u,Y and all inherited auxiliary/source coordinates.
They are not independent disturbances. Define

```
b = L' J_N t
V = L' J_N Z J_N L
h_j = q_a' E_j (I+T A^-1)b + lambda tr(E_j V)
P_self = t' J_N t + lambda tr(J_N Z J_N Z)
P_mix = b' A^-1 b + h' S_lambda^-1 h.
```

When the actual strict score margin holds, the exact signed completion is

```
W_0-W_N = |u-A^-1(T Y q_a+b)|_A^2
          +|y-S_lambda^-1 h|_S_lambda^2
          +W_o-P_self-P_mix.
```

In particular, with `z=(u,y)`, the new **mixed-work bound** retains half the
net conditional root loss:

`W_0-W_N >= (1/2) z' D z + W_o - P_self - 2 P_mix`.

Proof: the linear coupling vector is
`ell=(b,E_q' b+lambda [tr(E_j V)]_j)`, and exact block elimination gives
`ell' D^-1 ell=P_mix`. Complete the square with half of D. This bounds the
actual signed mixed terms **after** aggregation, preserving cancellation
between source work and internally generated covariance contributions.
It does not require an independent bound for each U, v, residual or gain.

For the AW-frame quotient, add the exact signed
`Phi_0-Phi_N-Gamma_0+Gamma_N` to the right side. No endpoint term, physical
S, source-score kernel component or gauge amplitude is discarded. The helper
requires this adjustment explicitly. The bound is a matrix inequality on a
justified lifted stratum after substituting its actual linear packet maps;
endogenous packet dependence must not be replaced by an independent maximizer.

**Remaining OPEN:** one fixed lambda and uniform positive three-row score
margin on the activated histories; domination of
`P_self+2 P_mix-Phi_0+Phi_N+Gamma_0-Gamma_N` by complementary process/S/magnetic
loss and the retained conditional half-loss. The previously proved c_aw
justifies the base inverses, not these net assertions. The inherited comparison
identity `E_loss=V_0-V_N+Supply-sum NIS` remains the score budget; no first-Live
restart or sampled NIS cap supplies a uniform bound. No radius is calculated.

Structures preserved: all 21-state/covariance coordinates, actual conditional
precision and closed-loop transports, held reduction, full suffix ports,
reference/noise/reset/AW/source chronology and endpoint gauge. Relaxations:
regular real branches and the explicit half-loss sufficient bound; no new
physical premise. D remains the unclosed uniform domination, not instability.


### No separate AW-frame penalty in the consistent packet coordinates

**PROVED — analytical:** for the actual state-dependent shear L_i, both
`L_i E_aw=E_aw` and `(L_i^-1 dL_i)E_aw=0`. The full differential therefore
preserves conditional C, dC, `E_aw' J eta`, u,Y and conditional joint storage.
Base/score conjugacy gives `M_tilde=L_N M L_0^-1`, `q_tilde=L_0^-T q`, so
**C,T,q_a and the mixed reader are unchanged**. This does not freeze dL.

The endpoint connection remains in the actual transformed packets:

```
t_tilde = L_N (t-P_N Omega_N' J_N e_N)
Z_tilde = L_N (Z+Omega_N P_N+P_N Omega_N') L_N'.
```

The complementary root storage changes by Phi_0. Hence the packet bound can
be evaluated entirely in the AW frame with the same conditional root matrix;
then only the signed quotient gauge adjustment is added. A second independent
Phi penalty would double-count the same frame work. If using input-frame
packets, the explicit Phi_0-Phi_N adjustment remains mandatory. Packet
dependence on conditional coordinates is not removed by this invariance.
This discharges separate AW-frame pricing within the consistent bound, not
uniform packet absorption or the full endpoint Schur margin.

## Complementary base packets are already paid for

**PROVED — analytical**, `app:conditional-packet-absorption`: on every regular
word covered by the conditional process lemma, let `H_i` be the full joint
mean/Fisher metric and `T_0(eta,D)=(M eta,M D M')` its actual base transport.
Take the conditional lift `R_c(u,Y)=(E u,E Y E')` and its full root-metric
orthogonal complement `R_o`. Define

```
D_b = diag(A, lambda G_F)
B   = R_c' T_0' H_N T_0 R_o
H_o = R_o' H_0 R_o
Q_o = R_o' T_0' H_N T_0 R_o.
```

Actual base nonexpansion gives `H_0-T_0' H_N T_0 >= 0`. Since the qualified
conditional process loss proves `D_b>0`, its exact Schur complement proves

`S_b = H_o-Q_o-B' D_b^-1 B >= 0`.

Consequently **base complementary packets satisfy
`P_self,b+P_mix,b <= W_o`**, with exact coefficient one, uniformly on this
regular domain. This includes all covariance cross entries, requires no new
covariance ceiling, and does not assume a positive full transverse gap.
It removes the need to charge propagation of the complementary root as if
it were newly generated forcing. It does not absorb the score or generated
ports, and does not replace the complete tangent by M.

For `v=R_c z+R_o w` and the complete remainder
`R(v,s)=actual_terminal_tangent-T_0 v`, the exact remaining balance is

```
W_0-W_N = |z-D_b^-1 B w|_D_b^2 + w' S_b w
          -2 <T_0 v,R(v,s)>_H_N - |R(v,s)|_H_N^2.
```

All `M D_0 q`, generated `U_i q_suffix`, reference/noise/reset/AW-face and
source contributions remain linked in R. In consistent AW coordinates add
only the signed quotient term `-Gamma_0+Gamma_N`. Frame work is already in
the actual coordinates and remainder. No physical gauge is deleted.

### Two algebraic scope checks change the next obligation

These are exact rational/symbolic tests of matrix implications, **not
shipping-reachable histories or numerical experiments**.

1. With `C=I`, `T=diag(1-epsilon,1/2,1/2)` and the linked score
   `q=(I-T)(0,a,0)'`, the directional critical weights are exactly
   `a^2(1-epsilon)/(4 epsilon(1+epsilon)), a^2/3, a^2/6`.
   The loss-score energy is only `a^2/2`. Thus the isotropic cancellation
   cannot be used across unequal loss directions. The fixed qualified
   `c_aw>0` remains valid; this does not assert arbitrarily small shipping
   loss or impossibility of choosing a common lambda.
2. The strictly contractive formal base
   `M=[[3/4,1/5],[1/5,3/4]]`, `P_0=P_N=I`, `q=0`, and complementary
   mean root `(0,1)` has true loss `159/400`, but the fixed half-charge
   lower bound is `-1173/21200`. Its exact base Schur remainder is positive,
   `3627/21200`. Thus fixed half-charge domination is not a necessary
   milestone. The bound remains mathematically valid.

The failed promotion of base loss to those sufficient margins is **D**.
Do not tighten the half-charge relaxation repeatedly or infer shipping
instability. The architecture review returns to the exact signed combined
work above, retaining the already-proved base absorption. No radius is solved.

| Remaining bound | Exact role | Status |
|---|---|---|
| Actual directional score reader for one lambda | Controls off-diagonal covariance-to-mean transfer across unequal losses | OPEN; inherited signed supply and source score not bounded uniformly |
| Complementary base Schur loss | Pays for base complementary packets | PROVED, semidefinite; strict transverse comparison remains OPEN |
| Combined remainder and quotient work | Must be dominated by retained squares/Schur loss with positive uniform remainder | OPEN; actual causal auxiliary/source maps retained |
| Branch/precision/entry and finite error | Activates regular identities and transfers tangent gap | OPEN; no all-time compactness or first-Live restart assumed |

Structures preserved: actual full-state base P/gain/Joseph/reset/AW chronology,
full Fisher metric, all generated suffix ports and consistent frame/gauge.
Relaxations: explicit formal matrix scope checks only; no shipping equation
or physical assumption changed. No admitted counterexample was found.

## Combined physical process score and causal generated work

**PROVED — analytical**, `app:process-source-score`. For each literal regular
prediction/addition with `C=F P F'+Q`, `e+=F e+s`, define

```
L = J-F' J+ F
score = J e-F' J+ e+.
```

When `Q t=s` is solvable, the explicit vector `z=e-P F' t` satisfies

`L z=score`, and `score' z=V-V+ + s' t >= 0`.

This controls the **combined** process/model-mismatch score. It improves the
older budget for only `L e`, followed by a separate source score. No OU law
is imposed on physical truth. An unsupported defect is rejected by the range
check and remains explicit, notably physical held-BA evolution or a chart
shift with zero process noise. No pseudoinverse removes it.

Across the actual covered word, with all unsupported local scores retained
as `q_u`, the exact inherited budget is

```
q = q_c+q_u
q_c q_c' <= B_W G
B_W = V_0-V_N + sum s' Q^dagger s
      +sum (delta' R^-1 delta-NIS) +sum_uncovered (V_next-V).
```

Every signed term belongs to the same history. Physical S contributes
`delta=-S_physical`; the estimator pseudo-update does not reset it. Compressing
this matrix bound to AW gives `q_c,a q_c,a' <= B_W A`. The actual directional
reader still uses `q_c,a+q_u,a`, including both signed cross terms. This does
not justify the failed isotropic-loss substitution or prove a common lambda.

### The derivative ports share one process Fisher loss

For SPD Q, put `F_1=dF`, `Q_1=dQ`, `s_1=ds` and form

```
Sigma = diag(P,Q)             T = [F,I]
X = [[dP,P F_1'],[F_1 P,Q_1]]
nu = s_1-Q_1 Q^-1 s
zeta = (eta-P F_1' Q^-1 s, nu)
Lambda = Sigma^-1-T' C^-1 T
f = T X Lambda (e,s).
```

Exact multiplication reproduces **all** terms of the actual prediction
`dP+=F dP F'+dF P F'+F P dF'+dQ` and
`de+=F de+dF e+ds`. In particular, `eta+=T zeta+f`.
The algebraic augmentation is not a new stochastic input or surrogate filter.
With `chi=V+s'Q^-1s-V+` and the same augmented Fisher loss `ell_P`,

`|f|_(C^-1)^2 <= (chi/2) ell_P`.

The exact signed joint balance is

```
W+-W = A_aux-ell_mean-lambda ell_P
       +2 <T zeta,f>_(C^-1)+|f|_(C^-1)^2,

A_aux = -2 eta' F_1' Q^-1 s
        +s' Q^-1 F_1 P F_1' Q^-1 s +nu' Q^-1 nu
        +2 lambda tr(Q^-1 F_1 P F_1')
        +lambda tr(Q^-1 Q_1 Q^-1 Q_1).
```

Thus the feedback square consumes `chi/2` of that same Fisher loss. The signed
cross term is **not** absorbed by this statement. Keep `A_aux`, the augmented
loss and this cross term grouped: separate maxima can destroy coordinate
cancellations. This is alternative accounting of the actual process step,
not another decrement to add to the old base loss.

### Causal coordinates and remaining threshold

All derivatives are formed through the literal auxiliary/source recurrence.
The private default tuner has its proved zero reverse MEKF edge; reference
and gate feedback do not inherit that deletion. Actual BG enters the AG
transition/noise; the same committed tuple enters both LIN transition/noise;
BA release/hold changes its actual process support. BG projection, pending AW
sync, scheduler and corrections retain their source order and separate maps.
Branch crossings require qualified finite-increment treatment.

For endpoint AW frames, transform the entire tuple:

```
dF_tilde = L+ (dF+Omega+ F-F Omega) L^-1
dQ_tilde = L+ (dQ+Omega+ Q+Q Omega+') L+'
ds_tilde = L+ (ds+Omega+ s).
```

The combined score budget is frame invariant. The full differential balance
retains the frame connections in its actual inputs; no second frame penalty
is added. Only the quotient gauge adjustment remains separate. All internally
generated covariance suffix scores survive downstream composition.

**Discharged:** separate loss/source-score pricing for supported process
forcing; independent pricing of the process covariance-induced square from
`dP,dF,dQ`; missing explicit causal process work in this accounting.

**Still OPEN:** a uniform bound on the inherited budget and uncovered scores;
affordable domination of the grouped signed auxiliary/cross work and the
remaining correction/reset/AW/gauge work by the exact Schur loss; activation
and finite-error/arithmetic transfer. No positive complete gap or radius is
claimed. The general theorem still assumes actual magnetic service; planar
all-time admission is separate.

Structures preserved: all literal states and process coefficients, physical
OU mismatch, causal auxiliary/source derivatives, supported/unsupported
forcing and endpoint frames. Relaxations: regular real-operation derivatives,
with SPD Q required only for the augmented identity. No new assumption on
physical histories, numerical exploration, or admitted counterexample.

## Uniform held score; exact remaining active-budget obligation

**PROVED — analytical**, `app:held-score-boundary`. On every reached CoG word
wholly before BA release, with constant feasible held nominal bias and frozen
actual `P_b`, all base maps and covariances split into active 18 and held 3.
Consequently the unsupported physical SLOW-bias score is exactly

```
q_u,b = P_b^-1 (b_s(t_N)-b_s(t_0)),    q_u,X = 0,
q_u,b' P_b q_u,b <= j_b min(D_s^2 T^2,4 B_s^2),  P_b^-1 <= j_b I.
```

This is uniform in time already spent in hold and in the number of operations.
It retains signed cancellation: a bias excursion returning to its initial
value has zero score. No sum of increment norms is used. With the unchanged
real core `.004` standard deviation inherited, `j_b=62500`; the 60-second
rate contribution is exactly 225. Check actual caller covariance and arithmetic
before using that profile-specific value; it is not an affordability result.

Reached tangents have block-diagonal dP and `U=diag(U_X,0)`, so both root and
generated covariance/suffix-score terms from this held score have **zero
active component**. This is also true in consistent AW coordinates. Physical
held bias still enters the accelerometer residual and `dr`, and `dP_b` still
enters `dR_eff=dR_acc+dP_b`. No active port is deleted. The result does not
cover release/A21, arbitrary masked covariance, BG projection or chart shifts.

The held tangent work itself is one signed boundary expression. With
`h=d(Delta b)-dP_b P_b^-1 Delta b`,

`W_b,N-W_b,0 = -2 eta_b,0' P_b^-1 h+h' P_b^-1 h`.

It cannot be charged to active Schur loss. Physical time-rate bounds do not
bound differentiation with respect to an arbitrary family parameter; source
tangents require their explicit causal coordinates/charge.

For the supported active process terms define, using actual same-operation
operands,

```
E_L = sum e' (J-F' J+ F) e
S_q = sum (s' Q^dagger s-s' J+ s)
C_W = sum e' F' J+ s.
B_W = E_L+S_q-2 C_W,    C_W^2 <= E_L S_q.
```

These equalities preserve the source signs. They imply that, **given a finite
uniform S_q cap**, a uniform B_W cap is equivalent to a uniform cap on E_L.
The required quantity is the process-loss projection of the inherited
comparison, not an independent box for all coordinates. Source action alone
does not bound it. The active signed telescope uses V_X, not V_X+V_b, and
retains correction defect minus NIS and every uncovered active energy change.

**Attempt outcome: D, incomplete sufficient proof.** The held unsupported
score and its repeated work are now bounded/regrouped, but no general
activated E_L cap, remaining unsupported AG/chart-score cap, or affordable
domination of grouped active process/correction/reference/AW/gauge work is
proved. First-Live energy is not a recurring bound; held LIN BIBO does not
close all those causal/precision directions. No claim that the remaining task
is merely scalar tightening, no uniform complete margin, and no radius.

Structures preserved: all states, actual reached hold/masks, effective noise,
suffix covariance terms, inherited physical SLOW/FAST and S, consistent frames.
Relaxations: regular real-operation held scope, with an explicit actual-P_b
comparison for its numeric specialization; no new physical assumption.

## Actual locked-Live reference work and paired source forcing

The actual causal reduction following the bordered audit is
`app:locked-live-causal-reference`, not another scalar-storage refinement.
For a locked-Live word with fixed inherited private/reference state and input,
the reference/offset stream is identical for all qualified MEKF-root
variations. Its homogeneous reverse port work is exactly zero. Refinement's
MEKF-dependent yaw write remains in the actual derivative, and the gravity
gate remains relevant to reacquisition outside this scope.

For source variations, the literal continuous update is
`B_ref=B_anchor+c(w-A*b)-c(w-A*b_anchor)`, where
`c(v)=(|v_xy|,0,v_z)`. Its global one-Lipschitz property proves
`|B_ref-B_anchor|<=|A*(b-b_anchor)|<=|b-b_anchor|` on the same statistics.
The default real accepted-bias cap and convex slew then give
`|B_ref|,|m-b|<=27*M_mag/20` after actual acquisition. No small-angle or
comparison-energy bound is assumed for this amplitude result.

Writing `C=Dc(w-A*b)`, `C_a=Dc(w-A*b_anchor)`, `Delta=b-b_anchor`, its
complete regular source derivative is

```
dB_ref = dB_anchor + (C-C_a)*(dw-dA*b_anchor-A*db_anchor)
         - C*(dA*Delta+A*dDelta).
```

The shared-statistics term cancels when `Delta=dDelta=0`; otherwise retain it.
The magnetic source residual and row changes remain paired and feed the full
`v_i-U_i J e_i+U_i q_suffix` expression. The homogeneous MEKF derivative is
still the complete causal one, not frozen factors. In (LR8), the terminal
source response, its signed cross term and `-Gamma_0+Gamma_N` all remain.
Uniform domination of those terms and a positive complete gap are OPEN.
In particular, on a construction-correlated lift `da_0=Lambda*v_0+zeta`, use
`M_lift=M+R_a*Lambda`; that inherited contribution remains endogenous (LR9).
The zero reverse recurrence edge does not set the total inherited derivative
to zero. A lift with no proved graph retains its full auxiliary coordinates.
This source-bound reduction neither pays the bordered restoration nor revives
that exhausted shortcut. Float accumulation and reacquisition remain open.

## Bordered comparison storage: score absorption is not uniform retention

`app:bordered-comparison-storage` gives a new analytical storage test, rather
than another independent score bound. For the actual comparison define

```
M = [[P,e],[e',kappa]],  c=kappa-e'P^-1 e>0,
F_M = tr(J dP J dP)+2 eta'J eta/c+(dc/c)^2.
```

Start the proof scalar once at `kappa_0=V_0+c_0`. It is not an estimator
variable and does not reset physical S or history. Supported prediction updates
`kappa+=s'Q^dagger s`; the lifted matrix is exactly a congruence plus a PSD
addition. An optimal correction updates
`kappa+=delta'R^-1 delta-NIS`, with the literal `delta=r+H e`; its lifted
matrix is exactly a joint covariance's Schur complement. The same Fisher
storage is nonincreasing for both substeps when these coefficient/source
operands are fixed along the variation. This absorbs their covariance-induced
score work jointly. Correction loss is the sum of the innovation Fisher loss
and the conditional-regression Fisher loss, both nonnegative.

The actual derivatives are not silently frozen: (BC4) retains the process
coefficient/source port, and (BC7) retains the full dH/dR/d(delta) connection
and noise work before norms. Reference, projection, reset, original AW face,
held effective noise, physical source and gauge contributions remain. The
fixed-operand implication is not a theorem for the entire shipping derivative.

The crucial obstruction is now exact: `c_N=c_0+B_W`. The metric's mean weight
`2/c_N` has no proved uniform lower bound. A formal rational word with
`F=1/2,Q=7/8,H=R=delta=1,s=1/2` returns `P=1/2,e=1` every word but grows
slack by `9/7`. This is an algebraic scope test, not a shipping counterexample.
Even bounded comparison energy does not stop this metric from degenerating.

Restoring slack to c_0 at a word boundary, with actual P,e unchanged, costs

```
2 B_W/(c_0*(c_0+B_W)) * eta_N' J_N eta_N
 - (dB_W)^2/(c_0+B_W)^2.
```

The negative derivative square stays linked; no independent extrema are used.
Restoration is not free. This candidate therefore does not bypass the needed
uniform comparison/Schur-work estimate. Classification D concerns the failed
promotion to uniform retention; the lift and conditional loss identities are
PROVED analytical. No radius, new physical assumption, full endogenous-work
absorption or theorem-count increase is claimed.

The exact restoration audit in `app:bordered-causal-restoration` closes a
possible loophole in this interpretation. For `lambda=c_0/2`, the full causal
word obeys

```
W_lambda,0-W_lambda,N = (c_0/2)*(A_B-F_B-R_W).
```

The negative `(dB_W)^2` in restoration cancels its terminal stored square.
It cannot be credited again to absorb causal work. Endpoint quotient energies
remain `-Gamma_0+Gamma_N` in the original storage units.

On an actual root/source lift, write covariance, mean-shear and slack tangent
maps as `X_N,Y_N,b_N`, respectively, in the endpoint precision coordinates.
Uniform equivalence with constant gamma requires exactly

```
(1-gamma*lambda) X_N'X_N + (2/c_N-gamma) Y_N'Y_N
  + b_N'b_N/c_N^2 >= 0.
```

This is a condition on the propagated causal image, not independently chosen
endpoint tangents. A two-dimensional fixed-coefficient extension of the scope
test propagates `de_0=e_2` to `de_n=4^-n e_2` with `dP_n=dc_n=0`; its exact
storage ratio is `2/(c_0+9n/7)`. Thus causal propagation alone does not imply
uniform equivalence. This example does not assert OU-III reachability or
magnetic admission. The actual shipping inequality and net signed absorption
remain OPEN. The bordered candidate has used its one motivated refinement;
further work needs a quantitative shipping mechanism, not another scalar
restoration or norm search.

## Fixed-input planar prediction: a genuine zero port, with its literal residue

`app:planar-even-process-work` proves analytically that the original even
process block has `dF_E=dQ_E=0` for same-input MEKF-root variations on the
planar stratum. The axis identities `W e_y=W^2 e_y=0` retain both literal
rotation coefficient branches and the Simpson noise construction; private
tuning is one-way. This is not true of the odd block or source/private-state
variations.

The nominal quaternion polynomial prevents the stronger claim `ds_E=0`.
Its actual normalized angle alpha has

`1-alpha'(x)=x^6(1920-80x^2+x^4)/(14745600*(w^2+v^2))`.

Thus the remaining original even auxiliary charge is exactly
`h*(1-alpha')^2/(q_g+q_b*h^2/12)*(d BG_y)^2`, below
`1e-27*(d BG_y)^2` on the stated regular real profile. The complete even
prediction still contains the combined signed work

`2 <F eta, F dP score+ds>_Jnext + |F dP score+ds|_Jnext^2`.

No separate maxima of these terms are used. Moving-frame coefficient ports
remain the exact connection differences until the whole word is composed.
This bound does not control inherited E_L, prove score/Schur domination, or
qualify the odd, reference, reset, AW, source and arithmetic contributions.
Those requested uniform bounds remain OPEN; no radius is computed.

## Direct coupled-matrix attack and the AW regression obstruction

The continuation in app:correlated-complete-gap partitions the complete
LR8/LR9 lift into endogenous root and genuine admissible physical/source
columns before forming [[G,H],[H',-S]]. All endpoint gauge cross entries
remain. An inherited auxiliary direction is endogenous unless its actual
history-origin map proves it comes from admissible exogenous forcing.
The zero reverse locked-reference edge does not establish that partition.
The full homogeneous G is tested first; completing a source square with
G alone leaves no contraction reserve. A later linked supply must use
G-theta J_root with 0<theta<c after a positive margin has been proved.

The signed AW regression work now uses the causal target-relative receipt
(CR5), not an independent AW budget. On the qualified fixed-input
covariance fibre after an actual active floor, mu=-sum w_j dI_a,j.
Each dI_cov is an exact directional reader of the SAME scalar correction's
Fisher loss, with sharp coefficient I_a*(2P_aa-I_a). Literal row/noise,
inherited target, staged-tuner and queued/process target-lag ports remain.
The face completion (CR8) keeps these readers coupled to dC and dT.
Its released-coordinate Schur test (CR9) fails analytically for nonzero
cross regression; class D, no reached shipping direction or instability.
The actual S/magnetic map (CR7) restores the omitted shared-tangent link.
Uniform domination of that full linked operator and its cross Schur cost
remain OPEN; no dependency or status count is promoted.

Appendix `app:odd-covariance-signed-gap` (OF1)--(OF1a) forms the entire
signed quadratic in the existing justified root/source coordinates. For each
process it keeps `f=T X Lambda z` and
`z'Lambda z=E_i+S_i-2C_i` inside the same polarized Fisher/cross/square
expression. The scalar-budget equivalence in (HB6) does not require bounding
E_L separately before forming that matrix. Complementary base propagation is
already paid by (CP4); generated `U_i q_suffix` and the complete source,
reference, reset, projection, AW and endpoint gauge derivatives remain.

On the qualified central-planar A21 pure odd-covariance partial derivative,
actual parity, zero odd residual and one-way private history prove zero mean
and non-AW generated self-ports. The score reader is zero on this block even
when the full E_L is nonzero. This fibre fixes inherited input/auxiliary
state; it is not substituted for their correlated variations on a general
construction lift. The 45-dimensional ambient covariance block is restricted
to the actual justified causal image.

The original AW_y floor retains cross covariance. In the existing conditional
coordinates its actual tangent is `dC_plus=-d(T B T')`. The exact same-storage
word matrix becomes `G_OO=lambda*(H_W-R_W'R_W)`, with actual transported
process/S/magnetic Fisher action and rows `d(T B T')/C_plus` at active faces.
The exact relative reader threshold is one, after subtracting the proposed
root-metric margin. No independent scalar maxima, extra base charge or new
metric is used.

**D — failed automatic absorption.** The exact formal correction--prediction--
pending-floor word (OF6)--(OF7) exceeds that threshold even with D_AA,0=0,
an optimal non-AW correction and positive process. Its offending direction is
`(-10,3,1,0,0)` in `(dB11,dB12,dB22,dT1,dT2)`. The signed gap is
`-691428110973/567390082009`. This three-state formal process is not
shipping-reached or MAGNETIC-SERVICE-admitted; no instability is inferred.
What is needed is the quantitative same-word regression reader comparison
with the actual integrated process and directional S/magnetic losses,
followed by the retained even/source/gauge cross Schur cost. Uniform margins
and the full machine-readable OPEN dependency remain unproved.
