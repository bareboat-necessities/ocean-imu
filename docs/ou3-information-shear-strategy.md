# Governing OU-III analytical continuation

This is the persistent strategy for PR #653 and subsequent handoffs. It
supersedes the order “first enclose every nominal coordinate / admit the planar
witness, then investigate general dissipation.” The shipping-faithfulness
protocol, physical assumptions, source profiles and all required CI remain.

**Order:** identify symmetry and kernels → justify a quotient or affine tube →
derive exact linked information dissipation → bound actual forcing → prove
only the domain/entry facts needed to activate that mechanism. General
stability assumes the existing MAGNETIC SERVICE premise on admitted histories.
Proving that the special planar family satisfies that premise is a separate
lemma. Neither task discharges the other.

Before a substantial calculation, state its lemma and quantifiers, the
finite-error obligation it serves, the literal source of dissipation or
cancellation, necessary proved/open bounds, and why recorded dead ends do not
apply. No new replay, secant, Monte Carlo, sampled grid or empirical contraction
campaign. Computation is exact algebra/rational arithmetic, algebraic derivative
regression, or rigorous enclosure of an already-derived compact inequality.

## Exact linked word balance — PROVED analytical identity

**Quantifiers and scope.** For every finite word on a fixed regular branch of
the real-operation shipping lift, every differentiable linked variation of its
inherited physical/estimator history, every SPD actual covariance and regular
innovation, and every fixed weight λ>0, use the inherited additive comparison
coordinate e and

\[
 J_i=P_i^{-1},\quad \eta_i=de_i-dP_iJ_ie_i,\quad
 W_i=\eta_i^TJ_i\eta_i+\lambda\operatorname{tr}(J_idP_iJ_idP_i).
\]

This is tangent storage, not yet a finite-error stability theorem. Nonlinear
remainders and binary32 defects remain separate obligations. The comparison is
not reseeded. An intermediate additive correction does not replace the literal
SO(3) injection. The covariance metric uses all symmetric entries, including
both off-diagonal contributions.

Write v=(η,vecsym dP), with the corresponding metric 𝒥. Split each actual tangent
operation *exactly* as v⁺=L v+p. Here L is the dissipative base specified below;
p is the remaining literal derivative, including endogenous dependence on v,
reference/gate history and physical source. No independence or boundedness of
p is asserted. A history-dependent frontend/reference tangent is propagated
causally or included in the augmented root; mean/P alone need not be Markov.

* Optimal unmasked acc/mag/S correction: A=I−KH, C=P−PHᵀS⁻¹HP,
  L(η,dP)=(Aη,A dP Aᵀ). Its mean loss is (Hη)ᵀS⁻¹(Hη).
  The exact shear port, with a=He+r, is
  C dHᵀR⁻¹a+K dH e−K dR R⁻¹a+K(dr+H de).
  Its covariance port is −K dH C−C dHᵀKᵀ+K dR Kᵀ.
  Thus fixed H,R and dr=−H de cancel covariance-induced dK r exactly;
  the other ports do not vanish. Solver bumps/repairs require their own
  branchwise terms; rejected corrections are identity operations.
* Reached held BA: use the proved active 18-state reduction with
  R_eff=R_acc+P_BA, retaining dR_eff=dR_acc+dP_BA and the held mean/source
  ports. The held covariance block has identity transport and zero base loss.
  No strict contraction of frozen BA is requested. This reduction does not
  cover arbitrary masked gains or unproved cross-covariance states.
* Prediction: L(η,dP)=(Fη,F dP Fᵀ), P⁺=FPFᵀ+Q. The existing prediction
  shear identity retains dF,dQ and e⁺−Fe, including physical OU mismatch.
  Its covariance port is dF P Fᵀ+F P dFᵀ+dQ. Coefficients are those of the
  same actual operation; the physical source is not forced to follow OU.
* Covariance reset: P⁺=GPGᵀ, e⁺=Ge+t, where t is the discrepancy between
  the literal injected comparison and G e. The base is the joint G congruence
  and is an isometry. Direct differentiation gives
  η⁺=Gη−GP dGᵀG⁻ᵀJe+dt−dP⁺(P⁺)⁻¹t.
  The covariance port is dG P Gᵀ+G P dGᵀ. This explicitly retains the
  difference between the finite SO(3) map and shipping covariance reset.
* Default AW: P⁺=P+U, U=E_aw Π₊(target−P_awaw)E_awᵀ≽0, e⁺=e.
  With base identity and the actual changed metric,
  η⁺=η+dP(J−J⁺)e−dU J⁺e and dP⁺=dP+dU.
  dU retains the pending target and spectral face dependence. It is not a
  fresh independent PSD choice. At a face crossing use directional/finite
  inclusion; this identity does not justify an ordinary derivative there.
* BA/BG projection: P is unchanged, e⁺=e+s, giving
  η⁺=η+ds−dP J s. Its base identity has zero loss; Euclidean ball projection
  is not asserted nonexpansive in the full covariance metric. BG projection
  precedes injection and BA projection follows it in the literal source.
* Reference, tuner and scheduler updates keep mean/P fixed until the actual
  affected operation, but propagate their auxiliary tangents and pending
  choices. The one-way default Complementary result does not remove reference
  or gate feedback. A watchdog reset or unsupported covariance replacement
  needs its literal jump map; it is not covered by the regular reset identity.

For prediction, correction and AW, put U₀=(P⁺)⁻¹/² B P¹/², where B is F,A,I
respectively. In each case D=U₀ᵀU₀ satisfies 0≼D≼I. For reset D=I. Set
E=I−D, y=P⁻¹/²η, Z=P⁻¹/²dP P⁻¹/². The base loss is

\[
 \ell(v)=y^TEy+\lambda\{2\operatorname{tr}(ZDZE)
                         +\operatorname{tr}(ZEZE)\}\ge0.
\]

Proof: the covariance output energy is tr(ZDZD); expand I=D+E.
The two traces are squared Frobenius norms of D¹/²ZE¹/² and E¹/²ZE¹/².
Consequently loss is zero **iff** Ey=0 and ZE=0, including singular D/E.
For optimal correction these conditions are exactly Hη=0 and H dP=0.
For reset all base directions have zero action. These are operation loss
kernels, not physical observability statements.

Let Q_i=𝒥_i−L_iᵀ𝒥_{i+1}L_i≽0. Expansion and telescoping, without norms, yield

\[
 W_N-W_0=-\underbrace{\sum_i v_i^TQ_i v_i}_{\mathcal A_W}
 +\underbrace{\sum_i[2(L_iv_i)^T\mathcal J_{i+1}p_i
                         +p_i^T\mathcal J_{i+1}p_i]}_{\mathcal F_W}.
\]

This is the requested balance with equality. The square charge is retained.
𝓕_W is signed linked work, **not** an independently bounded disturbance.
For the nonsheared alternative, the retained covariance/residual estimate is
||f_P||²≤NIS L_P/2 on covered corrections; its partial loss coefficient remains
(λ−NIS)/2, not λ/2. Variational NIS bounds upper-bound realized residual work,
not arbitrary-direction observability.

For exact prefix maps v_i=T_i a in one augmented root/history coordinate a,
put E_i=T_{i+1}−L_iT_i. Define

\[
 G_A=\sum_iT_i^TQ_iT_i,\qquad
 G_F=\sum_i\{(L_iT_i)^T\mathcal J_{i+1}E_i+
 E_i^T\mathcal J_{i+1}L_iT_i+E_i^T\mathcal J_{i+1}E_i\}.
\]

Then G_A−G_F=T_0ᵀ𝒥_0T_0−T_Nᵀ𝒥_NT_N, exactly. These T_i are actual
coupled derivatives, not products of frozen F/(I−KH)/G factors. The PSD action
kernel is ∩_i ker(Q_iT_i). The full gap can be indefinite even when that kernel
is trivial, because G_F has not been absorbed. Conversely pointwise positivity
of the full gap supplies no uniform margin over histories.

## Compatibility, quotient transport and the missing coercivity

**PROVED analytical distinction.** A same-record physical variation has zero
variation of the deterministic estimator and P. Its physical-error variation
can be nonzero. In A21 the physical planar BA tangent is constant while the
estimator BA transition multiplies it by φ_OU. Hence its comparison equation
carries the nonzero source (1−φ_OU)r_BA. It is not a homogeneous BA solution.
This already prevents identifying physical ambiguity with the homogeneous
word kernel. It does not exclude ambiguity or prove instability. On a regular
A21 prediction with positive full process covariance the nonzero tangent also
has positive base process loss, compensated in the forced balance.

For a branch on which an actual compatible tangent basis R_G has been proved,
use Π_G=R_G(R_GᵀJR_G)⁻¹R_GᵀJ and Π_perp=I−Π_G. Work on a fixed-rank stratum;
rank changes need separate strata or the compatible set itself. No planar line
is extrapolated to arbitrary admitted histories. With γ=Π_Gη and ξ=Π_perpη,

\[
 W_{\perp,N}-W_{\perp,0}=-\mathcal A_W+\mathcal F_W+
 \|\gamma_0\|_{J_0}^2-\|\gamma_N\|_{J_N}^2.
\]

This exact endpoint identity retains moving projectors without pretending they
commute with transport. If the mean block of the full forced word is
η_N=M_ηη η_0+M_ηP vecsym dP_0+b_W, then its transverse projection contains
Π_perp,N M_ηη R_G,0 α_0 as well as covariance and auxiliary/source ports.
It cannot be discarded. For a varying physical coordinate,
d(Π_perp e)=Π_perpη+Π_perp dPJe−dΠ_G e. Writing U=(R_GᵀJR_G)⁻¹,

\[
 d\Pi_G=dR_GUR_G^TJ+R_GU dR_G^TJ+R_GUR_G^TdJ
 -R_GU\,d(R_G^TJR_G)\,UR_G^TJ.
\]

Thus the shear quotient is not automatically the physical distance tangent;
these conversion terms, angle-to-precision scaling and nonlinear fibre
curvature g²(β⁶/36+β⁴/4) remain necessary. For general physical forcing define
compatibility through the physical same-record constraints and its causal
forced transport, rather than declaring it a homogeneous nullspace. For fixed
admissible source variation w, the tangent relation is affine in the root;
its nonlinear realization is a compatible tube.

**Precisely isolated OPEN dependency.** On each justified transverse root lift
B (including covariance and required auxiliary coordinates), define
Γ_i=T_iᵀ diag(Π_G,iᵀJ_iΠ_G,i,0)T_i and
J_perp,0=T_0ᵀ𝒥_0T_0−Γ_0. On the root subspace where BᵀJ_perp,0B is SPD,

\[
 \boxed{B^T(G_A-G_F-\Gamma_0+\Gamma_N)B
       \succeq c B^TJ_{\perp,0}B,\quad c>0.}
\]

The required c is one common constant for **every activated admitted word**,
with inherited source, covariance, clock, reference and branch state. Unstored
auxiliary directions with zero root storage must be retained as linked forcing
or given their own justified storage; they cannot be hidden in B. This is a
single linked domination obligation, not a product of independent suprema.
It is equivalent to a positive homogeneous transverse word gap on that lift.
It does not assert that a suitable B or c has been established for all histories.

Existing MAGNETIC SERVICE controls a transported, innovation-whitened physical
probe action. To discharge the boxed inequality one still needs a quantitative
comparison from those probes to the *actual coupled sheared* word action,
process/S control of the remaining transverse directions, and absorption of
G_F with the moving-gauge terms. The source contract's two-dimensional service
floor alone supplies none of those comparisons on the full joint tangent.
No new physical assumption is adopted. This missing linked comparison, rather
than a larger Cartesian reachable box, is the next mathematical target.

### Relative Schur form of the exact missing comparison

On a stratum where the two physical service root directions have independent
images in the justified transverse lift, complete them to coordinates (s,n).
If they do not, use a different stratum; do not invent a two-column lift.
Partition the *full linked transverse gap* G and root metric J in those
coordinates. For a proposed c>0 set

\[
 A_c=G_{ss}-cJ_{ss},\quad B_c=G_{sn}-cJ_{sn},\quad
 N_c=G_{nn}-cJ_{nn}.
\]

When N_c≻0, exact completion of the square gives

\[
 \binom{s}{n}^T(G-cJ)\binom{s}{n}
 =(n+N_c^{-1}B_c^Ts)^TN_c(n+N_c^{-1}B_c^Ts)
 +s^T(A_c-B_cN_c^{-1}B_c^T)s.
\]

Thus **N_c≻0 and A_c−B_cN_c⁻¹B_cᵀ≽0** suffice, and are equivalent to
G−cJ≽0 on this positive complementary-block branch. G includes all signed
shear work and endpoint gauge transport; J retains cross precision. Existing
MAGNETIC SERVICE bounds its specified probe Gram, not this Schur complement.
Even an identified principal service lower bound can be cancelled by the
B_cN_c⁻¹B_cᵀ term. The precise remaining analytical task is to control these
linked blocks on the same histories; no numerical c is supplied. This is the
existing relative-Schur identity applied to the new full shear gap, not a
new arbitrary covariance enclosure or a claim that the service premise fails.

## Positive planar measurement loss now available

[The moving-frame continuation](ou3-measurement-frame-loss.md) proves a
qualified `9/10` coupled magnetic operation loss, including its square charge,
and derives pitch/BG covariance ceilings from two strictly prior applied
observations without an inherited covariance upper bound. Its exact rational
pitch-row coercivity improves the word action; it is not a full transverse
contraction. Use the improved action and its explicit remaining signed work
in the boxed gap above. The note identifies the still-open activation facts
and distinguishes the nominal planar stratum from physical compatibility.
The controlling order and separate general/service-admission tasks are unchanged.

## Minimal activation obligations

| Required fact | Source and exact role | Current status |
|---|---|---|
| Existing actually-applied MAGNETIC SERVICE | Current source contract; supplies its specified physical probe lower action only | Assumed for the general theorem; special planar admission OPEN |
| Compatible set / fixed-rank basis | Exact physical compatibility equations; defines B and Π without deleting genuine transverse directions | Planar and quiet constructions proved on stated domains; general branchwise relationship OPEN |
| Same-word probe-to-shear comparison and signed-work absorption | The boxed linked gap; supplies one uniform c>0 | OPEN; no replay factor or private-Mahony constant substitutes |
| SPD P,S and bounded required metric comparisons | Literal process/innovation floors, existing nuisance comparisons; define storage and convert only needed rows/directions | Local identities proved; needed full inherited comparisons OPEN |
| Comparison energy and row geometry | eᵀJe cap gives shear equivalence max(2,1+2E/λ); same-P row/remainder forms bound dH and residual mismatch | Qualified first-Live central-fibre bound only; persistence OPEN |
| Reference alignment and actual gates | Connect physical service probes to nominal rows; qualify acquisition/refinement branches | Qualified real cone and captured-domain liveness retained; general activation/float accumulation OPEN |
| Injection chart and reset size | Bound t,dt,dG and projector conversion; exclude chart singularities | Exact reset calculus/remainder available; inherited prefix margins OPEN |
| AW target/face and S placement | Determine actual chronological ports, not independent phase boxes | Literal maps and one-way default tuner dependence proved; required uniform branch/finite-increment cover OPEN |
| Source/gauge/arithmetic charges | Bound linked word work after positive gap; retain physical S, SLOW+FAST primitives, C_Q α and fibre curvature | Exact source formulas retained; temporal/device/float32 qualifications OPEN |
| Activation and every-prefix containment | Finite bootstrap from construction through H18/release, without restarting | Qualified boundary and captured-domain results retained; general completion OPEN |

No separate absolute p or v box is needed merely to state this identity.
Their bounds enter only if the actual source, S residual, or chart conversion
uses them. No all-time compactness is assumed to prove its own retention.
The planar wrapper uses σ_a=.2 and adaptive S; AtomS3R uses σ_a=.12 at its reference rate and fixed S.
This note does not transfer a theorem between them.

## Complete-word score continuation

The [linked suffix-score proof](ou3-linked-word-score.md) is the current
regrouping of the word balance. Optimal additive corrections have zero local
score; prediction/AW loss and physical/source/chart mismatch remain. Covariance
variations created inside the word retain their suffix-score term. The exact
range-aware Fisher-loss charge is `(chi/2)*L_P`, with its loss-generated score
budgeted by the same-history comparison loss. This replaces independent
covariance-gain/residual charges, not the actual derivative. The uniform score
budget, generated-port absorption and full transverse action remain OPEN.
Do not count magnetic work twice or erase score forcing in an action kernel.

The [AW shear continuation](ou3-measurement-frame-loss.md) now removes the
accelerometer's nominal AW-row covariance port by an exact connection
difference. Work in one consistent joint storage: retain the transformed
mean mismatch, actual correction-frame jump and all OU/reset/AW-face ports.
Do not return to an independent per-operation AW-row norm charge, omit its
endpoint derivative, or sum magnetic losses from different storage frames.

## Scope and handoff rules

Only after the positive gap and actual charges are proved may an invariant
radius be solved. Prefix retention, H18/release, STILL/TRANSITION/MOVING
composition and full shipping arithmetic remain mandatory. Quiet absolute
entry is not a target; compatible-tube boundedness is. A distance bound does
not bound a compatible set's diameter or every absolute state error.

Separately, prove/disprove all-time applied service for the planar physical
family. Its finite floors near 7.02 are diagnostics. Full admission plus the
pair comparison rules out universal eventual retention in the specified
absolute V<.0225 ball on that profile; it neither proves divergence nor rules
out a larger absolute practical bound.

**Structures preserved:** all 21 states, actual P/K/Joseph and held masks,
physical OU/bias forcing, reset/projection/AW, frontend/reference/gates and
inherited clocks, shared history and moving physical compatibility.
**Relaxations introduced:** regular real-operation tangent scope, with every
remaining literal derivative carried as a signed port. This is an exact split,
not a surrogate trajectory. Nonsmooth, nonlinear and arithmetic closure is OPEN.
**Failure classification:** D for the still-unclosed sufficient uniform gap;
no shipping counterexample. Retained analytical identities and qualified local
results remain valid. Next: derive the linked service-to-shear comparison on
justified compatibility strata, then only its necessary activation bounds.
