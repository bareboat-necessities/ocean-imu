# OU-III: linked soft return with the full covariance baseline

## Role and status

This calculation supplies the scalar O2 return used by Corollary K in
`ou3-corrected-word-proof.md`, section 7. The desired complete-word inequality
remains `V_next <= rho V + c_d ||d||^2`, `rho < 1`, with capture, coercivity,
nonlinear remainders and every-prefix retention proved separately. Contraction
words remain 100 s; 17 s is the nuisance/root warm-up. Nothing here changes the
estimator, the coupled tau/sigma_aw/R_S/T_S chronology, or the physical domain.

The Schur scalar return and its fixed-point criterion are already present in
sections LK/AL/DC/DP of the corrected-word proof. The new contribution is the
**full-baseline adjacent-word rank-one identity and its sharp linked generalized
eigenvalue**, including the correlated two-dimensional closed form below.

This note corrects the subsequent LP10 and LP14--LP22 application in
`ou3-aw-adjoint-cancellation.md` and the claims that the synthetic adjacent-line
update establishes shipping kernel invariance. A rank-one information bound is
not a rank-one covariance. The next ceiling is not a new applied measurement.
The old scalar algebra is valid only for its explicitly synthetic model.
All examples below are exact matrix examples, not shipping-reachable witnesses.

## 1. The product does not cancel merely because Pi is shared

Use one fixed positive physical metric and its dual throughout. Express the
following matrices in the resulting normalized coordinates; do not normalize
an attitude/BA mixture differently at successive boundaries. For one realized
word and the actual next soft functional n, write

    d_W = n' Pi_W n,
    H_W = lambda_max(Phi_tilde_W' Pi_W^-1 Phi_tilde_W).

Then `d_W H_W` is a same-word product, whereas `dbar Hbar` also separates the
suprema over words. Neither is generally the direct return. Cauchy gives

    |n' Phi_tilde x|^2 <= (n' Pi n) x' Phi_tilde' Pi^-1 Phi_tilde x.

This is a lower comparison on the product, not a way to upper-bound it by the
left side. Nor may Pi be compressed before inversion: in coordinates S/Sperp,

    Pi = [[A,B],[B',C]],
    y in S => y' Pi^-1 y = y_S' (A-B C^-1 B')^-1 y_S.

The inverse of the compression A is wrong unless the cross block vanishes.
For Pi=[[1,4/5],[4/5,1]], S=span(e1), n=y=e1, the compressed condition number
is 1 but `(n'Pi n)(y'Pi^-1 y)=25/9`; this tends to infinity as the correlation
tends to one. Thus LP10's compression-only condition-number bound is invalid.

## 2. Exact scalar return and dual certificate

Freeze one literal realized word, with all actual source correlations and
nuisance elimination retained as in Theorem D. If the actual root satisfies
`nu' P0 nu <= c`, its information obeys `P0^-1 >= nu nu'/c`. Consequently

    G_c = J + nu nu'/c,
    P_soft(c) = Pi + Phi_tilde G_c^-1 Phi_tilde',
    D_W(c) = n' P_soft(c) n
           = n' Pi n + w' G_c^-1 w,       w=Phi_tilde' n.              (LS1)

Here G_c must be positive definite on the retained coordinates. The diffuse
quotient and the known-root covariance Pi have NOT disappeared. In particular,
no inverse of Pi is needed to evaluate this scalar formula.

The exact generalized-eigenvalue and dual identities are

    w' G_c^-1 w = sup_(z!=0) (w'z)^2/(z'G_c z)
                = max_z (2w'z-z'G_c z)
                = lambda_max(w w', G_c).                            (LS2)

Completing the square proves the middle equality; Cauchy in the G_c metric
proves the quotient formula, attained by z=G_c^-1 w when w!=0. Therefore

    D_W(c) <= c_next
    iff [[c_next-n'Pi n, w'], [w, G_c]] >= 0.                         (LS3)

This is the linked certificate, not a product of separate maxima. It applies
also when Pi is semidefinite, provided G_c is positive definite; any additional
regularity needed by the contraction theorem remains a separate premise.

For completeness retain the general, non-eigenvector Schur reduction. Set
`s=||nu||^2`, take the first basis vector `nu/sqrt(s)`, and write

    J=[[a,h'],[h,B_Q]],        w=(beta,g),        B_Q>0,
    j=a-h'B_Q^-1 h >=0,
    ell=beta-h'B_Q^-1 g,
    d_perp=n'Pi n+g'B_Q^-1 g.

Block inversion gives exactly

    D_W(c)=d_perp+c ell^2/(s+c j).                                   (LS4)

For a unit least eigenvector, s=1, j=lambda_1(J), and ell=n'Phi_tilde nu.
All quotient uncertainty remains in d_perp. A fixed ceiling obeys

    D_W(c)<=c
    iff j c^2+(s-ell^2-j d_perp)c-s d_perp >=0.                       (LS5)

When j=0 and d_perp>0, this is possible precisely when
`alpha=ell^2/s<1`, with `c>=d_perp/(1-alpha)`. If d_perp=0, alpha=1 is the
non-strict equality case and must not be called strict contraction. When j>0,
the positive root of LS5 is the exact threshold. The derivative is
`D_W'(c)=s ell^2/(s+c j)^2`; it is not the complete-state energy contraction.

Example: Pi=diag(100,1), Phi_tilde=I/2, J=diag(0,1), nu=n=e1. The SAME-WORD
product d_W H_W=25 fails the old sufficient test, but D_W(c)=100+c/4 and
D_W(200)=150<200. No source-uniform claim follows from this matrix example.

Conversely Pi=I, Phi_tilde=I, J=diag(0,1), nu=n=e1 gives D_W(c)=1+c, which
fails every finite ceiling. Updating a bare rank-one covariance with an added
next precision 1/c instead gives c/2. That latter number answers a different,
fictitious-observation problem and cannot prove LS3.

## 3. Keep the background while transporting the one soft excess

After the current word, the spectral version of LS1 splits into

    P_soft = B + p u u',      B>0, ||u||=1, p>=0.                     (LS6)

B includes Pi AND the diffuse-quotient term, with their cross-covariances.
The physical direction and p come from that same word; neither is reseeded.
More generally write the excess a y y' without normalizing y.

Let J_plus>=0 be the next word's REAL root observation information; it includes
only the actual applied data and their literal source elimination. Define

    C=(B^-1+J_plus)^-1,
    L=(I+B J_plus)^-1=I-C J_plus,
    H=J_plus-J_plus C J_plus >=0.                                    (LS7)

The exact conditioned covariance is

    [(B+a y y')^-1+J_plus]^-1
      = C + a (L y)(L y)'/(1+a y'H y).                               (LS8)

Proof: Sherman--Morrison first gives `(B+a yy')^-1`. Apply it again after
adding J_plus. The numerator is `C B^-1 y=L y`; the denominator simplifies
using `B^-1-B^-1 C B^-1=H`. Equivalently factor J_plus^(1/2) to see

    H=J_plus^(1/2)
      [I+J_plus^(1/2) B J_plus^(1/2)]^-1 J_plus^(1/2).

This also proves H>=0 and ker(H)=ker(J_plus), including singular J_plus.
For the complete next-word terminal covariance, replace C by
`Pi_plus+Phi_tilde_plus C Phi_tilde_plus'` and L y by
`Phi_tilde_plus L y` in LS8. The extra denominator is unchanged. Thus LS8
composes the full word, not merely a hand-selected measurement.

Importantly, the raw quotient gap is not the effective scalar action: B
appears inside H. It represents uncertainty that can mask a root perturbation.
The proof ceiling 1/c_next is not added to J_plus while checking that ceiling.
Doing so would establish invariance only after a measurement the filter never
performed and double-count the desired conclusion as information.

## 4. Sharp linked generalized-eigenvalue identity

Fix B, J_plus, p and the next scalar functional n. Let E be a full-column-rank
basis for the ALLOWED propagated soft subspace, not necessarily orthonormal.
For `u=Ez`, `||u||=1`, set

    v=E' L' n,
    M_p=E'E+p E'H E >0.

LS8 and one generalized Rayleigh quotient give

    sup_(u in Range(E), ||u||=1)
      n' [(B+p uu')^-1+J_plus]^-1 n
      = n'C n + p v' M_p^-1 v.                                     (LS9)

Indeed the excess is `p(v'z)^2/(z'M_p z)` after homogenization. It is attained
by z proportional to M_p^-1 v, with normalization in E'E. If v=0, every allowed
direction has zero excess in this functional. The result is basis-invariant:
replacing E by E R, with R invertible, does not change it.

For a single soft line LS9 is scalar. A two-dimensional admissible soft chart
requires only a 2x2 inverse. At a higher-multiplicity least eigenspace retain
its actual full dimension; nullity<=1 does not by itself bound the multiplicity
of a positive least eigenvalue. Maximizing over an enlarged subspace gives a
valid relaxation, not an assertion that its maximizer is shipping reachable.
When the actual direction is known, LS8 is sharper still and should be used.

There is an exact full-space dual identity:

    sup_(||u||=1) n'[(B+p uu')^-1+J_plus]^-1 n
      = n'[(B+p I)^-1+J_plus]^-1 n.                                 (LS10)

To prove it, put E=I in LS9 and use the matrix identity
`[(B+pI)^-1+J_plus]^-1=C+p L(I+pH)^-1 L'`. Thus an isotropic covariance
upper envelope is attained by a rank-one addition for this SINGLE functional;
it is not attained by that same rank-one direction for every functional.
This is the requested same-word duality: baseline, transfer and information
remain linked until after the maximization, without dbar, Hbar or a condition
number estimate for Pi.

## 5. Exact 2x2 cancellation and the true worst angle

In coordinates {next kernel n, next quotient q}, take unit n=e1, q=e2,
`J_plus=j q q'`, `j>=0`, and write

    B=[[b_nn,b_nq],[b_nq,b_qq]]>0,
    eta=j/(1+j b_qq),
    u=(cos(theta),sin(theta)).

Then the total next-kernel variance is exactly

    b_nn-eta b_nq^2
      + p (cos(theta)-eta b_nq sin(theta))^2
          /[1+p eta sin(theta)^2].                                 (LS11)

The numerator is correlated with the same covariance decrease in the first
term. LS9 gives the SHARP total maximum

    max_theta variance
      = b_nn+p-eta b_nq^2/(1+p eta)
      = b_nn+p-b_nq^2/(b_qq+p+1/j)       (j>0).                      (LS12)

At j=0 the continuous formula is b_nn+p. A maximizing angle satisfies
`tan(theta)=-eta b_nq/(1+p eta)`; with nonzero correlation and information,
coincident lines are not the maximizer. The excess alone can exceed p even
though the TOTAL remains <=b_nn+p. Bounding its increase separately from the
baseline decrease would again destroy the cancellation.

For B=[[1,-9/10],[-9/10,1]], j=p=1, the conditioned baseline is 119/200,
the largest excess is 227/200, and the sharp total is 173/100. The coincident
line gives only 319/200. The rational intermediate direction (12/13,5/13)
already has excess 1083/968>1 and total larger than the coincident line.
The old no-background angle monotonicity does not extend to a full covariance.

## 6. What remains to close the theorem

The current O2 obligation is still LS3 on every admitted same-history adjacent
100-s pair, or a proved blockwise counterpart. LS9/LS12 are sharp relaxations
for its actual background and allowed soft directions; they are not a numerical
certificate of a shipping ceiling. Finite deterministic gain and an ordered
quotient gap alone do not exclude the exact-kernel case alpha>=1 in LS5.

On a fixed compact adjacent-pair class, independent of the candidate ceiling,
with continuous LS4 data, s uniformly positive and finite d_perp, the existing
DC argument reduces large-c strict closure to the exact-kernel subset j=0:
`max_(j=0) ell^2/s<1`. This is conditional; it does not prove compactness,
boundary nullity, that strict inequality, or retention. No arbitrary successor
line may replace the one from the same physical continuation.

The next analytical target is the linked LS3/LS8 functional on literal
compatibility transport, including the background quotient covariance. Either
prove its uniform ceiling directly or derive a linked block return. Do not
retry the independent dbar Hbar product, insert 1/c_next as real information,
or label the algebraic counterexamples physical trajectories. Once O1 and O2
are both established, Corollary K supplies the linear word bound; nonlinear
and arithmetic margins still require their existing separate obligations.

## Validation

`LinkedSoftFullCovarianceTests` in
`tests/validation/test_ou3_rank_loss_interval_factor.py` checks twelve exact
rational identities/counterexamples, including singular information, Schur
return, chart scaling, attainment, full-space duality and LS12. These checks
validate algebra, not admissible-history coverage or outward-rounded constants.
The helper `linked_product_status()` keeps source_uniform_O2_closed=false.


## Literal compatibility-graph alpha and finite-block return

The exact-kernel large-c slope can be written without choosing a soft
eigenvector. On a rank-one exact compatibility stratum use the homogeneous
literal graph

`r_W(a)=(a,0_{LIN/AW},-A_W a)`,                           (CG-1)

where `a` is a nonzero homogeneous attitude coordinate and `A_W a` is the
BA vector forced by every literal magnetic/accelerometer/S/process zero-action
compatibility equation of that SAME realized word. `A_W` therefore depends
on the actual chronological nominal-force rows, resets, held/active BA
transport and the committed coupled tuner history; it is not freely chosen.

Fix the physical metric `M>0` once and define

`s_W=r_W' M r_W`, `e_W=r_W/sqrt(s_W)`.                  (CG-2)

For a same-history successor `W_+` let `r_+` be its literal graph
generator and `n_+=M r_+/sqrt(s_+)` the corresponding unit dual functional.
On an exact-kernel word `J_W r_W=0`. In the Schur basis whose first vector
is `e_W`, positive semidefiniteness gives both the first diagonal and the
cross block of `J_W` equal to zero. Thus `j_W=0` and the Schur coefficient
in DC-3 is not arbitrary:

`ell_W=e_W' Phi_tilde_W' n_+
       = r_+' M Phi_tilde_W r_W / sqrt(s_W s_+)`.          (CG-3)

Consequently the exact large-c slope is

`alpha_W=ell_W^2
 = |r_+' M Phi_tilde_W r_W|^2/(s_W s_+)`.                 (CG-4)

If one retains the unnormalized convention of LS4, the same statement is
`alpha=ell_raw^2/s_W`; CG-4 is the invariant normalized form. On a
zero-action compatibility trajectory `Phi_tilde_W r_W=T_W r_W`, because the
conditional root map and literal deterministic homogeneous transport agree
on the data-null root. Therefore

`alpha_W
 = |r_+' M T_W r_W|^2/
   [(r_W'Mr_W)(r_+'Mr_+)]`.                               (CG-5)

This is the literal compatibility-graph transfer, not a least-eigenvector
surrogate. Every dependence of `A_W,T_W,r_+` on physical acceleration,
attitude, biases, S pseudo-updates, magnetic service, covariance-generated
gains and the coupled `tau,sigma_aw,R_S,T_S` chronology remains inside CG-5.

Hence the one-boundary theorem question is precisely

`alpha_bar_1 =
 sup_(same-history exact-kernel pairs)
 |r_+' M T_W r_W|^2/(s_W s_+) < 1 ?`.                    (CG-6)

No current lemma proves CG-6. BA decay alone cannot: the next graph can change
its attitude/BA ratio so that `T_W r_W` is collinear with `r_+`. Conversely
compactness does not prove equality; it only ensures that the maximum is
attained once the closed same-history pair class is established.

### Full-baseline block return

If `alpha_bar_1=1`, do NOT conclude instability. The correct m-word object
is the complete composed Riccati word, not the numerical composition of the
scalar functions `D_W`. Scalarization after each word discards the baseline
matrix B and its cross-covariances, which LS8--LS12 show can change the next
return.

For a same-history block
`B_m=W_(j+m-1) o ... o W_j`, compose the literal chronological factors
first, preserving the actual carried covariance, source columns, tuner state
and shared boundary states. Let

`Pi_[j,m], J_[j,m], Phi_[j,m]`

be Theorem-D's known-root covariance, root information and conditional
root-to-terminal map of that WHOLE block. With endpoint graph generators
`r_j,r_(j+m)`, define

`G_[j,m](c)=J_[j,m]+r_j r_j'/c`                         (CG-7)

in the fixed normalized physical coordinates and

`D_[j,m](c)=
 n_(j+m)' Pi_[j,m] n_(j+m)
 +n_(j+m)' Phi_[j,m] G_[j,m](c)^-1
                Phi_[j,m]' n_(j+m)`.                     (CG-8)

CG-8 is exactly LS1 applied once to the superword. It automatically includes
all intermediate quotient information and all B cross-covariance
cancellations. If the block itself has an exact endpoint kernel, its large-c
slope is

`alpha_[j,m]=
 |r_(j+m)' M T_[j,m] r_j|^2/
 [(r_j'Mr_j)(r_(j+m)'Mr_(j+m))]`,                         (CG-9)

where `T_[j,m]=T_(j+m-1)...T_j` only on the zero-action compatibility
trajectory. If any intermediate word forces positive action on the carried
mode, that mode is not in the block kernel; the block Schur information
`j_[j,m]>0` and its contribution to `D_[j,m](c)/c` tends to zero instead.

Thus a sufficient block O2 condition is: for some finite m,

`sup_(same-history exact block-kernel chains) alpha_[j,m] <= 1-delta`
for a `delta>0`.                                          (CG-10)

Together with compact full-baseline `dperp_[j,m]`, CG-10 gives a finite
block ceiling by the same LS5/DC argument. One-word unit transfer is harmless
if it cannot persist through the block.

When every constituent word is exact-kernel and the transported graph remains
exactly on each next graph, write

`T_k r_k=lambda_k r_(k+1)`.                              (CG-11)

Then CG-9 factorizes exactly:

`alpha_[j,m]=prod_(k=j)^(j+m-1) alpha_k`.                 (CG-12)

If a constituent has `alpha_k=1` but a later one has strict loss, the block
contracts. If all `alpha_k=1`, the chain is an exact persistent
compatibility execution. Therefore the finite-block alternative is equivalent
to excluding an infinite same-history unit-transfer chain, not to excluding
unit transfer on every individual word.

The existing EC compactness theorem now applies to the literal graph quantity
CG-5: if no admissible infinite recurring execution satisfies
`T_k r_k=lambda_k r_(k+1)` with unit normalized transfer at every boundary,
then diagonal compactness yields some finite `m` and `delta>0` satisfying
CG-10. Conversely an infinite equality execution defeats every such finite
block.

This is as far as the present assumptions close analytically. The coupled
physical/tuner chronology has NOT yet excluded or constructed the infinite
equality execution. Its exact equations are the already-derived compatibility
zero dynamics/BV system, now with the endpoint quantity fixed by CG-5. The
next decisive calculation is therefore to test global continuation of that
literal graph under the coupled shipping zero dynamics; further arbitrary
soft-eigenvector or separated covariance bounds cannot decide O2.


## Literal compatibility transport inside the full block-PSD certificate

This is the direct substitution of CG/ZG into LS3.  It is the controlling O2
object; no separately maximized d_soft or H appears.

Fix an admissible same-history m-word block
`B=[j,j+m)` and compose the literal Theorem-D factors over the WHOLE block,
with the carried physical/estimator/covariance/tuner/scheduler history.  Write

`Pi_B=Ric_B(0)`, `J_B` for the complete nuisance-shortened root information,
and `Phi_B` for the conditional root-to-terminal map.  Let the literal
endpoint compatibility graph generators be

`r_0=(a_0,0,-A_0a_0)`, `r_1=(a_1,0,-A_1a_1)`.          (BP-1)

Use one fixed physical metric M and normalize
`e_0=r_0/sqrt(r_0'Mr_0)`, `n_1=Mr_1/sqrt(r_1'Mr_1)`.
Equivalently transform root coordinates once by M^(1/2); below the Euclidean
rank-one `e_0e_0'` means that fixed metric normalization, not a wordwise
renormalization.

Define

`w_B=Phi_B' n_1`, `d_B=n_1'Pi_B n_1`,
`G_B(c)=J_B+e_0e_0'/c`.                                  (BP-2)

Then the exact block return is

`D_B(c)=d_B+w_B'G_B(c)^-1 w_B`.                          (BP-3)

All known-root process covariance, diffuse quotient uncertainty, intermediate
S/accelerometer/magnetic information, AW sync, tuner-dependent process
factors and cross-covariances are already in Pi_B,J_B,Phi_B.  In particular
Pi_B is NOT replaced by a scalar background.

The scalar ceiling is exactly the bordered PSD condition

`K_B(c):=
 [[c-d_B, w_B'],
  [w_B,   J_B+e_0e_0'/c]] >=0`.                          (BP-4)

BP-4 is LS3 for the literal superword.  Multiplying the first row/column by
sqrt(c) gives the congruent form

`Khat_B(c)=
 [[c-d_B,       sqrt(c) w_B'],
  [sqrt(c) w_B, c J_B+e_0e_0']] >=0`.                    (BP-5)

This form is useful at rank boundaries because it contains no inverse.

### Short only the true quotient, after inserting the literal graph

Choose a fixed-metric orthonormal root basis `[e_0,E_Q]`.  Write

`J_B=[[j,h'];[h,Q]]`, `w_B=(beta,g)`, `Q>0`.         (BP-6)

Here Q is the actual complete block quotient information.  Do not replace it
by an independent lower floor before taking the Schur complement.  Shorting Q
in BP-4 gives the EXACT 2x2 certificate

`S_B(c)=
 [[c-dperp,             ell],
  [ell, j+1/c]] >=0`,                                    (BP-7)

where

`dperp=d_B+g'Q^-1g`,
`j=j_B:=j-h'Q^-1h>=0`,
`ell=ell_B:=beta-h'Q^-1g`.                              (BP-8)

Thus BP-4 is equivalent to

`c>=dperp`,
`(c-dperp)(j+1/c)-ell^2>=0`.                            (BP-9)

After multiplying by c,

`j c^2+(1-ell^2-j dperp)c-dperp>=0`.                    (BP-10)

This is LS5 with unit fixed-metric root normalization, now derived after the
literal block composition.  The important point is that dperp contains the
complete background Pi_B AND the actual quotient return `g'Q^-1g`; ell
contains the correlated quotient cancellation `h'Q^-1g`.  Neither may be
bounded independently without losing the linked cancellation.

### Exact compatibility face

If the WHOLE block has a nonzero exact root compatibility mode e_0, then

`J_B e_0=0`.                                              (BP-11)

Because J_B is PSD, BP-11 forces `j=0` and `h=0` in BP-6.  Hence

`ell=beta=e_0'Phi_B'n_1
 = r_1'M Phi_B r_0/sqrt[(r_0'Mr_0)(r_1'Mr_1)]`.           (BP-12)

On the zero-action graph trajectory `Phi_B r_0=T_Br_0`, so

`alpha_B:=ell^2
 = |r_1'M T_B r_0|^2/
   [(r_0'Mr_0)(r_1'Mr_1)]`.                              (BP-13)

The complete-background certificate reduces exactly to

`S_B(c)=
 [[c-dperp_B, ell_B],
  [ell_B,     1/c]] >=0`,                                 (BP-14)

or

`c(1-alpha_B)>=dperp_B`.                                 (BP-15)

Therefore:
- if `alpha_B<1`, the exact finite ceiling is
  `c>=dperp_B/(1-alpha_B)`;
- if `alpha_B=1` and `dperp_B>0`, NO finite scalar ceiling exists for that
  block;
- if `alpha_B=1,dperp_B=0`, BP-14 is only semidefinite equality and gives no
  strict contraction.

This conclusion retains the full background.  The obstruction at alpha=1 is
not an artifact of multiplying dbar and Hbar: the positive dperp_B is the
literal known-root/quotient covariance appearing in the same Schur
certificate.

### Intermediate information and the m-word advantage

If the carried endpoint graph direction is charged anywhere inside the block,
then after complete nuisance shorting `j_B>0` unless another exact block
null mode survives.  BP-10 then has positive leading coefficient, so every
fixed block eventually satisfies BP-4 for sufficiently large c.  The positive
root is

`c_*(B)=
 [-(1-ell^2-j dperp)
  +sqrt((1-ell^2-j dperp)^2+4j dperp)]/(2j)`,             (BP-16)

with the j->0 limit given by BP-15 when alpha<1.

Thus the finite-block mechanism is sharper than multiplying boundary alphas:
a mode may have unit endpoint overlap on an early word but acquire positive
actual information later; then the superword has j_B>0 and its large-c slope
is zero.  Only a genuine exact null trajectory through the ENTIRE block lands
on BP-11--BP-15.

For a wholly exact persistent chain, ZG gives
`T_k r_k=lambda_k r_(k+1)`.  Then BP-13 factorizes into the product of the
boundary alphas under the same fixed metric.  If all are one, BP-15 fails
whenever dperp_B>0.  If some boundary is lossy while the block remains exact,
alpha_B<1 and BP-15 supplies the finite ceiling.

### Uniform same-history block theorem

For a fixed m let C_m be the compact class of admissible same-history
m-word blocks with the literal coupled tuner/physical chronology.  Define the
continuous/shorted BP quantities on each closed event stratum.  A uniform
block ceiling exists if and only if the BP-10 positive roots are uniformly
bounded.  A sufficient and, on the exact-kernel face, necessary condition is

`sup_(B in C_m: j_B=0) alpha_B <1`,                      (BP-17)

together with the already required compact finite `dperp_B` and positive
quotient shorting on the retained coordinates.  Then choose

`c_m >= sup_(B in C_m) c_*(B)<infinity`                  (BP-18)

and BP-4 holds for every block.  Strict inequality can be retained by choosing
c_m above the attained supremum and preserving the existing nonlinear/
arithmetic margins.

If BP-17 fails for every finite m because there is an infinite exact
unit-transfer compatibility execution, the scalar block ceiling architecture
cannot close: every prefix lies on BP-15 with alpha=1 and positive background.
If no such infinite execution exists, the EC compactness theorem supplies
some finite m and delta>0 on the exact face; continuity of BP-10 then gives a
finite uniform c_m.  This is the precise bridge from the compatibility
continuation problem to the full-baseline linked Riccati certificate.

No theorem flag is promoted here: the existence/nonexistence of the infinite
exact equality execution remains OPEN.  What is closed is the algebraic
question of how its answer enters O2: through BP-4/BP-10, with the complete
background covariance retained.
