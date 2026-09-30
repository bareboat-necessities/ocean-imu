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
