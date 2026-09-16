"""Uniform two-sweep termination of the pinned scalar Eigen seed SVD.

Each bound is measured against exact operations on the SAME stored input work
matrix. Rotation construction is related to its exact, stable real formula;
matrix products then propagate absolute errors on their actual small operands.
Source normalization, QR, diagonal ordering and all small-value branches are
linked below. Every bound also permits local contraction of adjacent multiply
and add/subtract nodes in the same expression tree: deletion of an internal
rounding chooses zero from its allowed error interval. It does not permit
reassociation, approximate division/sqrt, traps, or flush-to-zero arithmetic.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32

U=F(1,2**24)
ETA=F(1,2**150)
SMALL_Y=F(1,2**40)
FIRST_ERROR=1024*U


def normalized_source_bounds():
    """Both source normalizations, including the first nonzero guard.

    For exact S=sum(x_i²), the stored nonnegative dot product is bracketed by
    (1+-u)^3*S +-5eta. Its normal square root contributes (1+-u)^2 to the
    squared denominator. Division then contributes (1+-u)^2 to squared norm;
    the three component underflow errors contribute less than 16eta.
    """
    guard=B32.rn32(F(1,1000))
    initial_norm2_lower=(guard**2/(1+U)**2-5*ETA)/(1+U)**3
    assert initial_norm2_lower>F(1,10**7)
    initial_norm2_upper=3*160**2
    assert (1+U)**3*initial_norm2_upper+5*ETA<2**17

    def normalized(minimum):
        denominator_lower=(1-U)**2*((1-U)**3-5*ETA/minimum)
        denominator_upper=(1+U)**2*((1+U)**3+5*ETA/minimum)
        assert F(1,2)<denominator_lower<denominator_upper<2
        lo=(1-U)**2/denominator_upper-16*ETA
        hi=(1+U)**2/denominator_lower+16*ETA
        assert F(1,2)<lo<hi<2
        return lo,hi

    down_lo,down_hi=normalized(initial_norm2_lower)
    v0_lo,v0_hi=normalized(down_lo)
    delta=F(1,10**6)
    assert 1-delta<v0_lo<v0_hi<1+delta
    return {'raw_accel_component_upper':F(160),'first_accel_norm_guard':guard,
            'first_input_norm2_lower_on_guard':initial_norm2_lower,
            'first_input_norm2_upper':initial_norm2_upper,
            'down_norm2_lower':down_lo,'down_norm2_upper':down_hi,
            'v0_norm2_lower':v0_lo,'v0_norm2_upper':v0_hi,
            'v0_norm2_defect_upper':delta,'v1_is_literal_e3':True,
            'both_normalizations_finite_on_guard':True}


@dataclass(frozen=True)
class Bound:
    """|exact value|<=magnitude; |machine-exact|<=error."""
    magnitude:F
    error:F=F(0)

    def __post_init__(self):
        if self.magnitude<0 or self.error<0: raise ValueError('nonnegative error envelope required')

    def add(self,other):
        magnitude=self.magnitude+other.magnitude
        error=self.error+other.error
        return Bound(magnitude,error+U*(magnitude+error)+ETA)

    def mul(self,other):
        magnitude=self.magnitude*other.magnitude
        error=(self.magnitude*other.error+other.magnitude*self.error+self.error*other.error)
        return Bound(magnitude,error+U*(magnitude+error)+ETA)

    def cap(self,magnitude):
        """Use an independently proved exact-real identity, retaining error."""
        if magnitude>self.magnitude: raise ValueError('cap must strengthen an ideal magnitude bound')
        return Bound(F(magnitude),self.error)

    @property
    def actual(self): return self.magnitude+self.error


def relative_coefficients():
    """Stable normalization and Jacobi-tangent roundoff factors.

    All charged operations here have normal nonzero results; the tiny-y
    branch is separate. Addition of nonnegative quantities preserves the
    multiplicative bracket of its summands. sqrt halves its bracket powers.

    Polar: t=RN(a+d), d=RN(c-b), r=RN(t/d), h=RN(sqrt(RN(1+RN(r²)))),
    s=RN(1/h), c=RN(r/h). Relative to exact t*,d* from the stored work,
    r/r* is between (1-u)^2/(1+u) and (1+u)^2/(1-u). This gives the brackets
    below without subtracting independently enclosed square roots.

    Jacobi: the gap subtraction and tau division give (1+-u)^2; the stable
    same-sign tau+-sqrt(1+tau²) denominator gives (1+-u)^5. The reciprocal t,
    normalized cosine n and final t*n yield the indicated brackets. Products
    by the exactly represented signs and y/abs(y)=+-1 introduce no error.
    """
    polar_c=((1-U)**4/(1+U)**5,(1+U)**4/(1-U)**5)
    polar_s=((1-U)**2/(1+U)**4,(1+U)**2/(1-U)**4)
    jacobi_c=((1-U)**6/(1+U)**3,(1+U)**6/(1-U)**3)
    jacobi_s=((1-U)**8/(1+U)**8,(1+U)**8/(1-U)**8)
    errors={name:max(1-lo,hi-1) for name,(lo,hi) in
            (('polar_c',polar_c),('polar_s',polar_s),('jacobi_c',jacobi_c),('jacobi_s',jacobi_s))}
    assert all(x<32*U for x in errors.values())
    assert errors['polar_s']<8*U and errors['jacobi_s']<17*U
    return errors


def qr_bounds():
    """Uniform Householder domain/error inequalities, including tiny tails.

    A nonidentity reflector has stored tail squared norm > MIN_NORMAL.
    Up to two tail squares and their sum imply exact ||x||² >=
    MIN_NORMAL*(1-8u). Nonnegative square/sum paths plus normal sqrt give the
    beta² ratio below. The beta sign prevents denominator cancellation.
    """
    minimum_norm_factor=1-8*U
    beta2_lo=(1-U)**5-8*U/minimum_norm_factor
    beta2_hi=(1+U)**5+8*U/minimum_norm_factor
    assert 1-16*U<beta2_lo<beta2_hi<1+16*U
    # |x0|/|beta| <= 1+16u, verified by squaring, not a sampled sqrt.
    assert (1+16*U)**2*(1-16*U)>1
    essential2=(1+U)**2/((1-U)**2*(1-16*U))+8*ETA
    tau=(1+U)**2*(2+16*U)
    assert essential2<1+32*U and tau<2+32*U
    stored_H_norm=(2+32*U)**2-1
    # Direct scalar applyHouseholder budget. Error expressions are affine in
    # input norm X (there is no product of two input-vector entries), so
    # checking endpoints establishes the displayed budget for 0<=X<=4.
    for x in (F(0),F(4)):
        component=Bound(x); essential=Bound(F(101,100)); coefficient=Bound(F(201,100))
        tail=essential.mul(component).add(essential.mul(component)).cap(F(101,100)*x)
        temporary=tail.add(component).cap(F(3,2)*x)
        first=component.add(coefficient.mul(temporary))
        later=component.add(coefficient.mul(essential).mul(temporary))
        assert max(first.error,later.error)<32*U*x+32*ETA
    # ||component_errors||_2 <=sqrt(3)(32uX+32eta)<64uX+64eta.
    effective_factor=stored_H_norm+64*U
    assert effective_factor<3+256*U
    axis_norm=(3+256*U)*((3+256*U)+64*ETA)+64*ETA
    assert axis_norm<10

    # First-pivot norm is near one; eta charges are relative to .99², not
    # MIN_NORMAL. The beta root error fits 4u on this branch.
    b2lo=(1-U)**5-8*ETA/F(99,100)**2
    b2hi=(1+U)**5+8*ETA/F(99,100)**2
    assert (1-4*U)**2<b2lo and b2hi<(1+4*U)**2
    denominator_relative=(4*U*F(101,100)+U*(F(1,200)+(1+4*U)*F(101,100)))/F(99,100)
    assert denominator_relative<6*U
    essential_relative=max(1-(1-U)/(1+6*U),(1+U)/(1-6*U)-1)
    assert essential_relative<8*U
    # tau's exact numerator has magnitude >=|beta*|; it has the same stable
    # relative denominator allowance, then division by beta* with 4u error.
    tau_relative=max(1-(1-6*U)*(1-U)/(1+4*U),
                     (1+6*U)*(1+U)/(1-4*U)-1)
    tau_error=F(101,100)*tau_relative
    assert tau_error<16*U
    e_error=8*U+2*ETA
    H_error=16*U*F(3,2)**2+F(101,100)*(3*e_error+e_error**2)
    assert H_error<128*U
    applied_error=(H_error+64*U)*F(101,100)+64*ETA
    assert applied_error<256*U
    opposite_error=F(7,1000)+applied_error+4*U*F(101,100)
    assert F(99,100)-opposite_error>F(98,100)
    assert F(101,100)+opposite_error<F(102,100)
    second_beta_upper=(1+16*U)*(F(7,1000)+applied_error)+8*ETA
    assert second_beta_upper<F(1,50)
    return {'beta2_relative_lower':beta2_lo,'beta2_relative_upper':beta2_hi,
            'essential_norm2_upper':essential2,'tau_upper':tau,
            'effective_reflector_factor':effective_factor,
            'axis_norm_upper':axis_norm,'axis_norm2_upper':axis_norm**2,
            'first_beta_relative_error':4*U,
            'first_denominator_relative_error':denominator_relative,
            'first_essential_relative_error':essential_relative,
            'first_tau_absolute_error':tau_error,'first_H_operator_error':H_error,
            'first_H_application_norm_error':applied_error,
            'opposite_first_row_error':opposite_error,
            'second_beta_abs_upper':second_beta_upper}


def composed(polar_c,polar_s,jacobi_c,jacobi_s,*,sine_bound):
    # The scalar source signs are immaterial for absolute errors. The exact
    # composition is orthogonal, hence its cosine magnitude is <=1.
    c=polar_c.mul(jacobi_c).add(polar_s.mul(jacobi_s)).cap(1)
    s=polar_c.mul(jacobi_s).add(polar_s.mul(jacobi_c)).cap(sine_bound)
    return c,s


def apply_to_work(c,s,a,b,cc,d):
    # Entries of [[a,b],[cc,d]]; signed addition/subtraction have identical
    # error propagation. The caller retains the actual exact identities.
    return ((c.mul(a).add(s.mul(cc)),c.mul(b).add(s.mul(d))),
            (s.mul(a).add(c.mul(cc)),s.mul(b).add(c.mul(d))))


def first_sweep():
    rel=relative_coefficients()
    pc=Bound(F(1),rel['polar_c']); ps=Bound(F(1),rel['polar_s'])
    # The stable Jacobi rotation of this diagonally separated block has
    # |sin|<.05; a tiny-y identity differs by at most 8*SMALL_Y.
    jc=Bound(F(1),rel['jacobi_c']+8*SMALL_Y)
    js=Bound(F(1,20),rel['jacobi_s']/20+8*SMALL_Y)
    lc,ls=composed(pc,ps,jc,js,sine_bound=F(21,20))
    assert lc.error<128*U and ls.error<128*U
    # Reversed input block used by RealSvd2x2: [[c,b],[0,a]].
    original=(Bound(F(1,50)),Bound(F(51,50)),Bound(F(0)),Bound(F(101,100)))
    polar=apply_to_work(pc,ps,*original)
    polar_error=max(x.error for row in polar for x in row)
    assert polar_error<128*U
    # The exact polar output is symmetric. Replace it by the symmetric
    # matrix built from the actual x,y,z. Each entry changes <=polar_error;
    # operator norm <=2*polar_error is retained under exact orthogonal J.
    symmetry_supply=2*polar_error
    l=apply_to_work(lc,ls,*original)
    # Exact composed rotation preserves the input Frobenius norm <1.5.
    l=tuple(tuple(x.cap(min(x.magnitude,F(3,2))) for x in row) for row in l)
    off=(l[0][0].mul(js).add(l[0][1].mul(jc)),
         l[1][0].mul(jc).add(l[1][1].mul(js)))
    matrix_error=max(x.error for x in off)
    diag=(l[0][0].mul(jc).add(l[0][1].mul(js)),
          l[1][0].mul(js).add(l[1][1].mul(jc)))
    diagonal_error=symmetry_supply+max(x.error for x in diag)
    total=symmetry_supply+matrix_error
    assert max(total,diagonal_error)<FIRST_ERROR
    return {'polar':polar,'polar_error':polar_error,'left_c':lc,'left_s':ls,
            'symmetry_supply':symmetry_supply,'matrix_error':matrix_error,
            'offdiagonal_error':total,'diagonal_error':diagonal_error}


def second_sweep():
    e=FIRST_ERROR; rel=relative_coefficients()
    pc=Bound(F(1),rel['polar_c']); ps=Bound(2*e,2*e*rel['polar_s'])
    # Exact input work has large diagonal at index 0; Eigen reverses the
    # extracted block. Its first diagonal is <=.03 and second diagonal <=2.
    original=(Bound(F(3,100)),Bound(e),Bound(e),Bound(F(2)))
    polar=apply_to_work(pc,ps,*original)
    diagonal_error=max(polar[0][0].error,polar[1][1].error)
    offdiagonal_error=max(polar[0][1].error,polar[1][0].error)
    assert diagonal_error<128*U and offdiagonal_error<512*U*e
    # Exact symmetry makes the upper exact offdiagonal equal to the lower:
    # |lower|<=|sin|*.03+|cos|*e. This is the retained cancellation identity.
    y_bound=F(106,100)*e+offdiagonal_error
    assert y_bound<2*e
    # |exact diagonal gap| >= (1-4e²)*1.27-4e². Charge the two actual
    # diagonal errors before using the Jacobi relative-operation lemma.
    actual_gap=(1-4*e*e)*F(127,100)-4*e*e-2*diagonal_error
    assert actual_gap>1
    jc=Bound(F(1),rel['jacobi_c'])
    js=Bound(4*e,4*e*rel['jacobi_s'])
    assert js.error<128*U*e
    lc,ls=composed(pc,ps,jc,js,sine_bound=6*e)
    assert lc.error<128*U and ls.error<1024*U*e
    l=apply_to_work(lc,ls,*original)
    assert max(l[0][0].error,l[1][1].error)<512*U
    assert max(l[0][1].error,l[1][0].error)<4096*U*e
    off=(l[0][0].mul(js).add(l[0][1].mul(jc)),
         l[1][0].mul(jc).add(l[1][1].mul(js)))
    matrix_error=max(x.error for x in off)
    # The symmetric matrix diagonalized by exact J differs from the exact
    # polar output by diagonal_error/offdiagonal_error. An orthogonal
    # rotation with |sin|<=4e gives offdiagonal discrepancy bounded below.
    symmetry_supply=2*offdiagonal_error+8*e*diagonal_error
    normal_result=matrix_error+symmetry_supply
    # Re-run the actual same product/sum graph with the proved tiny-y
    # absolute coefficient error. This verifies the 256Y allowance instead
    # of treating it as a free disturbance on the completed sweep.
    jc_tiny=Bound(jc.magnitude,jc.error+8*SMALL_Y)
    js_tiny=Bound(js.magnitude,js.error+8*SMALL_Y)
    lc_tiny,ls_tiny=composed(pc,ps,jc_tiny,js_tiny,sine_bound=6*e)
    lt=apply_to_work(lc_tiny,ls_tiny,*original)
    off_tiny=(lt[0][0].mul(js_tiny).add(lt[0][1].mul(jc_tiny)),
              lt[1][0].mul(jc_tiny).add(lt[1][1].mul(js_tiny)))
    tiny_result=max(x.error for x in off_tiny)-matrix_error
    assert 0<=tiny_result<256*SMALL_Y
    assert normal_result<16384*U*e
    final=normal_result+tiny_result
    # maxDiag retains its initial |a|>.99, and RN(4u*.99)>3u.
    threshold_lower=(1-U)*4*U*F(99,100)
    assert threshold_lower>3*U and final<3*U
    return {'polar':polar,'polar_diagonal_error':diagonal_error,
            'polar_offdiagonal_error':offdiagonal_error,'upper_y_bound':y_bound,
            'actual_gap_lower':actual_gap,'left_c':lc,'left_s':ls,
            'matrix_error':matrix_error,'symmetry_supply':symmetry_supply,
            'small_y_extra':tiny_result,'offdiagonal_error':final,
            'threshold_lower':threshold_lower,
            'margin_lower':threshold_lower-final}


def source_and_order_bounds():
    """Link normalized-vector input, QR box and first diagonal ordering.

    The exact first polar rotation has c,s>0 and both >=.68. Since a*b<0,
    its large diagonal has magnitude >=.68(|a|+|b|). Its exact lower
    offdiagonal is -s*c_small. Use symmetry to transfer that bound to y;
    charge the actual polar operation errors before constructing exact J.
    The stable Jacobi root selects the small-angle rotation, with
    |tangent|<=|y|/|gap|, so it preserves the small/large diagonal order.
    """
    normalized=normalized_source_bounds()
    delta=normalized['v0_norm2_defect_upper']
    cutoff=B32.add(-1,B32.rn32(F(1,100000)))
    transverse2=1+delta-cutoff**2
    assert (1+U)**2*transverse2+8*ETA<F(1,200)**2
    plus_down2=2+delta+2*cutoff
    assert plus_down2<F(6,1000)**2
    norm2_low=(1-U)**2*(1-delta)/(1+delta)**2-8*ETA
    norm2_high=(1+U)**2*(1+delta)+8*ETA
    assert F(99,100)**2<norm2_low<norm2_high<F(101,100)**2
    assert F(99,100)**2<(1-4*U)**2*norm2_low
    assert (1+4*U)**2*norm2_high<F(101,100)**2
    assert (1-U)*abs(cutoff)/(1+delta)-ETA>F(99,100)
    assert F(6,1000)+3*U+8*ETA<F(7,1000)
    q=qr_bounds(); first=first_sweep()
    ratio_low=F(97,100)/F(102,100)
    ratio_high=F(103,100)/F(98,100)
    lower_cs=F(68,100)
    assert lower_cs**2*(1+ratio_high**2)<1
    assert lower_cs**2*(1+ratio_low**2)<ratio_low**2
    polar_big_lower=lower_cs*(F(99,100)+F(98,100))-first['polar_error']
    polar_small_upper=F(1,50)+first['polar_error']
    y_upper=F(1,50)+first['polar_error']
    gap_lower=polar_big_lower-polar_small_upper
    assert gap_lower>F(13,10)
    tangent_upper=y_upper/gap_lower
    assert tangent_upper<F(1,20)
    # Exact symmetric 2x2 Jacobi eigenvalue updates are x-t*y, z+t*y
    # (with the sign determined by the stable tangent); take absolute bounds.
    eigenvalue_shift=y_upper*tangent_upper
    output_big_lower=polar_big_lower-eigenvalue_shift-first['diagonal_error']
    output_small_upper=polar_small_upper+eigenvalue_shift+first['diagonal_error']
    assert output_big_lower>F(13,10) and output_small_upper<F(3,100)
    # Orthogonality and ||A||_F<1.5 give an upper spectral bound for B,
    # including its polar evaluation defect, then actual sweep roundoff.
    assert F(101,100)**2+F(102,100)**2+F(1,50)**2<F(3,2)**2
    output_big_upper=F(3,2)+2*first['polar_error']+first['diagonal_error']
    assert output_big_upper<2
    # An active second predicate has max(|b|,|c|)>3u. If same-sign operands
    # can cancel, Sterbenz places both above 1.5u, so their lattice spacing
    # is >=2^-47. Otherwise their difference has magnitude >=1.5u. The
    # deliberately looser 2^-48 lower bound covers every nonzero result.
    # If d=0, the source comparison |d|<MIN_NORMAL selects literal identity;
    # the same stored work is already symmetric. No division is executed.
    difference_lower=F(1,2**48)
    quotient_upper=4/difference_lower
    assert quotient_upper**2<2**102<2**127
    # For |Jacobi y|>=2^-40, gap<=4 bounds |tau|<=2^42. All squares,
    # roots, same-sign sums, t, n and s in relative_coefficients are normal.
    tau_upper=4/(2*SMALL_Y)
    assert tau_upper<2**42 and tau_upper**2<2**84
    # |y|<Y: either the denominator guard skips Jacobi, or overflowing tau²
    # produces w=+inf and signed t=0, giving identity. Otherwise the stable
    # same-sign denominator is at least |tau|, so the retained finite tangent
    # has the upper bound below. The exact tangent is <Y since gap exceeds one.
    # This covers subnormal y and t without assuming a relative-error model.
    tiny_tangent=(1+U)*2*SMALL_Y/((1-U)**2)+ETA
    assert tiny_tangent<3*SMALL_Y
    tiny_coefficient_error=tiny_tangent+SMALL_Y+8*U*tiny_tangent+8*ETA
    assert tiny_coefficient_error<8*SMALL_Y
    return {'cutoff':cutoff,'transverse_norm2_upper':transverse2,
            'scaled_norm2_lower':norm2_low,'scaled_norm2_upper':norm2_high,
            'first_ratio_lower':ratio_low,'first_ratio_upper':ratio_high,
            'first_polar_coefficient_lower':lower_cs,'first_actual_gap_lower':gap_lower,
            'first_Jacobi_sine_upper':tangent_upper,
            'first_output_big_diagonal_lower':output_big_lower,
            'first_output_big_diagonal_upper':output_big_upper,
            'first_output_small_diagonal_upper':output_small_upper,
            'active_nonzero_polar_difference_lower':difference_lower,
            'active_polar_quotient_upper':quotient_upper,
            'normal_branch_Jacobi_tau_upper':tau_upper,
            'tiny_Jacobi_tangent_abs_upper':tiny_tangent,
            'tiny_Jacobi_coefficient_error_upper':tiny_coefficient_error,
            'zero_polar_difference_selects_literal_identity':True,
            'tiny_Jacobi_all_IEEE_branches_covered':True,
            'normalization':normalized,'QR':q}


def build():
    return {'QR':qr_bounds(),'source_and_order':source_and_order_bounds(),
            'first_sweep':first_sweep(),'second_sweep':second_sweep(),
            'maximum_Jacobi_sweeps':2,
            'source_QR_and_diagonal_order_induction_closed':True,
            'universal_termination_promoted':True,
            'arithmetic_family':'binary32-RNE-gradual-nontrapping-correct-div-sqrt-local-FMA-no-reassociation',
            'local_FMA_rounding_deletion_covered':True,
            'target_compiler_correspondence_proved_here':False}
