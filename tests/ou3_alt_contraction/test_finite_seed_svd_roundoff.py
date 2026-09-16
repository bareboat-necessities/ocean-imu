from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_seed_svd_roundoff as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_seed_eigen_svd as S
from test_finite_seed_eigen_svd import inputs


class Tests(unittest.TestCase):
    def test_relative_factors_and_operation_bounds_are_exact_rationals(self):
        r=X.relative_coefficients()
        self.assertLess(r['polar_s'],8*X.U)
        self.assertLess(r['jacobi_s'],17*X.U)
        for a,b in ((F(1)+X.U,F(-1)+X.U),
                    (F(1,2**100),F(1,2**49)),
                    (F(1,2**149),F(-1,2**149))):
            aa,bb=B.rn32(a),B.rn32(b)
            ea=X.Bound(abs(a),abs(aa-a)); eb=X.Bound(abs(b),abs(bb-b))
            self.assertLessEqual(abs(B.add(aa,bb)-(a+b)),ea.add(eb).error)
            self.assertLessEqual(abs(B.mul(aa,bb)-a*b),ea.mul(eb).error)

    def test_QR_and_diagonal_induction_bounds_keep_room_for_rounding(self):
        r=X.source_and_order_bounds()
        source=r['normalization']
        self.assertGreater(source['v0_norm2_lower'],1-F(1,10**6))
        self.assertLess(source['v0_norm2_upper'],1+F(1,10**6))
        self.assertTrue(r['zero_polar_difference_selects_literal_identity'])
        self.assertTrue(r['tiny_Jacobi_all_IEEE_branches_covered'])
        self.assertLess(r['tiny_Jacobi_coefficient_error_upper'],8*X.SMALL_Y)
        self.assertLess(r['QR']['axis_norm2_upper'],100)
        self.assertGreater(r['first_actual_gap_lower'],F(13,10))
        self.assertLess(r['first_Jacobi_sine_upper'],F(1,20))
        self.assertGreater(r['first_output_big_diagonal_lower'],F(13,10))
        self.assertLess(r['first_output_small_diagonal_upper'],F(3,100))
        self.assertLess(r['first_output_big_diagonal_upper'],2)

    def test_second_bound_is_below_the_actual_retained_threshold(self):
        r=X.second_sweep()
        self.assertLess(r['polar_offdiagonal_error'],512*X.U*X.FIRST_ERROR)
        self.assertLess(r['upper_y_bound'],2*X.FIRST_ERROR)
        self.assertGreater(r['actual_gap_lower'],1)
        self.assertLess(r['small_y_extra'],256*X.SMALL_Y)
        self.assertLess(r['offdiagonal_error'],3*X.U)
        self.assertGreater(r['threshold_lower'],3*X.U)
        self.assertGreater(r['margin_lower'],0)

    def test_native_input_cases_satisfy_the_uniform_induction_boxes(self):
        # Sanity check only; the quantified inequalities are the rational
        # derivation and its reviewed source-domain premises.
        for v0 in inputs():
            r=S.solve_bits(v0,(0,0,S.ONE))
            self.assertLess(sum(S.M.value(b)**2 for b in r.axis),100)
            m=r.work_history[1]
            self.assertLess(max(abs(S.M.value(m[0][1])),abs(S.M.value(m[1][0]))),X.FIRST_ERROR)
            self.assertGreater(abs(S.M.value(m[0][0])),F(13,10))
            self.assertLess(abs(S.M.value(m[1][1])),F(3,100))


if __name__=='__main__': unittest.main()
