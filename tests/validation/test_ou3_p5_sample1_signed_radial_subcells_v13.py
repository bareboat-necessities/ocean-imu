from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_signed_radial_subcells_v13 as V13


class Sample1SignedRadialSubcellsV13Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Fast fail-closed fixture.  The authoritative V12/V13 family remains
        # 24^3 with signed correction subdivision in the focused producer.
        cls.d=V13.build(source_pieces=4,source_cell_index=0,
                        p_pieces=2,tangent_pieces=2,axial_pieces=2,
                        residual_x_pieces=4,parallel_pieces=4)

    def test_signed_gain_components_are_finite(self):
        I=V13.FULL.I
        perp,par,detail=V13._signed_gain_components(
            a=I(0.02),Y=I(0.03),c0=I(-0.004),alpha=I(0.99),qaw=I(0.001),
            b=I(0.04),bz=I(0.02),det_first=I(0.0005),d=I(0.5),
            fy=I(-1.0),fz=I(9.0),r=I(0.001))
        for x in (*perp,*par):
            self.assertTrue(math.isfinite(x.lo))
            self.assertTrue(math.isfinite(x.hi))
        self.assertGreater(detail["scalar_Sx"][0],0.0)
        self.assertGreater(detail["two_by_two_det"][0],0.0)

    def test_coarse_fixture_is_fail_closed_until_v12(self):
        self.assertIs(self.d["V12_actual_innovation_PSD_S_prerequisite_required"],True)
        if not self.d["V12_prerequisite_passed"]:
            self.assertEqual(self.d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13"],"NOT_ESTABLISHED")
            self.assertEqual(self.d["evaluated_signed_subcells"],0)
            self.assertTrue(any("V12 prerequisite did not pass" in x for x in self.d["failures"]))

    def test_validate_accepts_fail_closed_or_closed_result(self):
        f=V13.validate(self.d)
        if self.d["V12_prerequisite_passed"]:
            self.assertEqual(f,[])
        else:
            # The producer deliberately reports the unmet V12 prerequisite as
            # validation evidence instead of pretending the signed cover ran.
            self.assertTrue(any("V12 prerequisite did not pass" in x for x in f))

    def test_no_promotion_or_shipping_limit_change(self):
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]),6.0)
        for k in (
            "source_replay_used","filter_changed","deployed_correction_limit_increased",
            "signed_cayley_q8_composed_here","complete_sample1_branch_closed_here",
            "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)

    def test_radial_contract_when_prerequisite_closes(self):
        if not self.d["V12_prerequisite_passed"]:
            self.skipTest("coarse V12 prerequisite intentionally not closed")
        self.assertIs(self.d["signed_one_plus_two_correction_cover_used"],True)
        self.assertIs(self.d["residual_ball_not_double_counted_across_blocks"],True)
        self.assertIs(self.d["V12_correction_perturbation_ball_retained"],True)
        st=self.d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        if st=="PASS":
            self.assertEqual(self.d["unclosed_radial_subcells"],0)
            self.assertLessEqual(self.d["max_radial_upper"],9.0)
            if self.d["above_6rad_subcells"]:
                self.assertGreater(self.d["minimum_radial_lower_above_6"],0.0)


if __name__=="__main__": unittest.main()
