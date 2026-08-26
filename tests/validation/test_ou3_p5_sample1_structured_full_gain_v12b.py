from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v12b as V12B


class StructuredFullGainV12BTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V12B.build(source_pieces=4,source_cell_index=0,
                         p_pieces=2,tangent_pieces=2,axial_pieces=2)

    def test_exact_psd_resolvent_has_no_small_denominator_subtraction(self):
        r=V12B._gain_perturbation_psd_floor(
            dP=1e-6,dH=2e-4,pnorm=0.03,hnorm=20.0,
            ktheta_norm=0.1,r_floor=8e-4)
        self.assertGreater(r["sample1_gain_operator_perturbation_upper"],0.0)
        self.assertAlmostEqual(r["actual_innovation_inverse_operator_upper"],1.0/8e-4,delta=1e-9)
        self.assertGreater(r["sample1_innovation_perturbation_upper"],0.0)

    def test_psd_shipping_and_jacobian_semantics_are_explicit(self):
        for k in (
            "actual_shipping_covariance_PSD_used",
            "actual_innovation_noise_floor_Ra_used_without_subtraction",
            "exact_gain_resolvent_identity_used",
            "sample1_force_Jacobian_perturbation_included",
            "V11_PSD_and_S_perturbation_magnitudes_retained",
        ):
            self.assertIs(self.d[k],True)
        self.assertIs(self.d["generic_interval_innovation_LDLT_floor_used"],False)
        self.assertIs(self.d["perturbation_subtracted_from_measurement_noise_floor"],False)

    def test_result_is_fail_closed(self):
        st=self.d["P5_SAMPLE1_PSD_S_PSD_RESOLVENT_V12B"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        self.assertGreater(self.d["evaluated_joint_cells"],0)
        self.assertTrue(math.isfinite(float(self.d["max_V12B_correction_norm_upper_rad"])))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")
        if st=="PASS":
            self.assertLess(float(self.d["max_V12B_correction_norm_upper_rad"]),9.0)

    def test_no_promotion_or_shipping_limit_change(self):
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]),6.0)
        for k in (
            "source_replay_used","filter_changed","deployed_correction_limit_increased",
            "complete_sample1_branch_closed_here","signed_cayley_q8_composed_here",
            "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)


if __name__=="__main__": unittest.main()
