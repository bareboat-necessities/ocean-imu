import unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import physical_defect as D

class PhysicalDefectTests(unittest.TestCase):
    def test_source_uniform_defect_forms_close(self):
        d=D.build();self.assertEqual(D.validate(d),[])
        self.assertTrue(d['source_uniform_defect_graph_closed'])
        self.assertEqual(d['physical_defect_forms']['H18_accelerometer'],'xi=e_ba+eta_acc')
        self.assertEqual(d['physical_defect_forms']['A21_accelerometer'],'xi=eta_acc')
        self.assertEqual(d['physical_defect_forms']['H18_S_zero'],'xi=0')

    def test_held_bias_is_independent_neutral_supply_not_motion_box(self):
        d=D.build()
        self.assertFalse(d['held_bias_bound_uses_startup_membership'])
        self.assertFalse(d['held_bias_bound_uses_covariance_consistency'])
        self.assertFalse(d['held_bias_bound_uses_desired_motion_basin'])
        self.assertGreater(d['held_bias_error_compactness_upper_mps2'],0.4)
        self.assertEqual(set(d['true_bias_norm_upper_by_family_mps2']),{'BIAS0','BIAS1','BIAS2'})

    def test_rank_three_energy_route_retained(self):
        d=D.build()
        self.assertEqual(d['positive_defect_port_dimension'],3)
        self.assertTrue(d['H18_acc_C_inverse_dominated_by_Racc_inverse'])
        self.assertFalse(d['independent_C_inverse_box_used'])
        self.assertFalse(d['source_uniform_complete_word_dissipation_closed_here'])

if __name__=='__main__':unittest.main()
