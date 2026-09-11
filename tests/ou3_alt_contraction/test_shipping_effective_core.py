import unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import shipping_effective_core as C

class ShippingEffectiveCoreTests(unittest.TestCase):
    def test_shipping_hold_semantics_close_core_form(self):
        d=C.build(); self.assertEqual(C.validate(d),[])
        self.assertTrue(d['all_required_shipping_parity_closed'])
        self.assertTrue(d['source_uniform_effective_core_form_closed'])
        self.assertTrue(d['source_uniform_C_positive_definite_from_R_positive_and_Pba_PSD'])

    def test_exact_H18_block_reduction(self):
        d=C.build(); b=d['exact_block_regression']
        self.assertTrue(b['acc_core_equals_R_plus_Pba'])
        self.assertTrue(b['mag_core_equals_R'])
        self.assertTrue(b['Szero_core_equals_R'])
        self.assertEqual(d['H18_acc_effective_core'],'C=R_acc+P_ba_ba')

    def test_core_positivity_does_not_import_unrelated_blockers(self):
        d=C.build()
        self.assertFalse(d['PE_vector_domain_lemma_needed_for_core_positivity'])
        self.assertFalse(d['Mahony_startup_needed_for_core_positivity'])
        self.assertFalse(d['covariance_tube_upper_bound_needed_for_core_positivity'])
        self.assertFalse(d['source_uniform_complete_word_dissipation_closed_here'])

if __name__=='__main__': unittest.main()
