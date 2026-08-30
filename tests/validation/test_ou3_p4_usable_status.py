from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_usable_status as S

class P4UsableStatusTests(unittest.TestCase):
    def test_partial_full_word_progress_does_not_promote_p4(self):
        t={'P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS':'PASS','modes':{m:{'complete_word_translation_margin_lower':1.3108509084723255e-29,'old_single_seed_translation_margin_lower':3.79233618e-35,'margin_widening_factor_lower':345657.88} for m in ('H','A')}}
        d={'modes':{m:{'horizons':{'8.0':{'worst_grid_point':{'translation_complete_word_generalized_margin_design':1.412150733739027e-15}}}} for m in ('H','A')}}
        p={'path_graph_ready':True,'old_worst_corner_has_internal_recurrent_cycle':True}
        out=S.build_from_payloads(t,d,p)
        self.assertEqual(S.validate(out),[])
        self.assertEqual(out['P4_USABLE_CERTIFICATE_STATUS'],'NOT_ESTABLISHED')
        self.assertTrue(out['meaningful_linear_improvement_established'])
        self.assertTrue(out['old_worst_corner_has_internal_recurrent_cycle'])
        for m in ('H','A'):
            self.assertEqual(out['modes'][m]['translation_complete_word_progress'],'PASS')
            self.assertFalse(out['modes'][m]['full_state_complete_word_validated'])
            self.assertFalse(out['modes'][m]['exact_finite_angle_complete_return_map_validated'])

if __name__=='__main__':unittest.main()
