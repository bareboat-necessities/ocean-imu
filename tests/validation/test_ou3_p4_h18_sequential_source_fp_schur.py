#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_sequential_source_fp_schur as S
class H18SequentialSourceFpSchurTest(unittest.TestCase):
    def test_contract(self):
        d=S.build();self.assertEqual(S.validate(d),[])
        self.assertEqual(d['persistent_dimension'],22)
        self.assertTrue(d['same_history_a_endpoint_carried_across_predictions'])
        self.assertTrue(d['local_source_and_fp_inputs_eliminated_immediately'])
        self.assertTrue(d['outward_interval_Schur_complement_available'])
        self.assertFalse(d['dense_whole_word_input_master_required'])
        self.assertFalse(d['P4_PASS'])
    def test_prediction_transition_has_fixed_local_dimension(self):
        A=S.eye(18);T=S.prediction_transition(A,S.I(2.0),S.I(.005),1e-6)
        self.assertEqual(S.shape(T),(22,52))
        sectors=S.prediction_source_sectors(52)
        self.assertEqual([n for n,_ in sectors],['physical_a0_ball','physical_a1_ball','physical_J012_moment_iqc'])
        self.assertTrue(all(S.shape(P)==(52,52) for _,P in sectors))
    def test_ordinary_event_does_not_add_physical_source_block(self):
        T=S.ordinary_transition(S.eye(18),1e-6);self.assertEqual(S.shape(T),(22,40))
    def test_schur_smoke_is_fixed_persistent_shape(self):
        s=S._smoke();self.assertEqual(tuple(s['required_state_shape']),(22,22));self.assertEqual(tuple(s['Kuu_shape']),(18,18));self.assertTrue(s['finite'])
if __name__=='__main__':unittest.main()
