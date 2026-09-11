import unittest
from types import SimpleNamespace

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import physical_word as W

I=Interval.point

def eye(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def prim(k):
    z=(I(0),I(0),I(0))
    return KERNEL.WavePrimitivePayload('g',f'p{k}',f'p{k+1}','live',z,z,z,z,z,z)
def line(k,parent,child,mode='H'):
    n=18 if mode=='H' else 21;p=prim(k)
    common=dict(source_token=child+':e0',predecessor_token=child,mode=mode,sample_index=k,event_ordinal=0,kind='prediction',state=[I(0) for _ in range(n)],P=eye(n),dt_s=I(.005),tau_applied_s=I(1.5),sigma_aw_mps2=I(.4),pseudo_elapsed_s=I(.1),radial_scale=I(1),estimator_source_token=child,estimator_predecessor_token=parent,estimator_generated_coefficients=True,wave_primitive=p)
    if mode=='A':common.update(true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4)
    c=COVER.SourceCoverCell(**common)
    s=SimpleNamespace(prefix_length=k+1,parent_source_cell_id=parent,source_cell_id=child,sample_coordinates=SimpleNamespace(omega_body_corrected=(I(0),I(0),I(0))))
    return SimpleNamespace(mode=mode,cells=(c,),selector=s)
def family(name='BIAS0'):return next(c for c in BIAS.contracts() if c.name==name)

class PhysicalWordTests(unittest.TestCase):
    def test_two_prefix_cocycle_preserves_source_and_bias_ancestry(self):
        d=W.compose_endpoint_lineage([line(0,'root','c1'),line(1,'c1','c2')],bias_contract=family())
        self.assertEqual(d['transitions_composed'],2);self.assertTrue(d['selector_ancestry_continuous']);self.assertTrue(d['primitive_transition_continuity']);self.assertTrue(d['one_generator_over_word']);self.assertTrue(d['one_Live_S_origin_over_word']);self.assertTrue(d['all_physical_blocks_tied_to_global_relation']);self.assertTrue(d['one_bias_parameter_token_over_word']);self.assertFalse(d['complete_word_length'])
    def test_broken_primitive_or_selector_continuity_fails_closed(self):
        a=line(0,'root','c1');b=line(1,'bad','c2')
        with self.assertRaises(ValueError):W.compose_endpoint_lineage([a,b],bias_contract=family())
        b=line(1,'c1','c2');object.__setattr__(b.cells[0].wave_primitive,'primitive_in_id','wrong')
        with self.assertRaises(ValueError):W.compose_endpoint_lineage([a,b],bias_contract=family())
    def test_universal_single_mode_word_induction_closes_but_hybrid_stays_blocked(self):
        d=W.induction_theorem();self.assertEqual(W.validate(d),[])
        self.assertTrue(d['finite_induction_premises_closed']);self.assertTrue(d['H18_complete_word_construction_closed']);self.assertTrue(d['A21_complete_word_construction_closed']);self.assertTrue(d['physical_prediction_forcing_attached']);self.assertTrue(d['physical_S_residual_attached']);self.assertTrue(d['all_bias_families_attached_to_single_mode_word']);self.assertTrue(d['same_history_COMPLETE_BRMM_single_mode_word_closed']);self.assertFalse(d['hybrid_H18_A21_word_closed']);self.assertFalse(d['storage_search_allowed'])

if __name__=='__main__':unittest.main()
