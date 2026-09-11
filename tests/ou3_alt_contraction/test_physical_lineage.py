import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'tools/stability')]
from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
from tools.stability.ou3_alt_contraction import physical_lineage as L

I=Interval.point

def eye(n):
    return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]

def primitive(S=(2.0,-1.0,.5)):
    z=(I(0),I(0),I(0))
    return KERNEL.WavePrimitivePayload(
        'generator-1','phi-in-7','phi-out-8','live-origin-3',
        z,z,tuple(I(x) for x in S),z,z,z)

def cell(mode,kind,ordinal,*,S=(2.0,-1.0,.5)):
    n=18 if mode=='H' else 21
    kw=dict(
        source_token=f'est:e{ordinal}', predecessor_token='est', mode=mode,
        sample_index=7,event_ordinal=ordinal,kind=kind,state=[I(0) for _ in range(n)],
        P=eye(n),dt_s=I(.005),tau_applied_s=I(1.7),sigma_aw_mps2=I(.4),
        pseudo_elapsed_s=I(.1),radial_scale=I(1),estimator_source_token='est',
        estimator_predecessor_token='root',estimator_generated_coefficients=True,
        wave_primitive=primitive(S))
    if kind=='S_zero':
        kw.update(R=eye(3),R_provenance=EVENTS.ACTUAL_RS_PROVENANCE)
    if kind=='accelerometer':
        kw.update(R=eye(3),f_hat=[I(0),I(0),I(-9.80665)],R_hat=eye(3))
    if kind=='magnetometer':
        kw.update(R=eye(3),m_body=[I(20),I(0),I(40)])
    if mode=='A':
        kw.update(true_bias=[I(.03),I(-.02),I(.01)],bias_projection_limit=.4)
    return COVER.SourceCoverCell(**kw)

def lineage(mode,cells):
    sample=SimpleNamespace(omega_body_corrected=(I(.01),I(-.02),I(.005)))
    selector=SimpleNamespace(sample_coordinates=sample)
    return SimpleNamespace(mode=mode,cells=tuple(cells),selector=selector)

class PhysicalLineageTests(unittest.TestCase):
    def test_H18_q15_and_physical_S_are_one_literal_cocycle(self):
        p=cell('H','prediction',0)
        s=cell('H','S_zero',1)
        lin=lineage('H',[p,s])
        d=L.compose_attached_sample(
            lin,phi_true=I(.99999),
            held_bias_error=[Interval(-.85,.85) for _ in range(3)],
            true_bias=[Interval(-.45,.45) for _ in range(3)])
        self.assertEqual(d['literal_event_kinds'],('prediction','S_zero'))
        self.assertTrue(d['all_literal_cells_consumed'])
        self.assertEqual(d['q15_prediction_blocks'],1)
        self.assertEqual(d['physical_S_blocks'],1)
        self.assertEqual(d['bias_driver_blocks'],1)
        self.assertTrue(d['same_history_source_ancestry_retained'])
        q=next(b for b in d['source_blocks'] if b.kind=='BRMM_q15_prediction')
        self.assertEqual(len(q.sectors),3)
        self.assertTrue(q.homogeneous_scale_fixed_one)
        self.assertIn('generator-1',q.token)
        self.assertIn('live-origin-3',q.token)
        # S=0 follows prediction, so its Joseph map must suffix-propagate the
        # earlier q15 forcing; the physical S response is appended afterwards.
        raw=L._q15_injection('H',I(1.7),I(.005))
        self.assertNotEqual(q.response,tuple(tuple(x for x in r) for r in raw))

    def test_A21_prediction_uses_one_bias_driver_for_error_and_truth(self):
        p=cell('A','prediction',0)
        d=L.compose_attached_sample(lineage('A',[p]),phi_true=I(.9999),tau_ba=I(1800))
        b=next(x for x in d['source_blocks'] if x.kind=='shared_bias_driver')
        B=b.response
        for i in range(3):
            self.assertTrue(B[18+i][i].contains(1.0))
            self.assertTrue(B[21+i][i].contains(1.0))
        self.assertEqual(d['q15_prediction_blocks'],1)

    def test_accelerometer_and_magnetometer_suffix_propagate_existing_sources(self):
        cells=[cell('H','prediction',0),cell('H','accelerometer',1),cell('H','magnetometer',2)]
        d=L.compose_attached_sample(
            lineage('H',cells),phi_true=I(.9999),
            held_bias_error=[Interval(-.85,.85) for _ in range(3)],
            true_bias=[Interval(-.45,.45) for _ in range(3)])
        self.assertEqual(d['literal_event_kinds'],('prediction','accelerometer','magnetometer'))
        self.assertEqual(d['q15_prediction_blocks'],1)
        self.assertEqual(len(d['J_joint24']),24)
        self.assertEqual(len(d['J_joint24'][0]),24)

    def test_contract_is_real_progress_but_still_blocks_storage(self):
        d=L.build()
        self.assertEqual(L.validate(d),[])
        self.assertTrue(d['source_uniform_attached_sample_cocycle_available'])
        self.assertTrue(d['joint_q15_quadratic_sector_attached_to_same_block'])
        self.assertFalse(d['multi_sample_COMPLETE_BRMM_word_composed'])
        self.assertFalse(d['storage_search_allowed'])
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':
    unittest.main()
