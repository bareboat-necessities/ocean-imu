#!/usr/bin/env python3
from __future__ import annotations
import copy,sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
from ou3_interval import Interval
import ou3_p4_brmm_physical_source_prefix_response as R
import ou3_p4_source_uniform_estimator_event_attachment as A
import ou3_brmm_complete_window_execution_kernel as K
import ou3_brmm_full_normal_live_word as W

I=R.I
def ident(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def sample():return K.SampleCoordinates(gyro_measurement=K.MAHONY.Vec3(I(.01),I(-.02),I(.005)),omega_body_corrected=(I(.01),I(-.02),I(.005)),specific_force=K.MAHONY.Vec3(I(.2),I(-.1),I(-9.75)),f_cog_body=(I(0),I(0),I(-9.80665)),R_wb=ident(3),due_S=True,aw_floor_requested=True,magnetometer_events_after_imu=(K.MagneticEvent((I(20),I(0),I(40))),))
def pair():
    js=A.JOINT._smoke_state();b=K.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=W.initialize_word('H',ident(18)),A=W.initialize_word('A',ident(21)),source_cell_id='root')
    rows=A.synchronize_sample(branch=b,joint_state=js,sample=sample(),state_in_H=[I(0) for _ in range(18)],state_in_A=[I(0) for _ in range(21)],radial_scale=Interval(0,1),true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4,tau_ba=I(1800),sample_index=0,next_cell_prefix='phys')
    if not rows:raise RuntimeError('attachment smoke empty')
    return rows[0]

class PhysicalSourcePrefixResponseTest(unittest.TestCase):
    def test_contract(self):
        d=R.build();self.assertEqual(R.validate(d),[]);self.assertTrue(d['H18_physical_source_prefix_response_materializable']);self.assertTrue(d['A21_all_BIAS0_BIAS1_BIAS2_physical_source_prefix_response_materializable']);self.assertFalse(d['hard_entry_radial_used_as_physical_source_amplitude']);self.assertFalse(d['P4_PASS'])
    def test_H18_one_prediction_appends_one_15D_witness(self):
        h,_=pair();rows=R.materialize([h],mode='H');self.assertTrue(rows);self.assertTrue(all(len(x.prediction_witness_ids)==1 for x in rows));self.assertTrue(all(R.shape(x.response)==(18,15) for x in rows))
    def test_A21_all_bias_families_keep_true_bias_rows_unforced_at_prediction(self):
        _,a=pair()
        for fam in ('BIAS0','BIAS1','BIAS2'):
            rows=R.materialize([a],mode='A',family=fam);self.assertTrue(rows);first=rows[0];self.assertEqual(R.shape(first.response),(24,15));self.assertTrue(all(v.lo==0 and v.hi==0 for row in first.response[21:24] for v in row))
    def test_H_mode_rejects_bias_family(self):
        h,_=pair()
        with self.assertRaises(ValueError):R.materialize([h],mode='H',family='BIAS0')

if __name__=='__main__':unittest.main()
