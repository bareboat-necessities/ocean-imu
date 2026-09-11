#!/usr/bin/env python3
from __future__ import annotations
import copy,sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
from ou3_interval import Interval
import ou3_p4_h18_binary32_prefix_response as FP
import ou3_p4_source_uniform_estimator_event_attachment as A
import ou3_brmm_complete_window_execution_kernel as K
import ou3_brmm_full_normal_live_word as W
I=FP.I
def ident(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def attached():
    js=A.JOINT._smoke_state();s=K.SampleCoordinates(gyro_measurement=K.MAHONY.Vec3(I(.01),I(-.02),I(.005)),omega_body_corrected=(I(.01),I(-.02),I(.005)),specific_force=K.MAHONY.Vec3(I(.2),I(-.1),I(-9.75)),f_cog_body=(I(0),I(0),I(-9.80665)),R_wb=ident(3),due_S=True,aw_floor_requested=True,magnetometer_events_after_imu=(K.MagneticEvent((I(20),I(0),I(40))),))
    b=K.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=W.initialize_word('H',ident(18)),A=W.initialize_word('A',ident(21)),source_cell_id='root')
    rows=A.synchronize_sample(branch=b,joint_state=js,sample=s,state_in_H=[I(0) for _ in range(18)],state_in_A=[I(0) for _ in range(21)],radial_scale=Interval(0,1),true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4,tau_ba=I(1800),sample_index=0,next_cell_prefix='fp-h18')
    if not rows:raise RuntimeError('attachment empty')
    return rows[0][0]
class H18Binary32PrefixResponseTest(unittest.TestCase):
    def test_contract_conditional_not_deployment(self):
        d=FP.build();self.assertEqual(FP.validate(d),[]);self.assertTrue(d['H18_every_prefix_binary32_response_materializable']);self.assertFalse(d['packet_count_times_worst_roundoff_used']);self.assertFalse(d['deployment_toolchain_qualification_required_for_mathematical_P4']);self.assertFalse(d['P4_PASS'])
    def test_one_event_appends_one_18D_block(self):
        h=attached();rows=FP.materialize([h]);self.assertEqual(len(rows),len(h.cells))
        for i,r in enumerate(rows,1):
            self.assertEqual(FP.shape(r.response),(18,18*i));self.assertEqual(r.normalized_input_dimension,18*i)
    def test_first_block_is_local_bound_times_identity(self):
        h=attached();r=FP.materialize([h])[0];d=r.local_bound
        for i in range(18):
            self.assertLessEqual(r.response[i][i].lo,d);self.assertGreaterEqual(r.response[i][i].hi,d)
if __name__=='__main__':unittest.main()
