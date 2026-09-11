import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_typed_sample_source_cell_binding as B

I = Interval.point


def zmat(r,c): return [[I(0.0) for _ in range(c)] for _ in range(r)]
def eye(n): return [[I(1.0 if i==j else 0.0) for j in range(n)] for i in range(n)]


class TypedSampleSourceCellBindingTests(unittest.TestCase):
    def _fixture(self):
        Rwb=eye(3)
        f=(I(0.2),I(-0.1),I(-9.75))
        sample=KERNEL.SampleCoordinates(
            gyro_measurement=MAHONY.Vec3(I(.01),I(-.02),I(.005)),
            omega_body_corrected=(I(.009),I(-.019),I(.004)),
            specific_force=MAHONY.Vec3(*f),
            f_cog_body=f,
            R_wb=Rwb,
            due_S=False,
            aw_floor_requested=False,
            magnetometer_events_after_imu=(),
        )
        P=eye(18); R=eye(3); H=zmat(3,18)
        ev=KERNEL.RiccatiEventCell(
            mode='H',kind='accelerometer',event_index_in_sample=0,
            P_before=P,P_after=P,H=H,R=R,
        )
        sel=SELECTORS.PrefixSelector(
            sample_index=7,prefix_length=8,parent_branch_ordinal=0,successor_ordinal=0,
            parent_source_cell_id='parent',source_cell_id='child',sample_coordinates=sample,
            active_schedule=None,actual_rs_std_xyz=(I(1),I(1),I(1)),
            H_before=None,H_after=None,A_before=None,A_after=None,
            H_events_this_sample=('accelerometer',),A_events_this_sample=(),
            H_event_cells=(ev,),A_event_cells=(),H_floor_case=None,A_floor_case=None,
        )
        c=COVER.SourceCoverCell(
            source_token='child:e0',predecessor_token='parent',mode='H',sample_index=7,event_ordinal=0,
            kind='accelerometer',state=[I(0)]*18,P=P,dt_s=I(.005),tau_applied_s=I(1),sigma_aw_mps2=I(1),
            pseudo_elapsed_s=I(0),R=R,f_hat=f,R_hat=Rwb,radial_scale=I(.5),
            estimator_source_token='child',estimator_predecessor_token='parent',estimator_generated_coefficients=True,
        )
        return sel,[c]

    def test_status_remains_fail_closed(self):
        d=B.build()
        self.assertEqual(B.validate(d),[])
        self.assertTrue(d['trusted_prefix_selector_sample_coordinates_consumed'])
        self.assertFalse(d['production_complete_601_sample_binding_closed_here'])
        self.assertFalse(d['P4_PASS'])
        self.assertEqual(d['P3_delta'],1e-18)

    def test_exact_accelerometer_binding_accepts(self):
        s,c=self._fixture()
        self.assertEqual(B.validate_event_cells_against_selector(s,c,mode='H'),[])

    def test_detached_f_hat_rejected(self):
        s,c=self._fixture(); c[0]=copy.copy(c[0])
        c[0]=COVER.SourceCoverCell(**{**c[0].__dict__,'f_hat':(I(.3),I(-.1),I(-9.75))})
        self.assertTrue(any('f_hat detached' in x for x in B.validate_event_cells_against_selector(s,c,mode='H')))

    def test_detached_R_hat_rejected(self):
        s,c=self._fixture(); bad=eye(3);bad[0][0]=I(.9)
        c[0]=COVER.SourceCoverCell(**{**c[0].__dict__,'R_hat':bad})
        self.assertTrue(any('R_hat detached' in x for x in B.validate_event_cells_against_selector(s,c,mode='H')))

    def test_detached_covariance_rejected(self):
        s,c=self._fixture(); bad=eye(18);bad[0][0]=I(2)
        c[0]=COVER.SourceCoverCell(**{**c[0].__dict__,'P':bad})
        self.assertTrue(any('P is not trusted event P_before' in x for x in B.validate_event_cells_against_selector(s,c,mode='H')))

    def test_wrong_sample_index_rejected(self):
        s,c=self._fixture()
        c[0]=COVER.SourceCoverCell(**{**c[0].__dict__,'sample_index':8})
        self.assertTrue(any('sample index detached' in x for x in B.validate_event_cells_against_selector(s,c,mode='H')))


if __name__=='__main__': unittest.main()
