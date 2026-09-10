import dataclasses
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_typed_sample_nonlinear_chain as CHAIN

I = Interval.point


def zmat(r, c):
    return [[I(0.0) for _ in range(c)] for _ in range(r)]


def eye(n):
    return [[I(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]


class TypedSampleNonlinearChainTests(unittest.TestCase):
    def _fixture(self, mode="H"):
        n = 18 if mode == "H" else 21
        Rwb = eye(3)
        f = (I(0.2), I(-0.1), I(-9.75))
        omega = (I(0.009), I(-0.019), I(0.004))
        sample = KERNEL.SampleCoordinates(
            gyro_measurement=MAHONY.Vec3(I(.01), I(-.02), I(.005)),
            omega_body_corrected=omega,
            specific_force=MAHONY.Vec3(*f),
            f_cog_body=f,
            R_wb=Rwb,
            due_S=False,
            aw_floor_requested=False,
            magnetometer_events_after_imu=(),
        )
        P = eye(n)
        R = eye(3)
        H = zmat(3, n)
        pred_cell = KERNEL.RiccatiEventCell(
            mode=mode, kind="prediction", event_index_in_sample=0,
            P_before=P, P_after=P, F=eye(n), Q=zmat(n, n),
        )
        acc_cell = KERNEL.RiccatiEventCell(
            mode=mode, kind="accelerometer", event_index_in_sample=1,
            P_before=P, P_after=P, H=H, R=R,
        )
        events = ("prediction", "accelerometer")
        sel = SELECTORS.PrefixSelector(
            sample_index=7, prefix_length=8, parent_branch_ordinal=0, successor_ordinal=0,
            parent_source_cell_id="parent", source_cell_id="child", sample_coordinates=sample,
            active_schedule=None, actual_rs_std_xyz=(I(1), I(1), I(1)),
            H_before=None, H_after=None, A_before=None, A_after=None,
            H_events_this_sample=events if mode == "H" else (),
            A_events_this_sample=events if mode == "A" else (),
            H_event_cells=(pred_cell, acc_cell) if mode == "H" else (),
            A_event_cells=(pred_cell, acc_cell) if mode == "A" else (),
            H_floor_case=None, A_floor_case=None,
        )
        state0 = [I(0.0)] * n
        tau_ba = I(1000.0) if mode == "A" else None
        pred = PRED.prediction_event(mode, state0, omega, I(.005), I(1.0), tau_ba=tau_ba)
        common = dict(
            predecessor_token="parent", mode=mode, sample_index=7,
            P=P, dt_s=I(.005), tau_applied_s=I(1.0), sigma_aw_mps2=I(1.0),
            pseudo_elapsed_s=I(0.0), radial_scale=I(.5),
            estimator_source_token="child", estimator_predecessor_token="parent",
            estimator_generated_coefficients=True,
        )
        if mode == "A":
            common.update(true_bias=(I(0.0), I(0.0), I(0.0)), bias_projection_limit=.5)
        cells = [
            COVER.SourceCoverCell(
                source_token="child:e0", event_ordinal=0, kind="prediction",
                state=state0, **common,
            ),
            COVER.SourceCoverCell(
                source_token="child:e1", event_ordinal=1, kind="accelerometer",
                state=pred["state_out"], R=R, f_hat=f, R_hat=Rwb, **common,
            ),
        ]
        return sel, cells, tau_ba

    def test_status_is_nonpromoting(self):
        d = CHAIN.build()
        self.assertEqual(CHAIN.validate(d), [])
        self.assertTrue(d["exact_finite_prediction_state_and_Jacobian_composed"])
        self.assertTrue(d["same_cell_Joseph_state_and_Jacobian_composed"])
        self.assertFalse(d["production_complete_601_sample_chain_materialized_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_prediction_and_joseph_compose_from_same_typed_sample(self):
        selector, cells, _ = self._fixture("H")
        out = CHAIN.materialize_sample_chain(selector, cells, mode="H")
        self.assertEqual(out["event_kinds"], ("prediction", "accelerometer"))
        self.assertEqual(out["prediction_count"], 1)
        self.assertEqual(out["accelerometer_update_count"], 1)
        self.assertTrue(out["literal_state_succession_exact"])
        self.assertFalse(out["prediction_identity_placeholder_used"])
        self.assertFalse(out["independent_event_Jacobian_box_used"])
        self.assertEqual(len(out["J_word"]), 18)
        self.assertEqual(len(out["state_out"]), 18)

    def test_A21_prediction_and_projection_compose_from_same_state_chain(self):
        selector, cells, tau_ba = self._fixture("A")
        out = CHAIN.materialize_sample_chain(selector, cells, mode="A", tau_ba=tau_ba)
        self.assertEqual(out["event_kinds"], ("prediction", "accelerometer"))
        self.assertEqual(out["prediction_count"], 1)
        self.assertEqual(out["accelerometer_update_count"], 1)
        self.assertEqual(len(out["J_word"]), 21)
        self.assertEqual(len(out["state_out"]), 21)

    def test_A21_rejects_missing_configured_bias_time_constant(self):
        selector, cells, _ = self._fixture("A")
        with self.assertRaisesRegex(ValueError, "configured tau_ba"):
            CHAIN.materialize_sample_chain(selector, cells, mode="A")

    def test_detached_next_event_state_is_rejected(self):
        selector, cells, _ = self._fixture("H")
        bad = list(cells[1].state)
        bad[6] = I(0.25)
        cells[1] = COVER.SourceCoverCell(**{**cells[1].__dict__, "state": bad})
        with self.assertRaisesRegex(ValueError, "detached from previous exact/outward state_out"):
            CHAIN.materialize_sample_chain(selector, cells, mode="H")

    def test_prediction_uses_typed_corrected_rate_not_raw_measurement(self):
        selector, cells, _ = self._fixture("H")
        out = CHAIN.materialize_sample_chain(selector, cells, mode="H")
        altered_sample = KERNEL.SampleCoordinates(
            gyro_measurement=selector.sample_coordinates.gyro_measurement,
            omega_body_corrected=(I(.109), I(-.019), I(.004)),
            specific_force=selector.sample_coordinates.specific_force,
            f_cog_body=selector.sample_coordinates.f_cog_body,
            R_wb=selector.sample_coordinates.R_wb,
            due_S=False,
            aw_floor_requested=False,
            magnetometer_events_after_imu=(),
        )
        altered = dataclasses.replace(selector, sample_coordinates=altered_sample)
        with self.assertRaisesRegex(ValueError, "detached from previous exact/outward state_out"):
            CHAIN.materialize_sample_chain(altered, cells, mode="H")
        self.assertEqual(out["prediction_count"], 1)


if __name__ == "__main__":
    unittest.main()
