from dataclasses import replace
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_frontend_state_step as mod  # noqa: E402


class Sea3FrontEndStateStepTest(unittest.TestCase):
    def _sample(self):
        return mod.Sample(
            mod.MAHONY.Vec3(
                mod.MAHONY.I(0.01), mod.MAHONY.I(-0.02), mod.MAHONY.I(0.005)
            ),
            mod.MAHONY.Vec3(
                mod.MAHONY.I(0.2), mod.MAHONY.I(-0.1), mod.MAHONY.I(-9.75)
            ),
        )

    def _advance(self, state):
        return mod.advance(
            state,
            self._sample(),
            gravity_ms2=mod.MAHONY.I(9.80665),
            two_kp=mod.MAHONY.I(0.2),
            two_ki=mod.MAHONY.I(0.02),
        )

    def test_contract_is_same_source_and_fail_closed(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["shipping_source_parity_pass"])
        self.assertTrue(d["same_SEA3_sample_drives_Mahony_tuner_WPE"])
        self.assertTrue(d["private_Mahony_live_entry_invariant_consumed"])
        self.assertTrue(d["current_Riccati_schedule_exported_before_current_measurement"])
        self.assertTrue(d["actual_applied_per_axis_RS_exported_from_same_active_schedule"])
        self.assertTrue(d["tuner_consumes_previous_WPE_state"])
        self.assertTrue(d["same_current_vertical_acceleration_consumed_by_tuner_and_WPE"])
        self.assertFalse(d["source_generator"])
        self.assertFalse(d["independent_vertical_acceleration_input_allowed"])
        self.assertFalse(d["independent_wave_frequency_input_allowed"])
        self.assertFalse(d["independent_tuner_schedule_input_allowed"])
        self.assertFalse(d["target_binary32_WPE_libm_roundoff_closed"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])

    def test_pending_candidate_becomes_current_schedule_before_measurement(self):
        st = mod._point_state()
        tc = mod.TUNER.constants()
        candidate = mod.TUNER.CandidateState(
            mod.TUNER.I(1.3), mod.TUNER.I(0.6), mod.TUNER.I(3.0)
        )
        tuner = replace(
            st.tuner,
            candidate=candidate,
            scheduler=mod.TUNER.SchedulerState(mod.TUNER.I(0.05), True),
        )
        st = replace(st, tuner=tuner)
        succ = self._advance(st)
        self.assertTrue(succ)

        expected = mod.TUNER.ActiveSchedule(
            candidate.tau,
            candidate.sigma,
            candidate.rs,
            mod.TUNER.pseudo_period(candidate.tau, tc),
        )
        for x in succ:
            self.assertEqual(x.active_schedule_for_current_riccati_sample, expected)
            self.assertEqual(x.state.tuner.active, expected)
            self.assertEqual(
                x.actual_rs_std_xyz_for_current_riccati_sample,
                tuple(mod.TUNER.active_rs_std_xyz(expected, tc)),
            )

    def test_tuner_reads_previous_wpe_frequency_then_wpe_advances(self):
        st = mod._point_state()
        previous = mod.WPE.frequency_hz(st.wpe)
        succ = self._advance(st)
        self.assertTrue(succ)
        for x in succ:
            self.assertEqual(x.tuner_frequency_previous_wpe, previous)
            self.assertEqual(
                x.state.wpe.accel_prev,
                x.vertical_acceleration_current_sample,
            )

    def test_composite_preserves_component_validity_branches(self):
        d = mod.build()
        self.assertTrue(d["WPE_validity_branches_retained"])
        self.assertTrue(d["timer_boundary_branches_retained"])
        self.assertTrue(d["successor_branch_product_is_not_a_source_generator"])
        self.assertTrue(d["future_cell_split_required_if_branch_history_correlation_matters"])
        self.assertTrue(d["component_validation_pass"])
        self.assertEqual(
            d["shipping_order"],
            "commit pending -> private Mahony -> current Riccati uses committed schedule -> tuner(old WPE frequency) -> WPE update",
        )


if __name__ == "__main__":
    unittest.main()
