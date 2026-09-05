#!/usr/bin/env python3
import copy
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_frontend_state_step as FRONTEND
import ou3_sea3_window_artifact_codec as CODEC
from ou3_interval import Interval


class Sea3WindowArtifactCodecTest(unittest.TestCase):
    def test_status_roundtrips_without_promoting(self):
        d = CODEC.build()
        self.assertEqual(CODEC.validate(d), [])
        self.assertTrue(d["strict_codec_ready"])
        self.assertTrue(d["requires_provider_acceptance_before_canonical_use"])
        self.assertFalse(d["establishes_source_reachability"])
        self.assertFalse(d["precomputed_aw_floor_increment_accepted"])
        self.assertFalse(d["P3_promoted"])

    def test_frontend_roundtrip_is_exact(self):
        state = FRONTEND._point_state()
        payload = CODEC.frontend_to_json(state, "front")
        self.assertEqual(CODEC.frontend_from_json(payload, "front"), state)
        with self.assertRaises(ValueError):
            CODEC.frontend_from_json(payload, "other")

    def test_interval_transport_does_not_rewiden_exact_zero(self):
        z = CODEC.interval_from_json([0.0, 0.0], "z")
        self.assertEqual(z, Interval.point(0.0))
        self.assertEqual(CODEC.interval_to_json(z), [0.0, 0.0])

    def test_interval_rejects_nonfinite_or_reversed_endpoints(self):
        for bad in ([1.0, -1.0], [0.0, math.inf], [math.nan, 1.0], [True, 1.0]):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    CODEC.interval_from_json(bad)

    def test_live_seed_requires_witness_symmetry_and_spd(self):
        PH = CODEC._diag(18, 2.0)
        PA = CODEC._diag(21, 2.0)
        payload = CODEC.live_seed_to_json(PH, PA, "seed")
        h, a = CODEC.live_seed_from_json(payload, "seed")
        self.assertEqual(h, PH)
        self.assertEqual(a, PA)

        bad = copy.deepcopy(payload)
        bad["P0_H_interval"][0][0] = [-1.0, -1.0]
        with self.assertRaises(ValueError):
            CODEC.live_seed_from_json(bad, "seed")

        bad = copy.deepcopy(payload)
        bad["P0_H_interval"][0][1] = [0.1, 0.1]
        with self.assertRaises(ValueError):
            CODEC.live_seed_from_json(bad, "seed")

    def test_sample_rejects_precomputed_floor_increment(self):
        d = CODEC.build()
        self.assertTrue(d["strict_codec_ready"])
        # Recreate the codec's fixture through its public serializers.
        z = Interval.point(0.0)
        import ou3_sea3_complete_window_execution_kernel as KERNEL
        import ou3_sea3_private_mahony_state_step as MAHONY
        sample = KERNEL.SampleCoordinates(
            gyro_measurement=MAHONY.Vec3(z, z, z),
            omega_body_corrected=(z, z, z),
            specific_force=MAHONY.Vec3(z, z, Interval.point(-9.8)),
            f_cog_body=(z, z, Interval.point(-9.8)),
            R_wb=[[Interval.point(1.0 if i == j else 0.0) for j in range(3)] for i in range(3)],
            due_S=False,
            aw_floor_requested=True,
        )
        physical, events = CODEC.sample_to_payload(
            sample,
            transition_witness_id="tr",
            joint_response_witness_id="resp",
        )
        events["aw_covariance_floor_increment"] = [[[0.0, 0.0]]]
        with self.assertRaises(ValueError):
            CODEC.sample_from_transition({
                "source_transition_witness_id": "tr",
                "joint_response_witness_id": "resp",
                "joint_physical_output": physical,
                "source_events": events,
            })


if __name__ == "__main__":
    unittest.main()
