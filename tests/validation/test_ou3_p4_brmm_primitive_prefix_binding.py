import dataclasses
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_brmm_primitive_prefix_binding as B


class PrimitivePrefixBindingTests(unittest.TestCase):
    def test_status_is_fail_closed(self):
        d = B.build()
        self.assertEqual(B.validate(d), [])
        self.assertTrue(d["same_provider_transition_owns_primitives_and_acceleration_moments"])
        self.assertTrue(d["cross_sample_primitive_out_equals_next_primitive_in_required"])
        self.assertFalse(d["production_complete_provider_lineages_bound_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_reference_metadata_closes(self):
        b = B._smoke_bindings()
        self.assertEqual(B.validate_binding_metadata(b), [])

    def test_detached_moment_witness_rejected(self):
        b = B._smoke_bindings()
        b[1] = dataclasses.replace(b[1], moment_witness_id="other")
        self.assertTrue(any("moments detached" in x for x in B.validate_binding_metadata(b)))

    def test_same_sample_transition_change_rejected(self):
        b = B._smoke_bindings()
        b[1] = dataclasses.replace(b[1], primitive_out_id="wrong")
        self.assertTrue(any("same-sample" in x for x in B.validate_binding_metadata(b)))

    def test_cross_sample_primitive_break_rejected(self):
        b = B._smoke_bindings()
        b[2] = dataclasses.replace(b[2], primitive_in_id="not-q1")
        self.assertTrue(any("primitive_out is not next primitive_in" in x for x in B.validate_binding_metadata(b)))

    def test_centered_origin_change_rejected(self):
        b = B._smoke_bindings()
        b[2] = dataclasses.replace(b[2], centered_S_origin_witness_id="new-origin")
        self.assertTrue(any("centered-S origin changed" in x for x in B.validate_binding_metadata(b)))

    def test_sample_skip_rejected(self):
        b = B._smoke_bindings()
        b[2] = dataclasses.replace(b[2], sample_index=2)
        self.assertTrue(any("skipped a sample" in x for x in B.validate_binding_metadata(b)))

    def test_duplicate_event_rejected(self):
        b = B._smoke_bindings()
        b[1] = dataclasses.replace(b[1], event_token=b[0].event_token)
        self.assertTrue(any("duplicate literal event token" in x for x in B.validate_binding_metadata(b)))

    def test_false_promotion_rejected(self):
        d = B.build()
        d["P4_PASS"] = True
        self.assertIn("P4_PASS not false", B.validate(d))


if __name__ == "__main__":
    unittest.main()
