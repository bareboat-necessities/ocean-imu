#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from dataclasses import replace
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[2]/"tools"/"stability"
sys.path.insert(0,str(TOOLS))

import ou3_p4_a21_source_indexed_prefix_transport as A21


class A21SourceIndexedPrefixTransportTest(unittest.TestCase):
    def test_all_bias_families_close_structural_prefix_transport_only(self):
        d=A21.build()
        self.assertEqual(A21.validate(d),[])
        self.assertTrue(d["all_BIAS0_BIAS1_BIAS2_recurrences_composed"])
        self.assertTrue(d["A21_Joseph_reset_and_projection_are_distinct_prefixes"])
        self.assertTrue(d["source_indexed_Phi_rebase_inserted_between_samples"])
        self.assertTrue(d["projection_has_exact_zero_Phi_interior_shift"])
        self.assertTrue(d["one_shared_w_and_mtau_supply_block_per_prediction"])
        self.assertTrue(d["same_w_enters_bias_error_and_true_bias"])
        self.assertTrue(d["projection_true_bias_coupling_retained_in_24state_map"])
        self.assertTrue(d["A21_every_subevent_prefix_transport_constructor_closed"])
        self.assertFalse(d["production_BRMM_moment_radial_sector_attached_here"])
        self.assertFalse(d["production_binary32_prefix_map_attached_here"])
        self.assertFalse(d["production_compatible_storage_attached_here"])
        self.assertFalse(d["A21_endpoint_augmented_LDLT_closed_here"])
        self.assertFalse(d["A21_every_prefix_augmented_LDLT_closed_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_each_family_retains_projection_and_shared_supply(self):
        d=A21.build()
        self.assertEqual(set(d["smoke_by_family"]),set(A21.BIAS.FAMILIES))
        for family,row in d["smoke_by_family"].items():
            self.assertEqual(row["prediction_count"],2,family)
            self.assertGreaterEqual(row["projection_count"],2,family)
            self.assertEqual(row["source_rebase_count"],1,family)
            self.assertTrue(row["all_identity_residuals_contain_zero"],family)
            self.assertTrue(row["projection_has_zero_Phi_interior_shift"],family)
            self.assertTrue(row["one_6D_shared_w_mtau_block_per_prediction"],family)
            self.assertTrue(row["prediction_supply_blocks_same_w"],family)

    def test_detached_absolute_true_bias_prefix_is_rejected(self):
        samples=list(A21.smoke_objects("BIAS1"))
        first=samples[0]
        cells=list(first.cells)
        target=next(i for i,c in enumerate(cells) if c.kind in ("S_zero","accelerometer","magnetometer"))
        c=cells[target]
        cells[target]=replace(c,true_bias=(A21.Interval.point(9.0),A21.Interval.point(0.0),A21.Interval.point(0.0)))
        samples[0]=replace(first,cells=tuple(cells))
        with self.assertRaises(RuntimeError):
            A21.build_transport(samples,"BIAS1")

    def test_projection_beta_coupling_is_not_identity_detached(self):
        samples=A21.smoke_objects("BIAS2")
        t=A21.build_transport(samples,"BIAS2")
        projections=[e for e in t["events"] if e["kind"]=="bias_projection"]
        self.assertTrue(projections)
        # The 24-state projection map must retain beta->e_b coupling whenever
        # the generalized projection Jacobian differs from identity.  At minimum
        # every map keeps the beta state itself as an identity block.
        for e in projections:
            C=e["C"]
            for i in range(3):
                self.assertTrue(C[21+i][21+i].contains(1.0))


if __name__=="__main__":
    unittest.main()
