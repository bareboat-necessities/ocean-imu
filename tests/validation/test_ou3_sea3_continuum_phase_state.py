#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TOOLS / "stability"))

import ou3_sea3_continuum_phase_state as PHASE


class ContinuumPhaseStateTest(unittest.TestCase):
    def test_continuum_phase_contract_closes_without_grid(self) -> None:
        d = PHASE.build()
        self.assertEqual(PHASE.validate(d), [])
        self.assertTrue(d["continuum_index_set_retained"])
        self.assertTrue(d["continuum_phase_coordinate_set_closed"])
        self.assertTrue(d["phase_continuous_propagation_closed"])
        self.assertFalse(d["finite_frequency_grid_used"])
        self.assertFalse(d["finite_direction_grid_used"])
        self.assertFalse(d["seeded_phase_realization_used"])
        self.assertFalse(d["phase_reset_on_lambda_transition_allowed"])
        self.assertFalse(d["hard_spectral_driver_set_closed"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])

    def test_rotation_preserves_circle(self) -> None:
        for omega in (0.0, 0.1, 1.0, 2.0 * math.pi * 6.0):
            for phase in (-2.7, -0.4, 0.0, 1.9):
                q0 = (math.cos(phase), math.sin(phase))
                q1 = PHASE.rotation_step(q0[0], q0[1], omega, 0.005)
                self.assertAlmostEqual(PHASE.norm2(q0), PHASE.norm2(q1), places=14)

    def test_invalid_frequency_or_step_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            PHASE.rotation_step(1.0, 0.0, -1.0, 0.005)
        with self.assertRaises(ValueError):
            PHASE.rotation_step(1.0, 0.0, 1.0, -0.005)


if __name__ == "__main__":
    unittest.main()
