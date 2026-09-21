import sys,unittest
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.signed_temporal import (
    atoms,
    balanced_gyro_weights,
    certificate,
    endpoint_annihilating_multiplier_constraints,
    endpoint_cancelled_source_bound,
    endpoint_jets_vanish,
    forced_adjoint_source_bound,
    gyro_zero_mean_companion,
    literal_signed_functional_bound,
    margin_to_Bstar_theorem,
    physical_bounds,
    projection_sector_gap,
    separated_reader_action_implication,
    source_margin_attempt,
    zero_terminal_homogeneous_adjoint_impossible,
)

class SignedTemporalTests(unittest.TestCase):
    def test_spline_atoms_and_endpoint_jets(self):
        k=(F(0),F(1),F(2),F(3))
        self.assertEqual(sum(atoms(k)),0)
        self.assertTrue(endpoint_jets_vanish(k))
    def test_homogeneous_adjoint_failure_is_exact(self):
        self.assertTrue(zero_terminal_homogeneous_adjoint_impossible(3))
    def test_physical_margin_positive(self):
        b=physical_bounds()
        self.assertGreater(b["physical_floor_margin"],0)
    def test_projection_gap_positive(self):
        self.assertGreater(projection_sector_gap(),0)
    def test_margin_to_Bstar_implication(self):
        t=margin_to_Bstar_theorem()
        self.assertTrue(t["implication_closed"])
        self.assertFalse(t["premise_margins_source_uniformly_certified"])
        b=separated_reader_action_implication(F(1,100),F(1,200),F(2),F(3),F(5),F(7),10)
        self.assertTrue(b["B_star_finite"])
        self.assertGreater(b["B_star_scalar_ceiling"],0)

    def test_source_margin_attempt_records_exact_gap(self):
        a=source_margin_attempt()
        self.assertGreater(a["physical_collinearity_reserve"],0)
        self.assertIsNone(a["forced_adjoint_signed_transfer_ceiling"])
        self.assertIsNone(a["Delta_col_lower"])
        self.assertFalse(a["new_physical_assumption_needed"])

    def test_literal_adjoint_eliminates_innovation_energy(self):
        b=forced_adjoint_source_bound()
        self.assertFalse(b["innovation_energy_needed"])
        self.assertTrue(b["BA_endpoint_bounded"])
        self.assertFalse(b["BG_endpoint_bounded"])
        self.assertFalse(b["AW_endpoint_bounded"])
        self.assertIsNone(b["finite_numeric_ceiling"])
        self.assertEqual(literal_signed_functional_bound(F(2),F(3),F(4),F(5)),F(32))

    def test_endpoint_annihilation(self):
        a=endpoint_annihilating_multiplier_constraints()
        self.assertTrue(a["four_S_spline_satisfies_AW_endpoint_conditions"])
        b=balanced_gyro_weights([1,2,4,8])
        self.assertEqual(sum(b),0)
        g=gyro_zero_mean_companion(b)
        self.assertTrue(g["endpoint_zero"])
        self.assertEqual(endpoint_cancelled_source_bound(2,3,4,F(1,100000),64),
                         F(6)+F(256,100000))

    def test_fail_closed(self):
        c=certificate()
        self.assertFalse(c["source_uniform_nominal_force_field_temporal_margin"])
        self.assertFalse(c["source_uniform_nominal_gyro_alias_temporal_margin"])
        self.assertFalse(c["uniform_historical_AG_readout_action"])
        self.assertFalse(c["full_21_covariance_upper"])
        self.assertFalse(c["rho0_certified"])

if __name__=="__main__": unittest.main()
