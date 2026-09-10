#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools';sys.path.insert(0,str(TOOLS));sys.path.insert(0,str(TOOLS/'stability'))
import ou3_brmm_hard_shaping_state as SHAPING
class HardShapingStateContractTest(unittest.TestCase):
 def test_contract_is_valid_but_provider_not_promoted(self):
  d=SHAPING.build();self.assertEqual(SHAPING.validate(d),[])
  self.assertTrue(d['reference_parameter_domain_compact']);self.assertTrue(d['compactness_is_not_an_open_obligation']);self.assertTrue(d['theorem_has_deterministic_shaping_contract']);self.assertTrue(d['theorem_has_explicit_hard_realization_set']);self.assertTrue(d['theorem_rejects_statistical_or_seeded_surrogates']);self.assertTrue(d['theorem_separates_probabilistic_random_sea_corollary']);self.assertTrue(d['complete_source_rejects_gaussian_word_generator']);self.assertTrue(d['correlated_outer_enclosure_route_used']);self.assertTrue(d['BRMM_to_correlated_outer_left_inclusion_closed']);self.assertFalse(d['global_physical_deployment_left_inclusion_closed_here']);self.assertFalse(d['exact_spectral_membership_oracle_required_for_P4'])
  b=d['sampled_behavior_target'];self.assertEqual(b['symbol'],'B^601_BRMM');self.assertEqual(b['correlated_outer_set_symbol'],'O^601_BRMM');self.assertTrue(b['compact']);self.assertTrue(b['validated_correlated_outer_enclosure_closed']);self.assertTrue(b['correlated_outer_left_inclusion_closed']);self.assertFalse(b['validated_membership_or_separation_oracle_closed'])
  o=d['joint_executor_coordinate_map'];self.assertTrue(o['closed']);self.assertTrue(o['raw_gyro_and_corrected_rate_distinct']);self.assertTrue(o['truth_attitude_and_nominal_R_hat_distinct']);self.assertTrue(o['truth_acceleration_and_nominal_a_w_hat_distinct']);self.assertFalse(o['sensor_forcing_hard_bound_closed_here']);self.assertFalse(o['BIAS0_assembled_sensor_qualification_closed_here'])
  self.assertTrue(d['hard_shaping_state_or_excitation_bound_closed']);self.assertFalse(d['complete_BRMM_family_materialized_here']);self.assertFalse(d['P3_promoted'])
 def test_forbidden_surrogates_cannot_close_xs(self):
  d=SHAPING.build()
  for k in ('power_spectrum_alone_is_hard_pathwise_bound','spectral_moments_alone_may_close_xs','gaussian_good_event_may_close_xs','replay_may_close_xs','seeded_128_frequency_generator_may_close_xs','finite_RAO_grid_may_close_xs','arbitrary_bounded_input_box_may_close_xs'):self.assertFalse(d[k],k)
 def test_outer_and_joint_output_routes_close_hard_shaping_without_deployment_admission(self):
  d=SHAPING.build();e=d['executable_ingredients'];self.assertTrue(e['correlated_outer_enclosure_closed']);self.assertTrue(e['BRMM_to_correlated_outer_left_inclusion_closed']);self.assertFalse(e['complete_BRMM_left_inclusion_closed']);self.assertTrue(e['joint_source_output_map_closed']);self.assertTrue(d['hard_shaping_state_or_excitation_bound_closed'])
if __name__=='__main__':unittest.main()
