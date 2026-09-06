import importlib.util
import math
from pathlib import Path
import struct
import sys
import unittest
ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
STABILITY = TOOLS / "stability"
sys.path.insert(0, str(TOOLS)); sys.path.insert(0, str(STABILITY))
spec = importlib.util.spec_from_file_location("ou3_source_domain_contract", STABILITY / "ou3_source_domain_contract.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
import ou3_sea3_complete_source as sea3_complete
import ou3_sea3_p1_compatibility as sea3_p1
import ou3_sea3_physical_admissibility as sea3_phys
import ou3_sea3_wave_period_spectral_identity as sea3_period_identity
def f32(value): return struct.unpack(">f", struct.pack(">f", float(value)))[0]
class SourceDomainContractTests(unittest.TestCase):
 def test_contract_uses_shipping_clamps_and_keeps_theorem_unpromoted(self):
  d=mod.build(mod.DEFAULT_HEADER); self.assertEqual(d["schema"],3); self.assertTrue(d["source_generated_not_trajectory_fit"]); self.assertTrue(d["source_complete_parameter_domain"]); self.assertFalse(d["validated_arithmetic"]); self.assertFalse(d["outward_rounded"]); self.assertEqual(d["implementation_scalar_semantics"]["type"],"IEEE754_BINARY32"); self.assertEqual(d["continuous_parameters"]["tau_aw_s"],[f32(.02),f32(12.)]); self.assertEqual(d["continuous_parameters"]["sigma_aw_mps2"],[f32(.05),f32(4.)]); self.assertEqual(set(d["discrete_source_branches"]["mode"]),{"H","A"})
 def test_constexpr_arithmetic_rounds_as_binary32_after_each_operation(self):
  text='''constexpr float A = 0.1f;\nconstexpr float B = A + A;\nconstexpr float C = B + A;\nconstexpr float DT = 1.0f / 200.0f;'''; a=f32(.1); b=f32(a+a); c=f32(b+a); dt=f32(f32(1.)/f32(200.)); self.assertEqual(mod.parse_const(text,"A"),a); self.assertEqual(mod.parse_const(text,"B"),b); self.assertEqual(mod.parse_const(text,"C"),c); self.assertEqual(mod.parse_const(text,"DT"),dt); self.assertNotEqual(dt,.005)
 def test_validated_parameter_box_outwardly_contains_every_source_endpoint(self):
  d=mod.build(mod.DEFAULT_HEADER); box=d["validated_parameter_box"]; self.assertTrue(box["validated_arithmetic"]); self.assertTrue(box["outward_rounded"]); self.assertEqual(box["theorem_promotion"],"NOT_ESTABLISHED"); self.assertFalse(box["continuous_word_enclosed"]); self.assertFalse(box["nonlinear_word_enclosed"])
  for name,b in d["continuous_parameters"].items(): lo,hi=box["continuous_parameters"][name]; self.assertEqual(lo,math.nextafter(b[0],-math.inf)); self.assertEqual(hi,math.nextafter(b[1],math.inf))
  for name,v in d["timing_constants_s"].items(): lo,hi=box["timing_constants_s"][name]; self.assertEqual(lo,math.nextafter(v,-math.inf)); self.assertEqual(hi,math.nextafter(v,math.inf))
 def test_configured_runtime_sampling_assumption_is_explicit_and_source_bound(self):
  d=mod.build(mod.DEFAULT_HEADER); r=d["configured_runtime_assumption"]; expected=f32(f32(1.)/f32(200.)); self.assertEqual(r["qualification"],"CONFIGURED_VALIDATION_RUNTIME_ASSUMPTION"); self.assertEqual(r["sample_period_contract"],"FIXED_SOURCE_NOMINAL"); self.assertEqual(r["imu_dt_s"],expected); self.assertFalse(r["api_enforces_this_bound"]); self.assertEqual(d["validated_parameter_box"]["configured_runtime"],r)
 def test_contract_names_every_hybrid_transition_required_for_deployment(self):
  d=mod.build(mod.DEFAULT_HEADER); self.assertEqual(set(d["hybrid_obligations"]),{"startup_handoff","held_to_active","magnetic_lock","magnetic_regauge_refinement","tilt_reset","tilt_relock","cooldown_reentry","periodic_aw_covariance_sync"}); self.assertEqual(d["periodic_aw_covariance_sync_proof"]["required_mode"],"PSD_NONEXPANSIVE")
 def test_physical_height_period_coupling_remains_fail_closed(self):
  d=sea3_phys.build(); self.assertEqual(sea3_phys.validate(d),[]); self.assertTrue(d["three_partition_contract"]["independent_H_r_and_T_p_rectangular_extrema_forbidden"]); self.assertEqual(d["repository_total_Hs_upper_m"],8.5); self.assertFalse(d["left_language_inclusion_closed"])
 def test_cartesian_sea_x_rao_domain_is_rejected_before_p1(self):
  d=sea3_p1.build(); self.assertEqual(sea3_p1.validate(d),[]); self.assertTrue(d["cartesian_product_refuted_by_analytical_witness"]); self.assertTrue(d["coupled_SEA3_domain_required"]); self.assertFalse(d["independent_cartesian_sea_x_RAO_domain_is_P1_sound"]); self.assertFalse(d["finite_window_realization_certificate_closed"]); self.assertFalse(d["L_actual_sea_subset_Lhat_SEA3_closed"])
 def test_wave_period_leak_subtraction_identity_stays_source_bound(self):
  d=sea3_period_identity.build(); self.assertEqual(sea3_period_identity.validate(d),[]); self.assertEqual(set(d["source_parity"]),set(sea3_period_identity.SOURCE_PARITY_KEYS)); self.assertTrue(all(d["source_parity"].values()))
 def test_complete_sea3_conditional_source_does_not_claim_physical_left_inclusion(self):
  d=sea3_complete.build(); self.assertEqual(sea3_complete.validate(d),[]); self.assertEqual(d["canonical_P3_source"],"COMPLETE_SEA3_NORMAL_LIVE_WORD"); self.assertTrue(d["P3_source_contract_ready"]); self.assertFalse(d["P3_source_family_materialized"]); self.assertFalse(d["global_physical_deployment_left_inclusion_closed_here"])
if __name__ == "__main__": unittest.main()
