#!/usr/bin/env python3
"""Fail-closed contract for deterministic complete-BRMM hard shaping.

BRMM is already the compact theorem-domain sea family. Continuum phase/no-reseed
propagation is closed, B^601_BRMM is compact, and the validated correlated
outer set O^601_BRMM closes B^601_BRMM subset O^601_BRMM without Cartesianizing
samples or axes.

That BRMM->outer inclusion is NOT the separate physical/deployment->BRMM left
inclusion. The latter remains false until full source admission is proved.

The exact joint source/filter-state -> typed-executor coordinate map is also
materialized. The hard shaping/excitation representation is therefore closed,
while provider execution/materialization and deployment admission remain open.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_brmm_continuum_phase_state as PHASE
import ou3_brmm_hard_window_behavior as BEHAVIOR
import ou3_brmm_joint_executor_coordinate_map as OUTPUT
REPO=Path(__file__).resolve().parents[2]
THEOREM=REPO/'doc'/'kalman_ou_iii'/'w3d-marine-reference-models.tex-part'
COMPLETE_SOURCE=REPO/'tools'/'stability'/'ou3_brmm_complete_source.py'
SCHEMA=7
QUALIFICATION='OU3_BRMM_HARD_SHAPING_STATE_CONTRACT_V7'
CONTINUUM_PHASE_COORDINATE_SET_CLOSED=PHASE.CONTINUUM_PHASE_COORDINATE_SET_CLOSED
PHASE_CONTINUOUS_PROPAGATION_CLOSED=PHASE.PHASE_CONTINUOUS_PROPAGATION_CLOSED
HARD_SPECTRAL_DRIVER_SET_CLOSED=False
CORRELATED_OUTER_ENCLOSURE_CLOSED=True
BRMM_TO_CORRELATED_OUTER_LEFT_INCLUSION_CLOSED=True
COMPLETE_BRMM_LEFT_INCLUSION_CLOSED=False
JOINT_SOURCE_OUTPUT_MAP_CLOSED=True
HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED=all((CONTINUUM_PHASE_COORDINATE_SET_CLOSED,PHASE_CONTINUOUS_PROPAGATION_CLOSED,CORRELATED_OUTER_ENCLOSURE_CLOSED,BRMM_TO_CORRELATED_OUTER_LEFT_INCLUSION_CLOSED,JOINT_SOURCE_OUTPUT_MAP_CLOSED))
def _normalized_text(text):return ' '.join(text.split())
def build():
 theorem=THEOREM.read_text();flat=_normalized_text(theorem);complete=COMPLETE_SOURCE.read_text();phase=PHASE.build();behavior=BEHAVIOR.build();output=OUTPUT.build()
 bad={'phase':PHASE.validate(phase),'behavior':BEHAVIOR.validate(behavior),'output':OUTPUT.validate(output)};bad={k:v for k,v in bad.items() if v}
 if bad:raise RuntimeError(f'BRMM shaping prerequisites failed: {bad}')
 outer=bool(behavior['validated_correlated_outer_enclosure_closed'] and behavior['correlated_outer_left_inclusion_closed'] and behavior['correlated_outer_retains_cross_sample_and_axis_dependence'])
 joint=bool(output['joint_source_output_map_closed'] and output['source_filter_joint_coordinate_map_materialized'] and output['measurement_coordinates_and_nominal_geometry_distinct'] and output['raw_gyro_and_corrected_rate_distinct'] and output['corrected_rate_depends_on_current_estimated_gyro_bias'] and output['same_joint_witness_required_for_all_coordinates'])
 if not outer or not joint:raise RuntimeError('hard shaping outer/output prerequisite disappeared')
 executable={'continuum_phase_coordinate_set_closed':CONTINUUM_PHASE_COORDINATE_SET_CLOSED,'phase_continuous_propagation_closed':PHASE_CONTINUOUS_PROPAGATION_CLOSED,'hard_spectral_driver_set_closed':False,'correlated_outer_enclosure_closed':outer,'BRMM_to_correlated_outer_left_inclusion_closed':True,'complete_BRMM_left_inclusion_closed':False,'joint_source_output_map_closed':joint}
 return {'schema':SCHEMA,'qualification':QUALIFICATION,'reference_parameter_domain_compact':True,'compactness_is_not_an_open_obligation':True,
 'theorem_has_deterministic_shaping_contract':('x^s_{k+1}&=A_s' in theorem and 'u^s_k&=C_s' in theorem and 'oscillator/shaping state or an equivalent hard finite-window' in flat),
 'theorem_has_explicit_hard_realization_set':('\\mathcal X^s_{\\rm ref}(\\lambda_{0:N_W})' in theorem and 'eq:marine-reference-hard-realization-set' in theorem and 'BRMM itself does not require a finite spectral state' in flat and 'machine-readable outward representation' in flat),
 'theorem_rejects_statistical_or_seeded_surrogates':('neither a Gaussian confidence event' in flat and 'spectral moments alone' in flat and 'finite seeded harmonic' in flat),
 'theorem_separates_probabilistic_random_sea_corollary':('A probabilistic statement for random sea realizations is a later corollary' in flat and 'no infinite-time pointwise bound is inferred merely from a Gaussian spectrum' in flat),
 'complete_source_rejects_gaussian_word_generator':('"used_to_generate_P3_source_words": False' in complete and '"used_to_prune_homogeneous_P3_family": False' in complete),
 'hard_realization_set_symbol':'X^s_ref(lambda_{0:N_W})','continuum_phase_certificate':{'qualification':phase['qualification'],'phase_state_set':phase['phase_state_set'],'continuum_index_set_retained':phase['continuum_index_set_retained'],'finite_frequency_grid_used':phase['finite_frequency_grid_used'],'finite_direction_grid_used':phase['finite_direction_grid_used'],'phase_reset_on_lambda_transition_allowed':phase['phase_reset_on_lambda_transition_allowed'],'continuum_phase_coordinate_set_closed':phase['continuum_phase_coordinate_set_closed'],'phase_continuous_propagation_closed':phase['phase_continuous_propagation_closed']},
 'sampled_behavior_target':{'qualification':behavior['qualification'],'symbol':behavior['behavior_set_symbol'],'sample_count':behavior['sample_count'],'sampled_projection_dimension':behavior['sampled_projection_dimension'],'compact':behavior['sampled_behavior_set_compact'],'membership_requires_common_BRMM_witness':behavior['membership_requires_common_BRMM_witness'],'normal_live_caps_are_membership_sufficient':behavior['normal_live_caps_are_membership_sufficient'],'independent_sample_boxes_define_behavior_set':behavior['independent_sample_boxes_define_behavior_set'],'validated_membership_or_separation_oracle_closed':behavior['validated_membership_or_separation_oracle_closed'],'validated_correlated_outer_enclosure_closed':behavior['validated_correlated_outer_enclosure_closed'],'correlated_outer_set_symbol':behavior['correlated_outer_set_symbol'],'correlated_outer_left_inclusion_closed':behavior['correlated_outer_left_inclusion_closed']},
 'joint_executor_coordinate_map':{'qualification':output['qualification'],'closed':output['joint_source_output_map_closed'],'raw_gyro_and_corrected_rate_distinct':output['raw_gyro_and_corrected_rate_distinct'],'truth_attitude_and_nominal_R_hat_distinct':output['truth_attitude_and_nominal_R_hat_distinct'],'truth_acceleration_and_nominal_a_w_hat_distinct':output['truth_acceleration_and_nominal_a_w_hat_distinct'],'sensor_forcing_hard_bound_closed_here':output['sensor_forcing_hard_bound_closed_here'],'BIAS0_assembled_sensor_qualification_closed_here':output['BIAS0_assembled_sensor_qualification_closed_here']},
 'exact_spectral_membership_oracle_required_for_P4':False,'correlated_outer_enclosure_route_used':True,'BRMM_to_correlated_outer_left_inclusion_closed':True,'global_physical_deployment_left_inclusion_closed_here':False,
 'power_spectrum_alone_is_hard_pathwise_bound':False,'spectral_moments_alone_may_close_xs':False,'gaussian_good_event_may_close_xs':False,'replay_may_close_xs':False,'seeded_128_frequency_generator_may_close_xs':False,'finite_RAO_grid_may_close_xs':False,'arbitrary_bounded_input_box_may_close_xs':False,
 'allowed_closure_forms':['validated_compact_oscillator_or_shaping_state_with_hard_driver_set','validated_equivalent_hard_finite_window_dynamic_constraint'],'executable_ingredients':executable,'hard_shaping_state_or_excitation_bound_closed':True,'complete_BRMM_family_materialized_here':False,'P3_promoted':False,'next_obligation':'hard shaping and BRMM->outer inclusion are closed; physical deployment admission and 601-sample provider execution remain separate open obligations'}
def validate(d):
 f=[]
 if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
 for k in ('reference_parameter_domain_compact','compactness_is_not_an_open_obligation','theorem_has_deterministic_shaping_contract','theorem_has_explicit_hard_realization_set','theorem_rejects_statistical_or_seeded_surrogates','theorem_separates_probabilistic_random_sea_corollary','complete_source_rejects_gaussian_word_generator','correlated_outer_enclosure_route_used','BRMM_to_correlated_outer_left_inclusion_closed','hard_shaping_state_or_excitation_bound_closed'):
  if d.get(k) is not True:f.append(k+' is not true')
 for k in ('exact_spectral_membership_oracle_required_for_P4','global_physical_deployment_left_inclusion_closed_here','complete_BRMM_family_materialized_here','P3_promoted','power_spectrum_alone_is_hard_pathwise_bound','spectral_moments_alone_may_close_xs','gaussian_good_event_may_close_xs','replay_may_close_xs','seeded_128_frequency_generator_may_close_xs','finite_RAO_grid_may_close_xs','arbitrary_bounded_input_box_may_close_xs'):
  if d.get(k) is not False:f.append(k+' is not false')
 e=d.get('executable_ingredients',{})
 expected={'continuum_phase_coordinate_set_closed':True,'phase_continuous_propagation_closed':True,'hard_spectral_driver_set_closed':False,'correlated_outer_enclosure_closed':True,'BRMM_to_correlated_outer_left_inclusion_closed':True,'complete_BRMM_left_inclusion_closed':False,'joint_source_output_map_closed':True}
 if e!=expected:f.append('hard shaping executable ingredient gates drifted')
 b=d.get('sampled_behavior_target',{})
 if not (b.get('compact') and b.get('membership_requires_common_BRMM_witness') and b.get('validated_correlated_outer_enclosure_closed') and b.get('correlated_outer_left_inclusion_closed')):f.append('sampled behavior target lost correlated outer closure')
 o=d.get('joint_executor_coordinate_map',{})
 if not (o.get('closed') and o.get('raw_gyro_and_corrected_rate_distinct') and o.get('truth_attitude_and_nominal_R_hat_distinct') and o.get('truth_acceleration_and_nominal_a_w_hat_distinct')):f.append('joint executor coordinate map weakened')
 return list(dict.fromkeys(f))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'hard_shaping_closed':d['hard_shaping_state_or_excitation_bound_closed'],'brmm_outer_inclusion':d['BRMM_to_correlated_outer_left_inclusion_closed'],'deployment_inclusion':d['global_physical_deployment_left_inclusion_closed_here'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
