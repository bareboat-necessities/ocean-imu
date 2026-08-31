.PHONY: build publication-sync evidence-contract test clean

build:
	python3 -m py_compile \
		../../tools/ou_validation.py \
		../../tools/ou_robustness.py \
		../../tools/ou_roundtrip_transition.py \
		../../tools/ou3_exact_replay.py \
		../../tools/ou3_numerical_certificate.py \
		../../tools/ou3_information_certificate.py \
		../../tools/ou3_information_completion.py \
		../../tools/ou3_information_enclosure_contract.py \
		../../tools/ou3_source_noise_certificate.py \
		../../tools/ou3_neighborhood_diagnostic.py \
		../../tools/ou3_neighborhood_radius_search.py \
		../../tools/ou3_certificate_completion.py \
		../../tools/ou3_validate_enclosure.py \
		../../tools/ou3_hybrid_aw_sync_proof.py \
		../../tools/ou3_source_word_theorem_contract.py \
		../../tools/ou3_translational_spread_uco.py \
		../../tools/ou3_interval_linear_algebra.py \
		../../tools/ou3_validated_kalman_interval.py \
		../../tools/ou3_source_reachable_matrix_p3.py \
		../../tools/ou3_source_reachable_matrix_p3_factored.py \
		../../tools/ou3_source_reachable_matrix_p3_direct.py \
		../../tools/ou3_p3_word_algebra.py \
		../../tools/ou3_explicit_information_word_certificate.py \
		../../tools/ou3_p4_group_algebra.py \
		../../tools/ou3_p4_exact_word_map.py \
		../../tools/ou3_p4_node_metrics.py \
		../../tools/ou3_p3_word_noise_accumulation.py \
		../../tools/ou3_p4_metric_defect_transport.py \
		../../tools/ou3_p4_nonlinear_word_certificate.py \
		../../tools/ou3_p4_p5_route_ceiling_certificate.py \
		../../tools/ou3_p5_startup_capture_certificate.py \
		../../tools/ou3_p5_outer_h_word_certificate.py \
		../../tools/ou3_p5_go_live_covariance_stage.py \
		../../tools/ou3_p5_first_s_gain_certificate.py \
		../../tools/ou3_p5_first_s_state_prefix_certificate.py \
		../../tools/ou3_p5_heading_handoff_contract.py \
		../../tools/ou3_p5_large_angle_sector_certificate.py \
		../../tools/ou3_p5_outer_information_geometry.py \
		../../tools/ou3_p5_exact_correction_transport.py \
		../../tools/ou3_p5_first_s_exact_prefix.py \
		../../tools/ou3_p5_cayley_eta_geometry.py \
		../../tools/ou3_p5_effective_vector_input.py \
		../../tools/ou3_p5_signed_cayley_cell.py \
		../../tools/ou3_p5_complete_word_transport.py \
		../../tools/ou3_p5_yaw_quotient_word_certificate.py \
		../../tools/ou3_p5_gravity_quotient_certificate.py \
		../../tools/ou3_p5_outer_h_bridge_certificate.py \
		../../tools/ou3_implementation_proof_manifest.py \
		../../tools/ou3_startup_stability_certificate.py \
		../../tools/ou3_implementation_word_language.py \
		../../tools/ou3_implementation_stability_gate.py \
		../../tools/ou3_result_contract.py \
		../../tools/ou_evidence_provenance.py \
		../../tools/ou_evidence_contract.py \
		../../tools/ou_replay_fingerprint.py \
		../../tools/ou_publication_sync.py \
		../../tools/rotate_record_heading.py \
		test_ou_validation.py \
		test_ou_robustness.py \
		test_ou_roundtrip_transition.py \
		test_ou3_numerical_certificate.py \
		test_ou3_information_certificate.py \
		test_ou3_information_completion.py \
		test_ou3_certificate_completion.py \
		test_ou3_source_noise_certificate.py \
		test_ou3_neighborhood_diagnostic.py \
		test_ou3_neighborhood_radius_search.py \
		test_ou3_validate_enclosure.py \
		test_ou3_hybrid_aw_sync_proof.py \
		test_ou3_source_word_theorem_contract.py \
		test_ou3_translational_uco_ucc.py \
		test_ou3_validated_kalman_interval.py \
		test_ou3_p3_word_algebra.py \
		test_ou3_explicit_information_word_certificate.py \
		test_ou3_p4_group_algebra.py \
		test_ou3_p4_word_map_and_metrics.py \
		test_ou3_p3_word_noise_accumulation.py \
		test_ou3_p4_metric_defect_transport.py \
		test_ou3_p4_nonlinear_word_certificate.py \
		test_ou3_p4_p5_route_ceiling.py \
		test_ou3_p5_startup_capture_certificate.py \
		test_ou3_p5_outer_h_word_certificate.py \
		test_ou3_p5_go_live_covariance_stage.py \
		test_ou3_p5_first_s_gain_certificate.py \
		test_ou3_p5_first_s_state_prefix_certificate.py \
		test_ou3_p5_heading_handoff_contract.py \
		test_ou3_p5_large_angle_sector_certificate.py \
		test_ou3_p5_outer_information_geometry.py \
		test_ou3_p5_exact_correction_transport.py \
		test_ou3_p5_first_s_exact_prefix.py \
		test_ou3_p5_cayley_eta_geometry.py \
		test_ou3_p5_effective_vector_input.py \
		test_ou3_p5_signed_cayley_cell.py \
		test_ou3_p5_complete_word_transport.py \
		test_ou3_p5_yaw_quotient_word_certificate.py \
		test_ou3_p5_gravity_quotient_certificate.py \
		test_ou3_p5_outer_h_bridge_certificate.py \
		test_ou3_startup_stability_certificate.py \
		test_ou3_implementation_word_language.py \
		test_ou3_implementation_stability_gate.py \
		test_ou3_result_determinism.py \
		test_ou_article_ablations.py \
		test_ou_legacy_law_cleanup.py \
		test_ou_paper_code_parity.py \
		test_ou_evidence_contract.py \
		test_ou_replay_fingerprint.py \
		test_ou_semiglobal_stability.py \
		test_ou_stability_phase_a.py \
		test_ou_stability_phase_b.py \
		test_ou_stability_phase_c.py \
		test_ou_stability_widening.py \
		test_ou_stability_large_angle_funnel.py \
		test_ou_stability_source_path_widening.py \
		test_ou_stochastic_concentration.py \
		test_source_archive_provenance.py \
		test_record_conventions.py \
		test_publication_references.py \
		test_deterministic_baseline_contract.py \
		test_readme_result_gallery.py \
		test_workflow_contract.py \
		test_reorg_compat_overrides.py \
		test_zz_publication_contract.py

publication-sync: build
	python3 ../../tools/ou_publication_sync.py \
		--validation-dir ../../reports/results/ou_validation \
		--direction-results ../../doc/kalman_ou_iii/w3d-wave-direction-results.tex-part \
		--deterministic-results ../../doc/kalman_ou_iii/w3d-sim-results-generated.tex-part

evidence-contract: publication-sync
	python3 ../../tools/ou_evidence_contract.py --auto

test: evidence-contract
	python3 -m unittest discover -v -p 'test_*.py'

clean:
	rm -rf __pycache__ ../../tools/__pycache__
