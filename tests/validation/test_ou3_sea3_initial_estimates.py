from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_initial_estimates as sea3  # noqa: E402


def test_sea3_initial_estimates_are_non_promoting_and_source_derived() -> None:
    payload = sea3.build_estimates(ROOT)

    assert payload["schema_version"] == "OU3_SEA3_INITIAL_ESTIMATES_V1"
    assert payload["proof_status"] == "exploratory_non_promoting"
    assert payload["certificate_promoted"] is False
    assert payload["trajectory_replay_used"] is False
    assert payload["source_generated"] is True
    assert payload["m_max"] == 3
    assert payload["p2_relationship"]["current_interface"] == (
        "OU3_P2_CORRELATED_STAGE_TRANSFER_V1"
    )
    assert payload["p2_relationship"]["intent"] == (
        "sound_refinement_not_replacement"
    )


def test_sea3_shipping_clamps_match_tightened_live_envelope() -> None:
    payload = sea3.build_estimates(ROOT)
    clamps = payload["shipping_clamps"]

    assert clamps["imu_dt_s"] == pytest.approx(0.005)
    assert clamps["tune_frequency_hz"] == pytest.approx(
        {"min": 0.03, "max": 1.2}
    )
    assert clamps["tau_aw_s"]["max"] == pytest.approx(12.0)
    assert clamps["sigma_aw_mps2"]["max"] == pytest.approx(4.0)
    assert clamps["r_s"] == pytest.approx({"min": 0.15, "max": 100.0})
    assert clamps["s_pseudo_update_period_s"]["max"] == pytest.approx(0.15)
    assert clamps["dynamic_ema_horizon_s"]["max"] == pytest.approx(35.0)
    assert clamps["normal_live_accelerometer_rejection_allowed"] is False


def test_sea3_screening_estimates_preserve_oscillatory_time_structure() -> None:
    payload = sea3.build_estimates(ROOT)
    estimates = payload["derived_screening_estimates"]

    assert estimates["single_component_phase_advance_per_imu_sample_deg"] == (
        pytest.approx({"min": 0.054, "max": 2.16})
    )
    assert estimates["single_component_phase_advance_over_max_s_gap_deg"] == (
        pytest.approx({"min": 1.62, "max": 64.8})
    )
    assert estimates["max_s_gap_imu_samples"] == 30
    assert estimates["wave_period_estimator_screening_proxy"][
        "promotion_use"
    ] is False

    dimensions = {
        row["oscillator_pairs_per_mode"]: row["sea_shaping_state_dimension"]
        for row in estimates["oscillator_shaping_state_examples"]
    }
    assert dimensions == {2: 12, 4: 24, 6: 36}


def test_committed_sea3_estimate_artifact_is_current() -> None:
    expected = sea3.build_estimates(ROOT)
    artifact = json.loads(
        (ROOT / "tools/ou3_sea3_initial_estimates.json").read_text(
            encoding="utf-8"
        )
    )
    assert artifact == expected
