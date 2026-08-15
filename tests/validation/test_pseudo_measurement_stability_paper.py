"""Publication contract for the unified pseudo-measurement stability paper."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PAPER = REPO_ROOT / "doc" / "kalman_ou_iii" / "pseudo-meas-stability.tex"


def paper_text() -> str:
    return PAPER.read_text(encoding="utf-8")


def test_unified_stability_paper_exists_and_covers_all_estimators():
    text = paper_text()
    assert "Pseudo-Measurement Anchoring and Input-to-State Stability" in text
    for token in (
        "OU--II",
        "OU--III",
        "TFG",
        "PII",
        "Bryne",
        "Fossen",
        "Johansen",
    ):
        assert token in text


def test_generalized_anchoring_and_local_iss_proofs_remain_explicit():
    text = paper_text()
    for label in (
        r"\label{thm:head-anchor}",
        r"\label{lem:stable-tail}",
        r"\label{thm:kalman-template}",
        r"\label{lem:ues-local-iss}",
        r"\label{cor:ou2}",
        r"\label{cor:ou3}",
    ):
        assert label in text
    assert "Vandermonde" in text
    assert "Hermite--Genocchi" in text
    assert "uniformly detectable" in text
    assert "quadratic remainder" in text


def test_tfg_claim_is_local_and_not_a_fake_global_inekf_theorem():
    text = paper_text()
    assert r"\label{thm:coord-ues}" in text
    assert r"\label{cor:tfg}" in text
    assert "not claimed to be group affine" in text
    assert "compact chart" in text
    assert "Jacobian-consistency" in text


def test_pii_scheduling_requires_an_explicit_rate_margin():
    text = paper_text()
    assert r"\label{thm:pii-scheduled}" in text
    assert r"\dot r(t)/r(t)" in text
    assert r"\label{eq:pii-rate-condition}" in text
    assert "optional bias" in text.lower()
    assert "small-gain" in text


def test_nlo_distinguishes_published_usges_from_scheduled_theta_extension():
    text = paper_text()
    assert r"\label{prop:nlo-iss}" in text
    assert r"\label{thm:nlo-scheduled}" in text
    assert "USGES" in text
    assert r"\dot L_\theta L_\theta^{-1}" in text
    assert r"\frac{\dot\theta}{\theta}" in text
    assert "published USGES theorem cannot simply be quoted" in text
    assert "conditional local ISS" in text


def test_virtual_zero_is_treated_as_a_disturbance_not_physical_identity():
    text = paper_text()
    assert r"\label{prop:stability-accuracy}" in text
    assert "anchoring model" in text
    assert "not a claim that the true" in text
    assert "stability" in text.lower() and "accuracy" in text.lower()
