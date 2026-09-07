import importlib.util
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MOD_PATH = ROOT / "tools" / "stability" / "ou3_p4_complete_sea3_residual_sector.py"
spec = importlib.util.spec_from_file_location("ou3_p4_complete_sea3_residual_sector", MOD_PATH)
assert spec and spec.loader
SECTOR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SECTOR)


def _norm2(v):
    return sum(float(x) * float(x) for x in v)


def _close_vec(a, b, tol=2e-14):
    assert len(a) == len(b)
    assert max(abs(float(x) - float(y)) for x, y in zip(a, b)) <= tol


def test_vector_remainder_exact_factorization():
    cases = [
        ([0.0, 0.0, 0.0], [0.2, -9.7, 1.1]),
        ([0.03, -0.02, 0.01], [0.2, -9.7, 1.1]),
        ([0.4, -0.3, 0.2], [-2.0, 1.5, 0.7]),
    ]
    for c, v in cases:
        direct = SECTOR.exact_vector_eta(c, v)
        factored = SECTOR.exact_vector_eta_factorized(c, v)
        _close_vec(direct, factored)


def test_accelerometer_remainder_exact_factorization_and_no_ba_term():
    c = [0.15, -0.08, 0.04]
    f = [0.3, -9.4, 1.2]
    # Proper orthogonal body/world rotation fixture.
    a = 0.31
    R = [
        [math.cos(a), -math.sin(a), 0.0],
        [math.sin(a), math.cos(a), 0.0],
        [0.0, 0.0, 1.0],
    ]
    da = [0.4, -0.2, 0.1]
    direct = SECTOR.exact_accelerometer_eta(c, f, R, da)
    factored = SECTOR.exact_accelerometer_eta_factorized(c, f, R, da)
    _close_vec(direct, factored)


def test_vector_quadratic_sector_bound():
    c = [0.17, -0.11, 0.06]
    v = [0.3, -0.5, 0.9]
    r = math.sqrt(_norm2(c))
    lam = 3.7
    eta = SECTOR.exact_vector_eta(c, v)
    qdiag = SECTOR.sector_quadratic_diagonal(r, math.sqrt(_norm2(v)), lam, include_aw=False)
    rhs = sum(qdiag[i] * c[i] * c[i] for i in range(3))
    assert lam * _norm2(eta) <= rhs * (1.0 + 2e-13) + 1e-15


def test_accelerometer_quadratic_sector_bound():
    c = [0.12, 0.09, -0.07]
    f = [0.2, -9.6, 0.8]
    a = -0.43
    R = [
        [math.cos(a), 0.0, math.sin(a)],
        [0.0, 1.0, 0.0],
        [-math.sin(a), 0.0, math.cos(a)],
    ]
    da = [0.35, -0.18, 0.22]
    r = math.sqrt(_norm2(c))
    lam = 5.2
    eta = SECTOR.exact_accelerometer_eta(c, f, R, da)
    qdiag = SECTOR.sector_quadratic_diagonal(r, math.sqrt(_norm2(f)), lam, include_aw=True)
    x = c + da
    rhs = sum(qdiag[i] * x[i] * x[i] for i in range(6))
    assert lam * _norm2(eta) <= rhs * (1.0 + 3e-13) + 1e-15


def test_contract_metadata_stays_fail_closed():
    d = SECTOR.build()
    assert SECTOR.validate(d) == []
    assert d["canonical_source"] == "COMPLETE_SEA3_NORMAL_LIVE_WORD"
    assert d["packet_count_multiplier_used"] is False
    assert d["global_correction_radius_used"] is False
    assert d["P4_promoted_here"] is False
