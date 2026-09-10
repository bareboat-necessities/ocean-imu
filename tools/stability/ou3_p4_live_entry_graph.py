#!/usr/bin/env python3
"""Exact fresh-Live entry/source relation, not a covariance confidence set.

The S integration constant is a *shared source coordinate*. It is not an
independent initial error. Re-anchoring S at the handoff is an exact reference
change; it changes neither a shipping state nor any sensor sample. Position is
NOT re-anchored. Full P4 endpoint/prefix/storage obligations remain separate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[2]
S_OFFSET = 12


def _vector(v: Sequence, n: int = 3) -> tuple[F, ...]:
    if len(v) != n:
        raise ValueError(f"expected {n} coordinates")
    return tuple(F(x) for x in v)


def reanchor_S(S: Sequence, S_at_handoff: Sequence) -> tuple[F, ...]:
    """Exact proof-reference operation; never write this into the estimator."""
    return tuple(x-y for x, y in zip(_vector(S), _vector(S_at_handoff)))


def S_residual(error_S: Sequence, physical_S: Sequence) -> tuple[F, ...]:
    return tuple(x-y for x, y in zip(_vector(error_S), _vector(physical_S)))


def shifted_S_energy(energy, integral_S: Sequence, origin: Sequence, duration) -> F:
    """Integral |S-c|^2 = integral |S|^2 - 2 c.integral(S) + T |c|^2.

    The signed cross term must come from the SAME history. This is not an
    authorization to choose an independent S/noise energy or discard origin.
    """
    c = _vector(origin)
    j = _vector(integral_S)
    t, energy = F(duration), F(energy)
    if t < 0 or energy < 0:
        raise ValueError("negative duration or energy")
    value = energy - 2*sum(x*y for x, y in zip(c, j)) + t*sum(x*x for x in c)
    if value < 0:
        raise ValueError("energy, integral and origin do not describe one history")
    return value


def physical_jet(p0: Sequence, v0: Sequence, h,
                 acceleration_moments: Sequence[Sequence]) -> dict:
    """Exact continuous-history jet, with no ZOH or finite harmonic substitution.

    moments[j] = integral_0^h (h-s)^j/j! a(s) ds, j=0,1,2.
    These are one history's moments, NOT three independent forcing ports.
    The returned graph is usable in a signed source ledger before enclosure.
    """
    if len(acceleration_moments) != 3:
        raise ValueError("three jointly generated acceleration moments required")
    h = F(h)
    if h < 0:
        raise ValueError("negative duration")
    p, v = _vector(p0), _vector(v0)
    j0, j1, j2 = (_vector(j) for j in acceleration_moments)
    return {
        "v": tuple(v[i]+j0[i] for i in range(3)),
        "p": tuple(p[i]+h*v[i]+j1[i] for i in range(3)),
        "S_from_handoff": tuple(h*p[i]+h*h*v[i]/2+j2[i] for i in range(3)),
    }


def prefix_S_bound(h, p_at_entry_bound, v_at_entry_bound, acceleration_bound,
                   *, uniform_position_bound=None, uniform_S_bound=None) -> F:
    """A consequence of the jet, not a replacement for its correlated graph.

    With a globally bounded original S, the re-anchored bound is 2*S_m.
    The working-domain charge is therefore retained, not erased at entry.
    """
    h, p, v, a = map(F, (h, p_at_entry_bound, v_at_entry_bound, acceleration_bound))
    if min(h, p, v, a) < 0:
        raise ValueError("nonnegative bounds required")
    bounds = [h*p+h*h*v/2+a*h**3/6]
    if uniform_position_bound is not None:
        pm = F(uniform_position_bound)
        if pm < 0: raise ValueError("negative position bound")
        bounds.append(h*pm)
    if uniform_S_bound is not None:
        sm = F(uniform_S_bound)
        if sm < 0: raise ValueError("negative S bound")
        bounds.append(2*sm)
    return min(bounds)


def center_prefix_map(error_map: Sequence[Sequence], origin_map: Sequence[Sequence]):
    """Exact row operation for an augmented endpoint/every-prefix master.

    Keep every graph column: deleting an origin column is legal only AFTER
    both the error and the physical S forcing have been transformed. This
    function deliberately neither deletes columns nor issues a P4 flag.
    """
    if len(error_map) not in (18, 21) or len(origin_map) != 3:
        raise ValueError("expected H18/A21 map and three origin rows")
    width = len(error_map[0])
    if width == 0 or any(len(r) != width for r in (*error_map, *origin_map)):
        raise ValueError("ragged/empty augmented graph")
    out = [list(map(F, row)) for row in error_map]
    for i in range(3):
        out[S_OFFSET+i] = [x-F(y) for x, y in zip(out[S_OFFSET+i], origin_map[i])]
    return out


def cancellation_certificate() -> dict:
    """Coefficientwise polynomial identity for ARBITRARY same-event gain.

    (I-K H_S) E_S + K = E_S.  Entries of K are formal independent symbols
    ONLY for checking this universal algebraic identity, never proof boxes
    or a substitute Riccati family. Non-S event H E_S=0 gives the analogous
    identity. Finite attitude correction/projection sees an unchanged
    residual, hence exactly the same nominal/reset/arithmetic execution.
    """
    # Sparse linear polynomials in formal K[i,j], with None the constant.
    # Actually form the matrix products, so a sign/selector/offset change is
    # detected rather than pre-asserting the expected cancellation.
    def add(a, b, scale=F(1)):
        out=dict(a)
        for key,value in b.items():
            out[key]=out.get(key,F(0))+scale*value
        return {key:value for key,value in out.items() if value}
    checked=0
    for n in (18,21):
        E=[[F(int(i==S_OFFSET+j)) for j in range(3)] for i in range(n)]
        H=[[F(int(k==S_OFFSET+j)) for k in range(n)] for j in range(3)]
        K=[[{(i,j):F(1)} for j in range(3)] for i in range(n)]
        for i in range(n):
            for j in range(3):
                value={None:E[i][j]}
                for ell in range(3):
                    he=sum(H[ell][k]*E[k][j] for k in range(n))
                    value=add(value,K[i][ell],-he)
                value=add(value,K[i][j])
                value=add(value,{None:E[i][j]},F(-1))
                assert not value
                checked+=1
        centered=center_prefix_map(E,[[int(i==j) for j in range(3)] for i in range(3)])
        assert all(x==0 for row in centered for x in row)
    return {
        "identity": "(I-K*H_S)*E_S + K = E_S",
        "formal_gain_coefficients_checked": checked,
        "coefficient_residuals_exactly_zero": True,
        "scope": "all completed events; any same-history gain and finite reset/projection",
        "source_shift": "S_true -> S_true+c; e_S -> e_S+c",
        "nominal_residual_and_binary32_execution": "unchanged bit-for-bit",
        "does_not_prove_remaining_motion_contraction": True,
    }


def quiet_plateau_witness() -> dict:
    """C2 periodic physical motion, used to falsify the entrance model ONLY.

    In each coordinate: +A on [0,D], A*(1-2*s(u)) on [D,2D],
    -A on [2D,3D], -A*(1-2*s(u)) on [3D,4D]; u in [0,1],
    s=10u^3-15u^4+6u^5. Half-period antisymmetry bounds every primitive.
    |s'|<=15/8 and |s''|<=6; sqrt(3)<7/4 encloses Euclidean norms.
    """
    A, D, norm3 = F(21,10), F(200), F(7,4)
    vm = norm3*2*A*F(15,8)/D
    am = norm3*2*A*6/(D*D)
    pm = norm3*A
    sm = 4*D*pm  # conservative global primitive bound from a zero-mean cycle
    tr, eq = F(10), F(9,1000)
    assert tr*am*am < eq and am < 4 and 2*vm < 2
    return {
        "position_each_axis_amplitude_m": str(A), "plateau_duration_s": str(D),
        "period_s": str(4*D), "velocity_norm_upper_mps": str(vm),
        "acceleration_norm_upper_mps2": str(am), "position_norm_upper_m": str(pm),
        "S_norm_global_upper_m_s": str(sm), "window_s": str(tr),
        "all_window_AC_energy_upper": str(tr*am*am), "quiet_threshold": str(eq),
        "all_continuous_window_starts_covered": True,
        "Q_only_no_minimum_excitation_needed": True,
        "vertical_Hs_upper_via_4_std_le_4A_m": str(4*A),
        "mean_over_every_full_cycle": "p=0, a=0",
        "session_origin_at_150_seconds": "S=(315,315,315) m*s; norm>300",
        "handoff_origin": "S_from_handoff=0 at entry; p is NOT zeroed",
        "role": "analytic BRMM-Q witness conditional on displayed budgets; not deployment qualification or a replacement source family",
    }


def build() -> dict:
    files = ["src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
             "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h",
             "src/kalman_common/SeaStateFusionFilterCommon.h"]
    return {
        "qualification": "OU3_EXACT_FRESH_LIVE_ENTRY_SOURCE_GRAPH_V1",
        "shipping_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        "boundary": "fresh outer begin -> first goLive, before next completed Live event",
        "nominal_v_p_S_aw_bg_ba": "exactly zero carried from construction",
        "goLive_resets_linear_values": False,
        "previously_driven_inner_filter_covered_by_zero_claim": False,
        "session_origin_errors": {"v":"v_true(t_L)","p":"p_true(t_L)","S":"S_true(t_L)","aw":"a_true(t_L)","bg":"bg_true(t_L)","ba":"beta_true(t_L)"},
        "canonical_reference_for_fresh_entry": "S_L(t)=S_true(t)-S_true(t_L); physical p is unchanged",
        "canonical_S_entry_error": 0,
        "independent_session_S_entry_error_ball_admitted": False,
        "position_origin_at_handoff_established_by_shipping": False,
        "physical_velocity_position_acceleration_bounds": "same BRMM history; no startup Cartesian replacement",
        "origin_scope": "one anchor for the whole Live continuation, NOT re-zeroed at every word",
        "fresh_covariance_diagonal_v_p_S": [1,400,2500],
        "fresh_covariance_cross_v_p_S": "zero at handoff; subsequent propagation retained",
        "aw_covariance_at_handoff": "same active Sigma_aw, all aw cross-covariances cleared",
        "BA_release": "state preserved; existing covariance retained except diagonal floor",
        "tuner_state": "learned history retained and committed, not reinitialized",
        "pseudo_scheduler": "zero progress while MEKF held; period retargeted without fictitious tick",
        "timeout_requires_usable_WPE": False,
        "cancellation": cancellation_certificate(),
        "quiet_witness": quiet_plateau_witness(),
        "working_S_bound": "min(h*P_m, h*|p_L|+h^2*|v_L|/2+A*h^3/6, 2*S_m), where admitted",
        "source_forcing_energy": "E_shift=E-2*c.integral(S)+T*|c|^2; shared-history signed charge required",
        "raw_error_recovery": "|e_original| <= |e_L| + |S_true(t_L)|",
        "uniform_original_origin_ultimate_bound_requires_uniform_origin_bound": True,
        "independent_S_gauge_extent": "no attitude-chart restriction; original-error bound still charges origin",
        "source_uniform_reachable_other_entry_coordinates_closed": False,
        "source_uniform_endpoint_and_every_prefix_LDLT_closed": False,
        "P4_MOTION_PASS": False, "P4_PASS": False, "P5_MAY_START": False,
    }


def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--native-log",type=Path)
    a=p.parse_args();d=build()
    if a.native_log:
        lines=a.native_log.read_text().splitlines()
        if lines.count("LIVE_ENTRY_AUDIT_PASS=true")!=1:
            raise ValueError("native audit did not pass")
        records=[]
        for line in lines:
            if line.startswith("LIVE_ENTRY "):
                fields=line.split();record={"magnetometer":fields[1]=="mag"}
                for field in fields[2:]:
                    key,value=field.split("=",1);record[key]=float(value)
                records.append(record)
        if len(records)!=4 or len({(r["magnetometer"],r["zero_proxy"]) for r in records})!=4:
            raise ValueError("missing native startup branches")
        d["native_observations"]=records
        d["native_observations_are_target_arithmetic_qualification"]=False
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print("EXACT_SHARED_S_ORIGIN_CANCELLATION_PASS=true P4_PASS=false P5_MAY_START=false")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
