"""Every sample-aligned 10 s window of the eight unchanged reference records.

Outward arithmetic bounds the binary64 CSV-point, zero-order-hold statistics.
It does not certify intersample physical motion, Normal-Live membership,
frontend/tuner attachment, a source family, or P3/P4. Lobe scores are lower
bounds from explicit disjoint interval witnesses; no witness is UNRESOLVED,
never a proof that no opposing intervals exist.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile

import numpy as np

FAMILIES = ("jonswap", "pmstokes")
HEIGHTS = (.27, 1.5, 4., 8.5)
EPS = np.finfo(float).eps


def down(x):
    return np.nextafter(x, -np.inf)


def up(x):
    return np.nextafter(x, np.inf)


def point(x):
    return np.asarray(x), np.asarray(x)


def add(a, b):
    return down(a[0]+b[0]), up(a[1]+b[1])


def sub(a, b):
    return down(a[0]-b[1]), up(a[1]-b[0])


def mul(a, b):
    values = [a[i]*b[j] for i in (0, 1) for j in (0, 1)]
    return down(np.minimum.reduce(values)), up(np.maximum.reduce(values))


def square(a):
    lo = np.where((a[0] <= 0) & (a[1] >= 0), 0., np.minimum(a[0]**2, a[1]**2))
    return np.maximum(0., down(lo)), up(np.maximum(a[0]**2, a[1]**2))


def norm2(a):
    q = square(a)
    out = point(np.zeros(q[0].shape[:-1]))
    for i in range(q[0].shape[-1]):
        out = add(out, (q[0][..., i], q[1][..., i]))
    return np.maximum(0., out[0]), out[1]


def sqrt(a):
    return np.maximum(0., down(np.sqrt(np.maximum(0., a[0])))), up(np.sqrt(np.maximum(0., a[1])))


def cumulative(a):
    """Any-order n-addition error bound, with outward absolute-sum majorant."""
    n = len(a[0])
    if n*EPS >= .01:
        raise ValueError("record too long for the cumulative rounding bound")
    g = up((n*EPS)/(1-n*EPS))
    result = []
    for endpoint, sign in zip(a, (-1, 1)):
        approx = np.cumsum(endpoint, axis=0)
        absolute = up(np.sum(np.abs(endpoint), axis=0)/down(1-g))
        error = up(g*absolute)
        edge = down(approx-error) if sign < 0 else up(approx+error)
        result.append(np.concatenate((np.zeros_like(edge[:1]), edge)))
    return tuple(result)


def intervals(c, start, stop, count):
    return sub((c[0][stop:stop+count], c[1][stop:stop+count]),
               (c[0][start:start+count], c[1][start:start+count]))


def dot(a, b):
    out = point(np.zeros(a[0].shape[:-1]))
    for i in range(a[0].shape[-1]):
        out = add(out, mul((a[0][..., i], a[1][..., i]), (b[0][..., i], b[1][..., i])))
    return out


def lobe_candidates(n, sample_hz):
    # Prefix/suffix splits plus local adjacent intervals avoid integer-period
    # cancellation in a whole-window-only search. These are witness candidates,
    # not an exhaustive optimization over all interval pairs.
    choices = set()
    for split in range(sample_hz, n, sample_hz):
        choices.add((0, split, split, n))
    hop = sample_hz//4
    for seconds in (.25, .5, 1., 2., 3., 4.):
        width = round(seconds*sample_hz)
        for start in range(0, n-2*width+1, hop):
            choices.add((start, start+width, start+width, start+2*width))
    return sorted(choices)


def windows(acceleration, *, sample_hz=200, duration_s=10):
    a = np.asarray(acceleration, dtype=float)
    n = int(sample_hz*duration_s)
    if a.ndim != 2 or a.shape[1] != 3 or len(a) < n or not np.all(np.isfinite(a)):
        raise ValueError("finite 3D acceleration and a complete window required")
    if sample_hz < 4 or sample_hz % 4 or n != sample_hz*duration_s:
        raise ValueError("integer window and quarter-second sample boundaries required")
    count = len(a)-n+1
    h = down(1./sample_hz), up(1./sample_hz)
    impulses = cumulative(mul(point(a), h))
    total = intervals(impulses, 0, n, count)
    raw_energy = intervals(cumulative(mul(norm2(point(a)), h)), 0, n, count)
    centered = sub(raw_energy, mul(norm2(total), (down(1./duration_s), up(1./duration_s))))
    centered = np.maximum(0., centered[0]), np.maximum(0., centered[1])
    jumps = sqrt(norm2(sub(point(a[1:]), point(a[:-1]))))
    tv = intervals(cumulative(jumps), 0, n-1, count)
    choices = lobe_candidates(n, sample_hz)
    best = np.zeros(count)
    witness = np.full(count, -1, dtype=np.int16)
    for i, (s1, t1, s2, t2) in enumerate(choices):
        first, second = intervals(impulses, s1, t1, count), intervals(impulses, s2, t2, count)
        # chi=1 here is only an audit setting. The stored intervals prove
        # norms>=J and dot<=-J^2 whenever this outward lower score is positive.
        score = np.maximum(0., np.minimum.reduce((norm2(first)[0], norm2(second)[0], -dot(first, second)[1])))
        improve = score > best
        witness[improve] = i
        best = np.maximum(best, score)
    return {"energy_lo": centered[0], "energy_hi": centered[1],
            "impulse_lo": total[0], "impulse_hi": total[1],
            "impulse_norm_lo": sqrt(norm2(total))[0], "impulse_norm_hi": sqrt(norm2(total))[1],
            "tv_lo": np.maximum(0., tv[0]), "tv_hi": tv[1],
            "lobe_J_lower": sqrt(point(best))[0], "lobe_witness_index": witness,
            "lobe_candidate_sample_offsets": np.asarray(choices),
            "sample_start_index": np.arange(count), "sample_hz": np.array(sample_hz),
            "window_samples": np.array(n)}


def stats(w, threshold, duration_s=10.):
    eq = duration_s*threshold**2
    eq_interval = mul(square(point(threshold)), point(duration_s))
    quiet, oscillatory = w["energy_hi"] <= eq_interval[0], w["energy_lo"] > eq_interval[1]
    unknown = ~(quiet | oscillatory)
    q_violate = quiet & (w["impulse_norm_lo"] > 2.)
    o_violate = oscillatory & (w["impulse_norm_lo"] > 2.)
    no_witness = oscillatory & (w["lobe_J_lower"] <= 0)
    def first(mask):
        indices = np.flatnonzero(mask)
        return int(indices[0]) if len(indices) else None
    return {"quiet_rms_threshold_mps2": threshold, "E_q": eq,
            "Q_windows": int(quiet.sum()), "O_windows": int(oscillatory.sum()),
            "classification_rounding_ambiguous": int(unknown.sum()),
            "Q_impulse_cap_2_violations": int(q_violate.sum()),
            "O_impulse_cap_2_violations": int(o_violate.sum()),
            "first_Q_impulse_violation_sample": first(q_violate),
            "first_O_impulse_violation_sample": first(o_violate),
            "O_lobe_search_unresolved": int(no_witness.sum()),
            "first_O_lobe_unresolved_sample": first(no_witness),
            "O_energy_min_lower": float(w["energy_lo"][oscillatory].min()) if oscillatory.any() else None,
            "O_lobe_J_min_witness_lower": float(w["lobe_J_lower"][oscillatory].min()) if oscillatory.any() else None,
            "lobe_search_failure_proves_nonexistence": False}


def audit_case(data, name, output):
    columns = np.genfromtxt(io.BytesIO(data), delimiter=",", names=True, dtype=float)
    if columns.ndim != 1 or len(columns) == 0:
        raise ValueError("empty reference record")
    a = np.column_stack([columns["acc_"+axis] for axis in "xyz"])
    gyro = np.column_stack([columns["gyro_"+axis] for axis in "xyz"])
    velocity = np.column_stack([columns["vel_"+axis] for axis in "xyz"])
    if not np.all(np.isfinite(gyro)) or not np.all(np.isfinite(velocity)):
        raise ValueError("nonfinite body rate or velocity")
    w = windows(a)
    np.savez_compressed(output/(Path(name).stem+".windows.npz"), **w)
    an, gn = sqrt(norm2(point(a))), sqrt(norm2(point(gyro)))
    radians_cap = np.pi/6
    report = {"record": name, "csv_sha256": hashlib.sha256(data).hexdigest(),
              "samples": len(a), "sliding_windows": len(w["energy_lo"]),
              "acceleration_norm_max_upper_mps2": float(an[1].max()),
              "body_rate_norm_max_upper_rad_s": float(gn[1].max()),
              "body_rate_norm_max_deg_s_diagnostic": float(gn[1].max()*180/np.pi),
              "acceleration_cap_4_definite_exceedances": int((an[0] > 4.).sum()),
              "body_rate_cap_30_definite_exceedances": int((gn[0] > up(radians_cap)).sum()),
              "window_impulse_norm_max_upper_mps": float(w["impulse_norm_hi"].max()),
              "AC_energy_min_lower": float(w["energy_lo"].min()),
              "AC_energy_max_upper": float(w["energy_hi"].max()),
              "TV_max_upper_ZOH_mps2": float(w["tv_hi"].max()),
              "thresholds": [stats(w, q) for q in (.03, .05)],
              "all_samples_audited_no_Normal_Live_selection": True,
              "complete_frontend_tuner_RS_word_admission": False}
    impulse_prefix = cumulative(mul(point(a), (down(.005), up(.005))))
    report["record_prefix_impulse_norm_max_upper_mps"] = float(sqrt(norm2(impulse_prefix))[1].max())
    report["record_velocity_norm_max_upper_mps"] = float(sqrt(norm2(point(velocity)))[1].max())
    report["infinite_horizon_zero_DC_certified_from_finite_record"] = False
    # A sampled left-rule residual is not a continuous kinematic inconsistency:
    # it includes quadrature and CSV quantization. Keep it descriptive.
    predicted_dv = np.cumsum(a[:-1]*.005, axis=0)
    report["sampled_velocity_primitive_residual_max_mps_diagnostic"] = float(
        np.linalg.norm(predicted_dv-(velocity[1:]-velocity[0]), axis=1).max())
    report["diagnostic_V_R_candidate_with_25_percent_margin"] = np.ceil(1.25*report["window_impulse_norm_max_upper_mps"]*10)/10
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {"qualification": "BRMM_EIGHT_REFERENCE_RECORD_SAMPLED_AUDIT_V1",
              "simulation_release": "oceanography-waves-lib/v1.1.3",
              "archive_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
              "clock": "200 Hz sample index; rounded CSV timestamps are not integration weights",
              "frame": "fixed geographic world CoG acceleration; constant NED rotation preserves these statistics",
              "arithmetic_scope": "outward bounds for binary64 CSV-point ZOH, sample-aligned windows",
              "TV_scope": "internal jumps of each half-open ZOH window; not an intersample physical TV upper bound",
              "continuous_motion_intersample_inclusion_certified": False,
              "positive_lobe_score_is_a_witness_not_an_optimality_certificate": True,
              "true_bias_or_sensor_noise_included_in_marine_acceleration": False,
              "physical_parameters_frozen_from_replays": False,
              "BRMM_SOURCE_UNIFORM_PASS": False, "BRMM_P3_PASS": False,
              "BRMM_P4_MOTION_PASS": False, "BRMM_P5_MOTION_MAY_START": False, "cases": []}
    with zipfile.ZipFile(args.archive) as archive:
        for family in FAMILIES:
            for height in HEIGHTS:
                expression = rf"wave_data_{family}_H{height:.3f}_L[^_]+_A[^_]+_P[^_]+\.csv"
                names = [n for n in archive.namelist() if re.fullmatch(expression, Path(n).name)]
                if len(names) != 1:
                    raise ValueError(f"expected exactly one {family} H={height} record, got {names}")
                result = audit_case(archive.read(names[0]), names[0], args.output_dir)
                report["cases"].append(result)
                (args.output_dir/"brmm-audit.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
                print("BRMM_CASE", json.dumps(result, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
