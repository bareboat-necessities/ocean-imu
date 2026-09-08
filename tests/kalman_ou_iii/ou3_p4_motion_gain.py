"""Connected finite-coefficient motion/supply feasibility, never P4 promotion.

Exact algebraic factorizations give e_next=A(e,source)e+B(e,source)u.
This is NOT a Jacobian, a fitted secant, or a replacement shipping observer.
All 21 rows, within-sample corrected bias, and actual R_S gains are retained.
The coefficients evaluated on ONE attached word define a numerical diagnostic.
Freezing them does not certify the nonlinear/source graph that generated them.

The augmented supply test retains the complete bias trajectory in its cost.
Forcing ports are relaxed to arbitrary energy inputs for this sufficient
coefficient test, with explicit units; no independent SEA3 boxes are created.
Failure of this larger test is not a nonlinear/source-admissible counterexample.
Success still requires correlated nonlinear/source coverage and useful budgets.
"""
from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

import numpy as np

import ou3_p4_connected_motion as C

I3 = np.eye(3)
I21 = np.eye(21)
# These are coordinate units, NOT physical bounds or fitted source weights.
SUPPLY_UNITS = {
    "bias": "integral |e_b/(1 m/s^2)|^2 dt/(1 s), sample entrances",
    "latent": "sum squared increments in [v/(1 m/s), p/(1 m), S/(1 m s), a_w/(1 m/s^2)]",
    "gyro": "sum |Cayley discretization defect/(1 rad)|^2",
    "accelerometer": "sum u_acc^T actual_R_acc^-1 u_acc",
    "magnetometer": "sum u_mag^T actual_R_mag^-1 u_mag",
    "S_zero": "sum S_true^T actual_R_S^-1 S_true",
}


def sym(a):
    return (a+a.T)/2


def quat_coefficients(s):
    if s < 1e-4:
        return 1-s/8+s*s/384, .5-s/48+s*s/3840
    t = np.sqrt(s)
    return np.cos(t/2), np.sin(t/2)/t


def gyro_secant(a, b):
    """Exact V(a,b)b=Cayley(R(a+b)R(a)^T), shipping quaternion branches.

    For the small-angle branch, polynomial divided differences are evaluated
    without cancellation, including equal norms and b=0. Other branches use
    long-double divided differences; the diagnostic refuses a near-coincident
    non-polynomial pair rather than inserting a derivative approximation.
    """
    a, b = np.asarray(a, dtype=np.longdouble), np.asarray(b, dtype=np.longdouble)
    s, t = a@a, (a+b)@(a+b)
    wn, kn = quat_coefficients(s)
    wt, kt = quat_coefficients(t)
    if max(s, t) < 1e-4:
        dw, dk = -1/8+(s+t)/384, -1/48+(s+t)/3840
    else:
        if abs(t-s) <= 1e-12*max(s, t):
            raise ValueError("near-coincident non-polynomial gyro secant needs validated divided differences")
        dw, dk = (wt-wn)/(t-s), (kt-kn)/(t-s)
    denominator = wt*wn+kt*kn*((a+b)@a)
    if denominator <= 0:
        raise ValueError("relative gyro rotation left the retained Cayley chart")
    matrix = (wn*kt*I3+kt*kn*C.skew(a)
              + np.outer(a, 2*a+b)*(wn*dk-kn*dw))*(2/denominator)
    return np.asarray(matrix, dtype=float)


def measurement_factor(row, e):
    """Exact finite residual, finite reset, and zero-true-bias projection."""
    if row["true_bias"] != [0, 0, 0]:
        raise ValueError("this point capture requires its original zero true-bias root")
    residual, u, h = C.residual_graph(row, e)
    if row["kind"] in ("accelerometer", "magnetometer"):
        # (E(c)-I)f = -(I-[c]/2)^-1 [f] c, exactly.
        h[:, :3] = np.linalg.solve(I3-.5*C.skew(e[:3]), h[:, :3])
        if row["kind"] == "accelerometer":
            h[:, 15:18] = C.rotation(e[:3])@C.mat(row, "R_hat")
    k = C.mat(row, "K", 21, 3)
    d = k@residual
    w, q = quat_coefficients(d[:3]@d[:3])
    scale = 2*q/w
    denominator = 1+scale*(e[:3]@d[:3])/4
    if denominator <= 0:
        raise ValueError("finite correction left the retained Cayley chart")
    a, b = I21-k@h, -k.copy()
    reset_gain = scale*(I3+.5*C.skew(e[:3]))@k[:3]
    a[:3] = (I21[:3]-reset_gain@h)/denominator
    b[:3] = -reset_gain/denominator
    before_projection = e[18:]-d[18:]
    norm = np.linalg.norm(before_projection)
    projection = min(1., row["projection_radius"]/norm) if norm else 1.
    a[18:] *= projection
    b[18:] *= projection
    return a, b, u, C.mat(row, "R"), float(denominator)


def prediction_factor(before, row, dt):
    e, _ = C.SOURCE.error_and_covariance(before, 21)
    a = I21.copy()
    a[6:18, 6:18] = C.mat(row, "F_LL", 12)
    a[18:, 18:] *= np.exp(-dt/row["tau_b"]) if row["active"] else 1.
    step = -dt*np.asarray(row["omega_hat"])
    v = gyro_secant(step, dt*e[3:6])*dt
    bg = v@e[3:6]
    rn = C.deployed_rotation(step)
    denominator = 1-bg@(rn@e[:3])/4
    if denominator <= 0:
        raise ValueError("prediction bias/attitude product left chart")
    ac, ag = (I3+.5*C.skew(bg))@rn/denominator, v/denominator
    base = ac@e[:3]+ag@e[3:6]
    truth_defect = C.mat(row, "R_true")@(
        C.deployed_rotation(-dt*np.asarray(row["gyro_measured"]))@C.mat(before, "R_true")).T
    gyro_input = C.cayley(truth_defect)
    source_denominator = 1-gyro_input@base/4
    if source_denominator <= 0:
        raise ValueError("physical gyro forcing left chart")
    left = (I3+.5*C.skew(gyro_input))/source_denominator
    a[:3] = 0
    a[:3, :3], a[:3, 3:6] = left@ac, left@ag
    b = np.zeros((21, 15))
    b[6:18, :12], b[:3, 12:] = np.eye(12), I3/source_denominator
    latent = np.asarray(row["linear_true"])-C.mat(row, "F_LL", 12)@before["linear_true"]
    return a, b, np.r_[latent, gyro_input], np.eye(15), min(denominator, source_denominator)


def metric(row):
    _, p = C.SOURCE.error_and_covariance(row, 21)
    return sym(np.linalg.solve(p, I21)[:18, :18])


def build_word(root, rows, points, mode):
    events = [r for r in rows if r["word"] == mode]
    original = [r for r in points if r["word"] == mode]
    if not original or not events or events[0]["stage"] != "prediction_enter":
        raise ValueError("missing connected word root")
    e0, _ = C.SOURCE.error_and_covariance(original[0], 21)
    running = e0.copy()
    steps, counts, defects = [], Counter(), C.ParityChecks()
    previous = original[0]
    entrance = pending = None
    for row in events:
        stage = row["stage"]
        if stage == "prediction_enter":
            entrance = row
            a, b, u, covariance, den = I21.copy(), np.zeros((21, 0)), np.zeros(0), np.zeros((0, 0)), 1.
        elif stage == "prediction":
            a, b, u, covariance, den = prediction_factor(entrance, row, root["dt"])
        elif stage == "aw_floor":
            a, b, u, covariance, den = I21.copy(), np.zeros((21, 0)), np.zeros(0), np.zeros((0, 0)), 1.
        elif stage == "measurement":
            pending = row
            continue
        elif stage == "projection":
            error, _ = C.SOURCE.error_and_covariance(pending, 21)
            a, b, u, covariance, den = measurement_factor(pending, error)
            if row["kind"] == "S_zero" and pending["R"] != pending["R_S"]:
                raise ValueError("derived gain lost actual R_S")
            counts[row["kind"]] += 1
        else:
            continue  # Injection/reset are included in the completed correction.
        actual, _ = C.SOURCE.error_and_covariance(row, 21)
        before, _ = C.SOURCE.error_and_covariance(previous, 21)
        defects.row = row
        defects("finite_factorization", a@before+b@u, actual, 1+np.linalg.norm(before, np.inf))
        running = a@running+b@u
        defects("composed_word_factorization", running, actual, 1+np.linalg.norm(actual, np.inf))
        steps.append({"A": a, "B": b, "u": u, "Uinv": covariance,
                      "metric": metric(row), "bias_cost": stage == "prediction_enter",
                      "index": row["index"], "stage": stage, "kind": row["kind"],
                      "denominator": float(den)})
        counts[stage] += 1
        previous = row
    if counts["prediction"] != 600 or counts["accelerometer"] != 600:
        raise ValueError("incomplete finite-factor word")
    return steps, metric(original[0]), e0, defects, counts


def supply_update(p, amount):
    """Exact Schur elimination of a bias-energy cost, NOT a bias estimator.

    Equals (P^-1+amount*E_b^T E_b)^-1. Joseph form avoids cancellation;
    P here is the 21-state inverse of minimum input energy, not shipping P.
    """
    noise = I3/amount
    k = np.linalg.solve(p[18:, 18:]+noise, p[18:, :]).T
    f = I21.copy()
    f[:, 18:] -= k
    return sym(f@p@f.T+k@noise@k.T)


def gain_test(steps, m0, dt, factor, gain, *, retain=False):
    """Exact-real elimination of the complete augmented quadratic supply.

    Minimize factor*W0 + gain*(D_b + explicitly normalized forcing energies)
    over all root/port coordinates producing e_j. Its inverse Hessian is P_j.
    The maximum W_j/cost is lambda_max(L_M^T P_j,HH L_M). Thus only forcing
    coordinates are eliminated; all 21 error states/cross terms remain.
    """
    if factor <= 0 or gain <= 0 or dt <= 0 or not steps[0]["bias_cost"]:
        raise ValueError("positive supply parameters and entrance bias cost required")
    p = np.zeros((21, 21))
    p[:18, :18] = np.linalg.solve(factor*m0, np.eye(18))
    p[18:, 18:] = I3/(gain*dt)
    root_p = p.copy()
    ratios, history = [], []
    first = True
    for step in steps:
        before = p.copy() if retain else None
        if step["bias_cost"]:
            if not first:
                p = supply_update(p, gain*dt)
            first = False
        else:
            a, b = step["A"], step["B"]
            p = sym(a@p@a.T+(b@step["Uinv"]@b.T)/gain)
        np.linalg.cholesky(p)
        l = np.linalg.cholesky(step["metric"])
        ratio = float(np.linalg.eigvalsh(sym(l.T@p[:18, :18]@l))[-1])
        if not np.isfinite(ratio) or ratio < 0:
            raise ValueError("nonfinite/negative augmented gain eigenvalue")
        ratios.append(ratio)
        if retain:
            history.append((before, p.copy()))
    return {"endpoint": ratios[-1], "prefix": max(ratios),
            "worst_prefix_step": int(np.argmax(ratios)), "ratios": ratios,
            "P": p, "root_P": root_p, "history": history}


def high_precision_direction(steps, m0, dt, factor, gain, test):
    """80-digit complete-word check of a maximizing COEFFICIENT direction.

    Decimal propagation starts from the actual binary64 coefficient matrices.
    It is not a high-precision nonlinear observer, a source membership check,
    an outward upper bound, or a new direction/metric search.
    """
    l = np.linalg.cholesky(steps[-1]["metric"])
    values, vectors = np.linalg.eigh(sym(l.T@test["P"][:18, :18]@l))
    adjoint = np.r_[l@vectors[:, -1], [0.]*3]/np.sqrt(values[-1])
    inputs = []
    for i in reversed(range(len(steps))):
        step = steps[i]
        if step["bias_cost"]:
            if i:
                state = test["history"][i][1]@adjoint
                adjoint[18:] -= gain*dt*state[18:]
            inputs.append(np.zeros(0))
        else:
            inputs.append(step["Uinv"]@step["B"].T@adjoint/gain)
            adjoint = step["A"].T@adjoint
    x0 = test["root_P"]@adjoint
    inputs.reverse()
    with localcontext() as ctx:
        ctx.prec = 80
        def dec(a):
            return np.asarray([Decimal.from_float(float(x)) for x in np.asarray(a).flat],
                              dtype=object).reshape(np.shape(a))
        x = dec(x0)
        previous_energy = x[:18]@dec(m0)@x[:18]
        cost = Decimal.from_float(factor)*previous_energy
        channels = Counter()
        changes = Counter()
        for step, u in zip(steps, inputs):
            if step["bias_cost"]:
                energy = Decimal.from_float(dt)*(x[18:]@x[18:])
                channels["bias"] += energy
            else:
                du = dec(u)
                if len(u):
                    energy = du@dec(np.linalg.solve(step["Uinv"], np.eye(len(u))))@du
                    channels[step["stage"]+":"+step["kind"]] += energy
                else:
                    energy = Decimal(0)
                x = dec(step["A"])@x+dec(step["B"])@du
            cost += Decimal.from_float(gain)*energy
            current_energy = x[:18]@dec(step["metric"])@x[:18]
            changes[step["stage"]+":"+step["kind"]] += current_energy-previous_energy
            previous_energy = current_energy
        final = x[:18]@dec(steps[-1]["metric"])@x[:18]
        return {"root_direction_full_21": x0.tolist(), "supply_cost_80_digit": str(cost),
                "endpoint_motion_energy_80_digit": str(final),
                "ratio_80_digit": str(final/cost),
                "channel_energies_80_digit": {k: str(v) for k, v in channels.items()},
                "signed_motion_energy_changes_80_digit": {k: str(v) for k, v in changes.items()},
                "SEA3_nonlinear_admissibility_of_maximizer_established": False}


def experiment(steps, m0, dt):
    # A fixed dyadic candidate search on the MATRIX inequality, never on the
    # observed trajectory ratio. Unit normalizations and storage are unchanged.
    candidates = []
    for power in range(0, 41, 4):
        gain = float(2**power)
        test = gain_test(steps, m0, dt, 1., gain)
        candidates.append({"gain": gain, "endpoint_at_factor_one": test["endpoint"]})
        if test["endpoint"] < 1-1e-8:
            break
    factor = (1+test["endpoint"])/2 if test["endpoint"] < 1-1e-8 else 1.
    endpoint = gain_test(steps, m0, dt, factor, gain, retain=True)
    witness = high_precision_direction(steps, m0, dt, factor, gain, endpoint)
    if abs(float(witness["ratio_80_digit"])-endpoint["endpoint"]) > 1e-7*max(1., endpoint["endpoint"]):
        raise ValueError("80-digit maximizing-direction check disagrees with augmented eigenvalue")
    prefix_gain = gain
    for power in range(int(np.log2(gain)), 41, 4):
        prefix_gain = float(2**power)
        prefix = gain_test(steps, m0, dt, 2., prefix_gain)
        if prefix["prefix"] < 1-1e-8:
            break
    accepted = factor < 1 and endpoint["endpoint"] < 1 and prefix["prefix"] < 1
    worst = steps[prefix["worst_prefix_step"]]
    source_energy = sum(float(s["u"]@np.linalg.solve(s["Uinv"], s["u"]))
                        for s in steps if len(s["u"]))
    duration = dt*sum(s["bias_cost"] for s in steps)
    # Lower bound on the bound obtained by the triangle-bound/geometric-series
    # composition: B_true>=0 and every uniform source budget covers this point.
    # This is NOT a lower bound on the actual filter's error or a uniform gain.
    floor_bound = gain*(duration*.4**2+source_energy)/(1-factor) if factor < 1 else None
    return {"decision": "FROZEN_COEFFICIENT_SUPPLY_FEASIBLE" if accepted else "FROZEN_COEFFICIENT_SUPPLY_NOT_CLOSED",
            "endpoint_factor": factor, "endpoint_channel_gains": dict.fromkeys(SUPPLY_UNITS, gain),
            "endpoint_augmented_ratio": endpoint["endpoint"], "endpoint_margin": 1-endpoint["endpoint"],
            "prefix_factor": 2., "prefix_channel_gains": dict.fromkeys(SUPPLY_UNITS, prefix_gain),
            "every_prefix_augmented_ratio_max": prefix["prefix"],
            "limiting_prefix": {k: worst[k] for k in ("index", "stage", "kind")},
            "candidate_matrix_tests": candidates, "maximizing_coefficient_direction": witness,
            "point_normalized_source_energy": source_energy,
            "candidate_composed_bound_at_least_if_gains_extended_uniformly": floor_bound,
            "composed_bound_lower_bound_is_not_observed_error_floor": True,
            "nonlinear_graph_or_source_family_certified": False,
            "useful_uniform_accuracy_floor_certified": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--attachment", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    trace_path = Path(str(args.prefix)+".events.jsonl")
    attachment = json.loads(args.attachment.read_text())
    digest = hashlib.sha256(trace_path.read_bytes()).hexdigest()
    if attachment["event_trace_sha256"] != digest or not attachment["read_only_trace_recovers_baseline_bit_for_bit"]:
        raise ValueError("gain experiment detached from connected attachment evidence")
    for suffix in (".root.json", ".inputs.csv", ".prefixes.jsonl"):
        if hashlib.sha256(Path(str(args.prefix)+suffix).read_bytes()).hexdigest() != attachment["baseline_capture_sha256"][suffix]:
            raise ValueError("source/root/prefix input detached: "+suffix)
    rows = [json.loads(line) for line in trace_path.read_text().splitlines()]
    points = [json.loads(line) for line in Path(str(args.prefix)+".prefixes.jsonl").read_text().splitlines()]
    root = json.loads(Path(str(args.prefix)+".root.json").read_text())
    report = {"experiment": "CONNECTED_FINITE_COEFFICIENT_MOTION_SUPPLY", "event_trace_sha256": digest,
              "supply_units": SUPPLY_UNITS, "coefficient_dependencies_frozen_only_for_point_test": True,
              "forcing_ports_relaxed_only_for_sufficient_coefficient_test": True,
              "independent_SEA3_sample_boxes_created": False, "bias_error_decay_required": False,
              "P4_MOTION_PASS": False, "P5_MOTION_MAY_START": False, "modes": {}}
    for mode in ("H18", "A21"):
        if attachment["modes"][mode]["decision"] != "CONNECTED_POINT_ATTACHMENT_PASS":
            raise ValueError("connected subevents have not passed: "+mode)
        steps, m0, _, defects, counts = build_word(root, rows, points, mode)
        item = {"counts": dict(counts), "factorization_defects": dict(defects.defects),
                "factorization_failure": defects.failures[:1],
                "minimum_chart_denominator": min(s["denominator"] for s in steps)}
        report["modes"][mode] = item
        if not defects.failures:
            item.update(experiment(steps, m0, root["dt"]))
        else:
            item["decision"] = "FINITE_FACTORIZATION_ATTACHMENT_FAIL"
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        print("CONNECTED_MOTION_GAIN", mode, json.dumps(item, allow_nan=False), flush=True)
        if defects.failures:
            raise SystemExit("finite factorization failed; augmented experiment not authorized")


if __name__ == "__main__":
    main()
