#!/usr/bin/env python3
"""Extract the OU-III continuous-source proof domain from implementation guards.

This producer intentionally does not infer bounds from the eight reference
trajectories. It parses the shipping implementation's safety/clamp constants
and records every discrete branch and hybrid transition the validated backend
must cover.

The implementation constants are C++ ``float`` constexpr expressions. For a
proof-domain boundary it is not sufficient to evaluate those expressions as
Python binary64. This parser therefore evaluates every literal and arithmetic
operation as IEEE-754 binary32, using exact rationals between operations and an
explicit nearest/ties-to-even rounding step. The returned Python float is an
exact binary64 representation of the deployed binary32 value.

Validated continuous-word building blocks included here are:
- exact-rational Taylor/Lagrange bounds for the OU scalar primitives;
- the exact 4x4 [v,p,S,a] IntegratedOUChain transition enclosure;
- a positive-kernel exact-rational enclosure of the mathematical OU Qd;
- exact-rational interval covariance prediction Phi P Phi^T + Qd; and
- the scalar-axis mathematical S=0 pseudo-measurement Riccati update.

The scalar pseudo update is a building block for the shipping 3-D update. The
actual filter forms the 3x3 innovation P_SS+R_S, solves it with LDLT, and applies
a Joseph update, so cross-axis coupling and shipping floating-point factorization
remain explicit obligations before theorem promotion.

The public wrapper accepts every positive finite binary32 ``dt``. That set is
technically finite -- [2^-149, FLT_MAX] -- but it has no operational safety
upper guard. A deployment theorem therefore still needs a source/configuration
supported-step bound that is small enough for finite, physically meaningful
transition and covariance enclosures; nominal 200 Hz is not silently assumed.
"""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction
import json
import math
import re
import struct
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_HEADER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"

REQUIRED = (
    "MIN_TUNE_FREQ_HZ", "MAX_TUNE_FREQ_HZ", "MIN_TAU_S", "MAX_TAU_S",
    "MAX_SIGMA_A", "MIN_R_S", "MAX_R_S",
    "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT", "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
    "MAG_DELAY_SEC", "ONLINE_TUNE_WARMUP_SEC",
)
HYBRID_OBLIGATIONS = (
    "startup_handoff", "held_to_active", "magnetic_lock",
    "magnetic_regauge_refinement", "tilt_reset", "tilt_relock",
    "cooldown_reentry", "periodic_aw_covariance_sync",
)
CONST_RE = re.compile(r"constexpr\s+float\s+([A-Za-z_]\w*)\s*=\s*([^;]+);", re.MULTILINE)
_FLOAT32_MAX_BITS = 0x7F7FFFFF


def _strip_float_suffixes(expr: str) -> str:
    return re.sub(r"(?<=\d)[fF]\b", "", expr)


def _bits_to_positive_fraction(bits: int) -> Fraction:
    if not 0 <= bits <= _FLOAT32_MAX_BITS:
        raise ValueError(f"not a finite positive binary32 pattern: 0x{bits:08x}")
    exponent = (bits >> 23) & 0xFF
    mantissa = bits & 0x7FFFFF
    if exponent == 0:
        return Fraction(mantissa, 1 << 149)
    significand = (1 << 23) | mantissa
    power = exponent - 127 - 23
    return Fraction(significand << power, 1) if power >= 0 else Fraction(significand, 1 << (-power))


def _round_fraction_binary32(value: Fraction) -> Fraction:
    if value == 0:
        return Fraction(0, 1)
    sign = -1 if value < 0 else 1
    x = abs(value)
    try:
        candidate = struct.unpack(">I", struct.pack(">f", float(x)))[0]
    except OverflowError as exc:
        raise RuntimeError(f"binary32 constant overflow for {value}") from exc
    candidate &= 0x7FFFFFFF
    if candidate > _FLOAT32_MAX_BITS:
        raise RuntimeError(f"binary32 constant overflow for {value}")
    choices = [(_bits_to_positive_fraction(bits), bits)
               for bits in (candidate - 1, candidate, candidate + 1)
               if 0 <= bits <= _FLOAT32_MAX_BITS]
    if not choices:
        raise RuntimeError(f"cannot round binary32 constant {value}")
    exact, _ = min(choices, key=lambda item: (abs(item[0] - x), item[1] & 1))
    return exact if sign > 0 else -exact


def _literal_fraction(expr: str, node: ast.Constant) -> Fraction:
    token = ast.get_source_segment(expr, node) or repr(node.value)
    try:
        return Fraction(token)
    except (ValueError, ZeroDivisionError) as exc:
        raise RuntimeError(f"cannot parse constexpr literal {token!r}") from exc


def _eval_constexpr32(node: ast.AST, expr: str, text: str, stack: tuple[str, ...]) -> Fraction:
    if isinstance(node, ast.Expression):
        return _eval_constexpr32(node.body, expr, text, stack)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return _round_fraction_binary32(_literal_fraction(expr, node))
    if isinstance(node, ast.Name):
        return parse_const_fraction(text, node.id, stack)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        x = _eval_constexpr32(node.operand, expr, text, stack)
        return x if isinstance(node.op, ast.UAdd) else -x
    if isinstance(node, ast.BinOp):
        a = _eval_constexpr32(node.left, expr, text, stack)
        b = _eval_constexpr32(node.right, expr, text, stack)
        if isinstance(node.op, ast.Add): exact = a + b
        elif isinstance(node.op, ast.Sub): exact = a - b
        elif isinstance(node.op, ast.Mult): exact = a * b
        elif isinstance(node.op, ast.Div):
            if b == 0: raise ZeroDivisionError("constexpr division by zero")
            exact = a / b
        else: raise RuntimeError(f"unsupported constexpr operator: {ast.dump(node.op)}")
        return _round_fraction_binary32(exact)
    raise RuntimeError(f"unsupported constexpr expression node: {ast.dump(node)}")


def parse_const_fraction(text: str, name: str, stack: tuple[str, ...] = ()) -> Fraction:
    if name in stack:
        raise RuntimeError(f"cyclic implementation constant alias: {' -> '.join((*stack, name))}")
    expressions = {n: expr for n, expr in CONST_RE.findall(text)}
    if name not in expressions:
        raise RuntimeError(f"cannot extract implementation constant {name}")
    expr = " ".join(_strip_float_suffixes(expressions[name]).split())
    try:
        return _eval_constexpr32(ast.parse(expr, mode="eval"), expr, text, (*stack, name))
    except (SyntaxError, ZeroDivisionError) as exc:
        raise RuntimeError(f"cannot evaluate implementation constant {name}: {expr!r}") from exc


def parse_const(text: str, name: str, stack: tuple[str, ...] = ()) -> float:
    value = float(parse_const_fraction(text, name, stack))
    if not math.isfinite(value): raise RuntimeError(f"implementation constant {name} is non-finite")
    return value


def parse_aw_sigma_floor(text: str) -> float:
    pat = re.compile(r"const\s+float\s+sigma_floor\s*=\s*std::max\(\s*([0-9.+\-eE]+)f?\s*,\s*band_noise_floor_sigma_\(\)\s*\)\s*;")
    m = pat.search(text)
    if not m: raise RuntimeError("cannot extract deployed a_w stationary-std floor")
    return float(_round_fraction_binary32(Fraction(m.group(1))))


def _outward_point(x: float) -> list[float]:
    x = float(x)
    if not math.isfinite(x): raise RuntimeError(f"non-finite source-domain endpoint {x!r}")
    return [math.nextafter(x, -math.inf), math.nextafter(x, math.inf)]


def _outward_box(lo: float, hi: float) -> list[float]:
    lo, hi = float(lo), float(hi)
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
        raise RuntimeError(f"invalid source interval [{lo!r}, {hi!r}]")
    return [math.nextafter(lo, -math.inf), math.nextafter(hi, math.inf)]


def _fraction_down(q: Fraction) -> float:
    y = float(q)
    return math.nextafter(y, -math.inf) if Fraction.from_float(y) > q else y


def _fraction_up(q: Fraction) -> float:
    y = float(q)
    return math.nextafter(y, math.inf) if Fraction.from_float(y) < q else y


def _exp_neg_point_rational(x: Fraction, order: int = 96) -> tuple[Fraction, Fraction]:
    if x < 0: raise ValueError("exp(-x) enclosure requires x>=0")
    if order < 1: raise ValueError("Taylor order must be positive")
    term = total = Fraction(1, 1)
    for n in range(1, order + 1):
        term *= -x
        term /= n
        total += term
    rem = x ** (order + 1) / math.factorial(order + 1)
    return total - rem, total + rem


def validated_ou_primitives(h_bounds, tau_bounds, order: int = 96) -> dict:
    h_lo, h_hi = map(Fraction.from_float, map(float, h_bounds))
    t_lo, t_hi = map(Fraction.from_float, map(float, tau_bounds))
    if h_lo < 0 or h_hi < h_lo: raise ValueError("invalid nonnegative step interval")
    if t_lo <= 0 or t_hi < t_lo: raise ValueError("invalid positive tau interval")
    x_lo, x_hi = h_lo / t_hi, h_hi / t_lo
    alpha_lo, _ = _exp_neg_point_rational(x_hi, order)
    _, alpha_hi = _exp_neg_point_rational(x_lo, order)
    alpha_lo, alpha_hi = max(Fraction(0), alpha_lo), min(Fraction(1), alpha_hi)
    pa_lo = max(Fraction(0), x_lo + alpha_lo - 1)
    pa_hi = max(Fraction(0), x_hi + alpha_hi - 1)
    sa_lo = max(Fraction(0), Fraction(1, 2)*x_lo*x_lo - x_hi - alpha_hi + 1)
    sa_hi = max(Fraction(0), Fraction(1, 2)*x_hi*x_hi - x_lo - alpha_lo + 1)
    t2lo, t2hi = t_lo*t_lo, t_hi*t_hi
    t3lo, t3hi = t2lo*t_lo, t2hi*t_hi
    return {
        "validated_arithmetic": True, "outward_rounded": True,
        "backend": f"EXACT_RATIONAL_TAYLOR_ORDER_{order}_LAGRANGE_REMAINDER",
        "h_s": [_fraction_down(h_lo), _fraction_up(h_hi)],
        "tau_s": [_fraction_down(t_lo), _fraction_up(t_hi)],
        "x_h_over_tau": [_fraction_down(x_lo), _fraction_up(x_hi)],
        "alpha": [_fraction_down(alpha_lo), _fraction_up(alpha_hi)],
        "phi_pa_s2": [_fraction_down(t2lo*pa_lo), _fraction_up(t2hi*pa_hi)],
        "phi_Sa_s3": [_fraction_down(t3lo*sa_lo), _fraction_up(t3hi*sa_hi)],
    }


def _mul_nonnegative_bounds(a, b):
    if a[0] < 0 or b[0] < 0: raise ValueError("negative lower bound")
    return [math.nextafter(float(a[0])*float(b[0]), -math.inf),
            math.nextafter(float(a[1])*float(b[1]), math.inf)]


def validated_phi_axis4(h_bounds, tau_bounds, order: int = 96) -> dict:
    p = validated_ou_primitives(h_bounds, tau_bounds, order)
    h_lo, h_hi = map(float, h_bounds); t_lo, t_hi = map(float, tau_bounds)
    zero, one = [0.0, 0.0], [1.0, 1.0]
    h = _outward_box(h_lo, h_hi)
    half_h2 = [math.nextafter(0.5*h_lo*h_lo, -math.inf), math.nextafter(0.5*h_hi*h_hi, math.inf)]
    one_minus_alpha = [max(0.0, math.nextafter(1-p["alpha"][1], -math.inf)), math.nextafter(1-p["alpha"][0], math.inf)]
    phi_va = _mul_nonnegative_bounds(_outward_box(t_lo, t_hi), one_minus_alpha)
    M = [[one, zero, zero, phi_va], [h, one, zero, p["phi_pa_s2"]],
         [half_h2, h, one, p["phi_Sa_s3"]], [zero, zero, zero, p["alpha"]]]
    return {"validated_arithmetic": True, "outward_rounded": True,
            "implementation": "IntegratedOUChain<T,3>::transition",
            "state_order": ["v", "p", "S", "a_w"], "scalar_primitives": p, "Phi_interval": M}


def _kernel_interval_on_cell(r_lo, r_hi, tau_bounds, order):
    M = validated_phi_axis4((_fraction_down(r_lo), _fraction_up(r_hi)), tau_bounds, order)["Phi_interval"]
    ans=[]
    for i in range(4):
        blo=max(0.0,float(M[i][3][0])); bhi=max(blo,float(M[i][3][1]))
        ans.append((Fraction.from_float(blo), Fraction.from_float(bhi)))
    return ans


def _integrate_kernel_bounds(h, tau_bounds, cells, order):
    if h < 0 or cells <= 0: raise ValueError("invalid integration domain")
    lower=[[Fraction(0) for _ in range(4)] for _ in range(4)]
    upper=[[Fraction(0) for _ in range(4)] for _ in range(4)]
    if h == 0: return lower, upper
    width=h/cells
    for k in range(cells):
        g=_kernel_interval_on_cell(width*k, width*(k+1), tau_bounds, order)
        for i in range(4):
            for j in range(i,4):
                lo=width*g[i][0]*g[j][0]; hi=width*g[i][1]*g[j][1]
                lower[i][j]+=lo; upper[i][j]+=hi
                if i!=j: lower[j][i]+=lo; upper[j][i]+=hi
    return lower, upper


def validated_qd_axis4_kernel(h_bounds, tau_bounds, sigma2_bounds, cells: int = 24, order: int = 96) -> dict:
    h_lo,h_hi=map(Fraction.from_float,map(float,h_bounds)); t_lo,t_hi=map(Fraction.from_float,map(float,tau_bounds)); s_lo,s_hi=map(Fraction.from_float,map(float,sigma2_bounds))
    if h_lo<0 or h_hi<h_lo or t_lo<=0 or t_hi<t_lo or s_lo<0 or s_hi<s_lo: raise ValueError("invalid Qd domain")
    int_lo,_=_integrate_kernel_bounds(h_lo,tau_bounds,cells,order); _,int_hi=_integrate_kernel_bounds(h_hi,tau_bounds,cells,order)
    qc_lo=Fraction(2)*s_lo/t_hi; qc_hi=Fraction(2)*s_hi/t_lo
    Q=[]
    for i in range(4):
        row=[]
        for j in range(4): row.append([_fraction_down(qc_lo*int_lo[i][j]), _fraction_up(qc_hi*int_hi[i][j])])
        Q.append(row)
    return {"validated_arithmetic":True,"outward_rounded":True,"claim":"MATHEMATICAL_INTEGRATED_OU_PROCESS_COVARIANCE_ENCLOSURE",
            "backend":"POSITIVE_KERNEL_CELL_INTERVAL_EXACT_RATIONAL_ACCUMULATION","cells":cells,"state_order":["v","p","S","a_w"],
            "h_s":[_fraction_down(h_lo),_fraction_up(h_hi)],"tau_s":[_fraction_down(t_lo),_fraction_up(t_hi)],"sigma2_aw":[_fraction_down(s_lo),_fraction_up(s_hi)],
            "Qd_interval":Q,"mathematical_integral_enclosed":True,"shipping_binary32_closed_form_enclosed":False,"shipping_psd_cleanup_enclosed":False,"theorem_promotion":"NOT_ESTABLISHED"}


def _iv_from_float(b):
    lo,hi=map(float,b)
    if not(math.isfinite(lo) and math.isfinite(hi)) or lo>hi: raise ValueError("invalid interval")
    return Fraction.from_float(lo),Fraction.from_float(hi)
def _iv_add(a,b): return a[0]+b[0],a[1]+b[1]
def _iv_sub(a,b): return a[0]-b[1],a[1]-b[0]
def _iv_mul(a,b):
    p=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]); return min(p),max(p)
def _iv_div_positive(a,b):
    if b[0] <= 0: raise ValueError("denominator interval must be positive")
    return _iv_mul(a,(Fraction(1,1)/b[1],Fraction(1,1)/b[0]))
def _iv_matrix_from_float(M): return [[_iv_from_float(x) for x in row] for row in M]
def _iv_transpose(M): return [list(row) for row in zip(*M)]
def _iv_matmul(A,B):
    if not A or not B or len(A[0])!=len(B): raise ValueError("matrix dimensions")
    out=[]
    for i in range(len(A)):
        row=[]
        for j in range(len(B[0])):
            acc=(Fraction(0),Fraction(0))
            for k in range(len(B)): acc=_iv_add(acc,_iv_mul(A[i][k],B[k][j]))
            row.append(acc)
        out.append(row)
    return out


def validated_covariance_predict_axis4(P_interval,h_bounds,tau_bounds,sigma2_bounds,cells: int=24,order: int=96)->dict:
    if len(P_interval)!=4 or any(len(r)!=4 for r in P_interval): raise ValueError("P_interval must be 4x4")
    phi=validated_phi_axis4(h_bounds,tau_bounds,order); qd=validated_qd_axis4_kernel(h_bounds,tau_bounds,sigma2_bounds,cells,order)
    A=_iv_matrix_from_float(phi["Phi_interval"]); P=_iv_matrix_from_float(P_interval); Q=_iv_matrix_from_float(qd["Qd_interval"])
    APA=_iv_matmul(_iv_matmul(A,P),_iv_transpose(A)); pred=[]
    for i in range(4):
        row=[]
        for j in range(4):
            lo,hi=_iv_add(APA[i][j],Q[i][j]); row.append([_fraction_down(lo),_fraction_up(hi)])
        pred.append(row)
    return {"validated_arithmetic":True,"outward_rounded":True,"claim":"MATHEMATICAL_OU_AXIS_COVARIANCE_PREDICTION_ENCLOSURE",
            "state_order":["v","p","S","a_w"],"Phi":phi,"Qd":qd,"P_prior_interval":P_interval,"P_predicted_interval":pred,
            "mathematical_prediction_enclosed":True,"shipping_covariance_prediction_enclosed":False,"theorem_promotion":"NOT_ESTABLISHED"}


def validated_scalar_s_pseudo_update_axis4(P_interval,R_bounds)->dict:
    """Enclose exact scalar S=0 covariance update for one [v,p,S,a] axis."""
    if len(P_interval)!=4 or any(len(r)!=4 for r in P_interval): raise ValueError("P_interval must be 4x4")
    P=_iv_matrix_from_float(P_interval); R=_iv_from_float(list(R_bounds))
    if R[0] <= 0: raise ValueError("R_S must be strictly positive")
    innov=_iv_add(P[2][2],R)
    if innov[0] <= 0: raise ValueError("P_SS+R_S lacks a positive lower bound")
    post=[]
    for i in range(4):
        row=[]
        for j in range(4):
            corr=_iv_div_positive(_iv_mul(P[i][2],P[2][j]),innov)
            lo,hi=_iv_sub(P[i][j],corr)
            row.append([_fraction_down(lo),_fraction_up(hi)])
        post.append(row)
    return {"validated_arithmetic":True,"outward_rounded":True,
            "claim":"MATHEMATICAL_SCALAR_AXIS_S_ZERO_RICCATI_UPDATE_ENCLOSURE",
            "state_order":["v","p","S","a_w"],"measurement":"S=0","R_S_interval":list(R_bounds),
            "innovation_variance_interval":[_fraction_down(innov[0]),_fraction_up(innov[1])],
            "P_prior_interval":P_interval,"P_posterior_interval":post,
            "mathematical_scalar_update_enclosed":True,"shipping_3d_ldlt_joseph_update_enclosed":False,"theorem_promotion":"NOT_ESTABLISHED"}


def build(header: Path) -> dict:
    text=header.read_text(); c={name:parse_const(text,name) for name in REQUIRED}; sigma_floor=parse_aw_sigma_floor(text)
    continuous={"wave_tune_frequency_hz":[c["MIN_TUNE_FREQ_HZ"],c["MAX_TUNE_FREQ_HZ"]],"tau_aw_s":[c["MIN_TAU_S"],c["MAX_TAU_S"]],
                "sigma_aw_mps2":[sigma_floor,c["MAX_SIGMA_A"]],"R_S_base":[c["MIN_R_S"],c["MAX_R_S"]],
                "pseudo_update_period_s":[c["PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"],c["PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"]]}
    timing={"mag_delay":c["MAG_DELAY_SEC"],"online_tune_warmup":c["ONLINE_TUNE_WARMUP_SEC"]}
    parameter_box={"qualification":"SOURCE_DERIVED_OUTWARD_ROUNDED_PARAMETER_BOX","validated_arithmetic":True,"outward_rounded":True,
                   "arithmetic_backend":"EXACT_BINARY32_SOURCE_PLUS_BINARY64_NEXTAFTER_OUTWARD",
                   "continuous_parameters":{n:_outward_box(b[0],b[1]) for n,b in continuous.items()},
                   "timing_constants_s":{n:_outward_point(v) for n,v in timing.items()},"continuous_word_enclosed":False,"nonlinear_word_enclosed":False,"theorem_promotion":"NOT_ESTABLISHED"}
    min_positive_f32=float(_bits_to_positive_fraction(1)); max_finite_f32=float(_bits_to_positive_fraction(_FLOAT32_MAX_BITS))
    return {"schema":2,"claim":"OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT","source_generated_not_trajectory_fit":True,"source_complete_parameter_domain":True,
            "validated_arithmetic":False,"outward_rounded":False,"implementation_header":str(header.relative_to(REPO)),
            "implementation_scalar_semantics":{"type":"IEEE754_BINARY32","rounding":"ROUND_TO_NEAREST_TIES_TO_EVEN_EACH_OPERATION","evaluation":"EXACT_RATIONAL_THEN_BINARY32_ROUND"},
            "continuous_parameters":continuous,"timing_constants_s":timing,"validated_parameter_box":parameter_box,
            "accepted_update_step_domain_s":{"lower_closed":min_positive_f32,"upper_closed":max_finite_f32,"type_level_finite_upper_bound":True,
                "operational_safety_upper_guard":False,"proof_usable_supported_upper_bound":False,
                "implementation_observation":"SeaStateFusionFilter_OU_III::updateCore_ accepts every positive finite binary32 dt; the IEEE-754 type supplies a finite maximum but the implementation supplies no operationally meaningful maximum supported step",
                "theorem_effect":"a source-complete type-level domain exists, but its FLT_MAX upper endpoint is not a finite-state operational domain for the transition/covariance formulas; theorem promotion requires a source/configuration supported-step guard or equivalent caller contract"},
            "validated_ou_primitive_backend":{"available":True,"requires_proof_usable_h_upper":True,"backend":"EXACT_RATIONAL_TAYLOR_WITH_LAGRANGE_REMAINDER",
                "transition_matrix_backend":"VALIDATED_INTERVAL_4X4_INTEGRATED_OU_CHAIN","process_covariance_backend":"POSITIVE_KERNEL_CELL_INTERVAL_EXACT_RATIONAL_ACCUMULATION",
                "process_covariance_mathematical_integral_enclosed":True,"process_covariance_shipping_float_path_enclosed":False,
                "covariance_prediction_backend":"EXACT_RATIONAL_INTERVAL_PHI_P_PHIT_PLUS_Q","covariance_prediction_mathematical_path_enclosed":True,
                "scalar_S_zero_update_backend":"EXACT_RATIONAL_SCALAR_RICCATI_UPDATE","scalar_S_zero_mathematical_path_enclosed":True,
                "shipping_3d_S_zero_update_enclosed":False,"theorem_promotion":"BLOCKED_BY_NO_PROOF_USABLE_ACCEPTED_DT_GUARD_AND_SHIPPING_QD_FLOAT_PATH"},
            "discrete_source_branches":{"mode":["H","A"],"accelerometer_gate":["accepted","rejected"],"magnetometer_gate":["not_due","accepted","rejected"],
                "S_zero_pseudo":["not_due","due"],"magnetic_gauge":["unlocked","locked","refined"],"tilt_recovery":["normal","reset","relock","cooldown_reentry"],
                "aw_covariance_sync":["not_due","due_psd_increment"]},"hybrid_obligations":list(HYBRID_OBLIGATIONS),
            "periodic_aw_covariance_sync_proof":{"required_mode":"PSD_NONEXPANSIVE","operation":"P_plus=P_minus+E_a Delta_plus E_a^T with Delta_plus>=0",
                "metric_consequence":"inverse-covariance information energy is nonexpansive"},
            "promotion_rule":"source-boundary scalar arithmetic, exact per-axis transition, mathematical OU process covariance, one-step covariance prediction and scalar-axis S=0 Riccati update are enclosed; deployment promotion still requires a proof-usable accepted-dt guard/contract, shipping binary32 Qd/PSD-cleanup and 3-D LDLT/Joseph enclosure, full H/A words, nonlinear SO(3) bounds and remaining jump certificates"}


def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--header",type=Path,default=DEFAULT_HEADER); ap.add_argument("--output",type=Path,required=True); args=ap.parse_args()
    payload=build(args.header.resolve()); args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)); print(json.dumps(payload,indent=2,sort_keys=True)); return 0


if __name__=="__main__": raise SystemExit(main())
