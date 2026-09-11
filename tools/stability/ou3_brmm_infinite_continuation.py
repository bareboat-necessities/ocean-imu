#!/usr/bin/env python3
"""Exact obstruction to bounded all-18 errors under the finite-window source.

This is a necessary-condition certificate, NOT an endpoint/P4 certificate.
The OLD finite-window-only source bounded p and every short-window Delta S, but not S_L for
an entire continuation.  The admitted quiet subfamily p=d, v=a=0 has S_L=h*d.
All d in one closed position cell drive the SAME shipping sensor execution.
Consequently e_p=d-p_hat and e_S,L=h*d-S_hat at every completed event.

The exact polynomial graph below retains the shared physical/error column;
it never replaces a reachable P/H/R/K or tuner tuple by independent boxes.
The information argument is valid for any common execution, including the
binary32 implementation: physical p and S are not executable inputs.  It does
not prove source admission of every Normal-Live branch or exhibit an unstable
nominal trajectory.  It refutes the *additional* finite ultimate-bound and
finite working-tube conclusion for this explicitly defined source family.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import ou3_brmm_physical_wave_source as WAVE
import ou3_p4_bias0_family as B0
import ou3_p4_bias1_family as B1
import ou3_p4_bias2_family as B2

ROOT = Path(__file__).resolve().parents[2]
QUALIFICATION = 'OU3_BRMM_INDEFINITE_S_NECESSARY_CONDITION_V1'
# Pins bind unchanged shipping sensor/entry semantics. They are
# NOT a substitute for the proof in w3d-brmm-stability-theorem.tex-part.
AUDITED = {
    'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h': '1008e931734f226f93f52e46a1408503a7c1be9b364c0756da64a19788ece5ed',
    'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h': 'c9ed955ef998992e33b253ed0d5d49673852b87413a8c560479a9b482786d64b',
    'src/kalman_common/SeaStateFusionFilterCommon.h': 'f76b6266ab4f403d2cce61058a79f1fdb5bb55aab1355c62ce2adaf516d7ea9f',
}

# Sparse exact polynomials over Q, with sorted tuples of variable names as
# monomials. No sampling, interval widening, or floating-point positivity test.
def constant(x):
    x = F(x)
    return {(): x} if x else {}


def symbol(name):
    return {(name,): F(1)}


def add(a, b, factor=F(1)):
    out = dict(a)
    for monomial, value in b.items():
        out[monomial] = out.get(monomial, F(0)) + factor*value
    return {k: v for k, v in out.items() if v}


def mul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(sorted(ka+kb))
            out[key] = out.get(key, F(0)) + va*vb
    return {k: v for k, v in out.items() if v}


def encoded(p):
    return {'*'.join(k) or '1': str(v) for k, v in sorted(p.items())}


def ambiguity_map(n, h):
    """Exact coefficients of the single shared d in e_p and e_S,L."""
    if n not in (18, 21, 24):
        raise ValueError('expected motion18, error21, or joint error/true-bias24')
    return [[constant(int(i == 9+j)) if i < 12 else
             (h if i == 12+j else {}) for j in range(3)] for i in range(n)]


def event_certificate(n, *, physical_S_sign=1):
    """Coefficientwise event identities, for all h>=0 and all closed radial d.

    Formal gain symbols are used only to establish an identity valid for ANY
    gain; they are not free gains for a stability theorem. Production gains
    remain one common output of the identical P/H/R execution.
    """
    h, dt = symbol('h'), symbol('dt_physical')
    B = ambiguity_map(n, h)
    # Prediction moves the PHYSICAL jet by its actual duration; nominal
    # prediction and its rounding cancel between identical executions.
    predicted = [[add(B[i][j], dt if i == 12+j else {})
                  for j in range(3)] for i in range(n)]
    expected = ambiguity_map(n, add(h, dt))
    identities = [add(predicted[i][j], expected[i][j], -1)
                  for i in range(n) for j in range(3)]
    residuals = {}
    for event in ('S_zero', 'accelerometer', 'magnetometer'):
        if event == 'S_zero':
            residual = [[add(B[12+i][j], h if i == j else {}, -physical_S_sign)
                         for j in range(3)] for i in range(3)]
        else:
            # Neither vector model observes p or S. All other H entries may
            # be arbitrary, but multiply zero in this shared source column.
            residual = [[{} for _ in range(3)] for _ in range(3)]
            for r in range(3):
                for j in range(3):
                    for k in range(n):
                        H = {} if 9 <= k < 15 else symbol(f'H_{event}_{r}_{k}')
                        residual[r][j] = add(residual[r][j], mul(H, B[k][j]))
        residuals[event] = [[encoded(x) for x in row] for row in residual]
        for i in range(n):
            for j in range(3):
                corrected = B[i][j]
                for r in range(3):
                    corrected = add(corrected, mul(symbol(f'K_{event}_{i}_{r}'), residual[r][j]), -1)
                identities.append(add(corrected, B[i][j], -1))
    # A finite common nominal output, not a differential/reset approximation:
    # (truth_base + B*d - y_common) - (truth_base - y_common) = B*d.
    # This covers floor/reset/projection/release/reinitialization/identity
    # branches even when they rewrite nominal states or covariance.
    for i in range(n):
        y = symbol(f'common_finite_output_{i}')
        for j in range(3):
            image = add(add(B[i][j], y, -1), y)
            identities.append(add(image, B[i][j], -1))
    return {
        'dimension': n,
        'coefficient_map_B_h': [[encoded(x) for x in row] for row in B],
        'prediction': 'B(h+dt_physical), independent of nominal roundoff',
        'Joseph_residual_root_coefficients': residuals,
        'finite_common_output_events': [
            'aw_floor', 'finite_attitude_reset', 'bias_projection', 'H18_to_A21',
            'late_north_or_tilt_reinitialization', 'not_due', 'rejected'],
        'polynomial_identities_checked': len(identities),
        'nonzero_polynomial_residuals': [encoded(x) for x in identities if x],
    }


def source_energy(h, word, displacement):
    """Integral |S_L|^2 on a word; not a substitute for the actual ISS port cost."""
    h, word, d = map(F, (h, word, displacement))
    if h < 0 or word <= 0:
        raise ValueError('nonnegative elapsed time and positive word required')
    return d*d*((h+word)**3-h**3)/3


def paired_errors(h, displacement, common_S_estimate):
    h, d, estimate = map(F, (h, displacement, common_S_estimate))
    if h < 0:
        raise ValueError('negative elapsed time')
    return h*d-estimate, -h*d-estimate


def separate_zero_bias_witnesses():
    records = {}
    # Invoke/validate each family independently. No BIAS0/2 inference from BIAS1.
    for name, module in (('BIAS0', B0), ('BIAS1', B1), ('BIAS2', B2)):
        d = module.build()
        failures = module.validate(d)
        if failures or not d[name+'_SOURCE_ADMISSION_PASS']:
            raise ValueError(f'{name} prerequisite failed: {failures}')
        records[name] = {
            'qualification': d['qualification'],
            'own_family_admission': d[name+'_SOURCE_ADMISSION_PASS'],
            'true_bias': '0', 'driver': '0',
            'recurrence_residual': encoded(add(mul(symbol(name+'_phi'), constant(0)), constant(0))),
            'construction': {
                'BIAS0': 'GM, turn-on, thermal, strain and non-GM components each zero',
                'BIAS1': 'exponential root, sinusoid amplitude and deterministic component each zero',
                'BIAS2': 'zero bounded-variation drift, including its own phi_true=1 member',
            }[name],
        }
    return records


def build():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in AUDITED}
    changed = [p for p in AUDITED if hashes[p] != AUDITED[p]]
    if changed:
        raise ValueError('re-audit the source/observation semantics after changes: '+repr(changed))
    # Historical definition at main f7123bd8: p'=v, v'=a, S'=p,
    # bounded p/v and short-window increments only. The constant d is
    # independently checked against a conservative subset of its position cap.
    # Do NOT import current BRMM admission here and call this an admitted witness.
    d = F(1, 8)
    pm = F(7)  # strictly below the old sqrt(3)*8.5/2 position cap
    wave = WAVE.build()
    rejection = WAVE.constant_history_admission((d, 0, 0),
        wave['analytic_examples_not_full_source_qualification']['continuum_band'])
    if not 0 < d < pm:
        raise ValueError('quiet witness must be inside the unchanged position domain')
    events = [event_certificate(n) for n in (18, 21, 24)]
    exact = all(not x['nonzero_polynomial_residuals'] for x in events)
    # Parallelogram identity in the S coordinate, with arbitrary common output.
    h, sh = symbol('h'), symbol('S_hat')
    hd = mul(constant(d), h)
    plus, minus = add(hd, sh, -1), add(mul(constant(-1), hd), sh, -1)
    lhs = add(mul(plus, plus), mul(minus, minus))
    rhs = add(mul(constant(2), mul(sh, sh)), mul(constant(2*d*d), mul(h, h)))
    gap = add(lhs, rhs, -1)
    if gap or not exact:
        raise ValueError('exact ambiguity certificate failed')
    return {
        'qualification': QUALIFICATION,
        'audited_sha256': hashes,
        'source_scope': 'OLD finite-window-only definition at f7123bd8; not current physical COMPLETE-BRMM admission',
        'old_source_definition': "p_dot=v; v_dot=a; S_dot=p; bounded p/v; bounded every short-window Delta S; no generator potential",
        'constant_nonzero_position_zero_velocity_history_admitted': False,
        'corrected_source_rejection': rejection,
        'old_witness_excluded_by_corrected_physical_theorem': wave['old_witness_excluded_by_corrected_physical_theorem'],
        'classification_under_intended_physical_semantics': 'E',
        'quiet_source_cell': {
            'parameter': 'd_vec in closed ball ||d_vec||<=1/8 m, retained once for the entire history',
            'p': 'd_vec', 'v': '0', 'a': '0', 'omega': '0', 'rotation': 'I',
            'S_session': 't*d_vec', 'S_L': '(t-t_L)*d_vec',
            'fresh_e_S_L': '0', 'gyro': '0', 'specific_force': '(0,0,-g)',
            'magnetic_field_uT': [20, 0, 40], 'sensor_noise': '0',
            'AC_energy_every_window': '0', 'impulse_every_window': '0',
            'gravity_direction_forcing': '0', 'P_m_m': str(pm),
            'three_second_DeltaS_norm_upper_m_s': str(3*d),
            'global_S_norm_bound': None,
        },
        'separate_bias_witnesses': separate_zero_bias_witnesses(),
        'literal_event_certificates': events,
        'parallelogram_identity_residual': encoded(gap),
        'sum_squared_S_errors': encoded(rhs),
        'max_pair_squared_error_lower_bound': encoded(mul(constant(d*d), mul(h, h))),
        'metric_consequence': 'max(W_plus,W_minus)>=m_minus*h^2*||d_vec||^2 in fixed units, even with different source-dependent metrics',
        'continuous_physical_S_energy_per_3s_word': {'h^2': str(3*d*d), 'h': str(9*d*d), '1': str(9*d*d)},
        'S_working_radius_300_unavoidable_pair_exit_after_s': str(F(300)/d),
        'general_working_radius_R_exit_after_s': 'R/||d_vec||; at least one fixed-sign history has unbounded limsup',
        'classification': 'B',
        'bounded_all18_indefinite_target_refuted_under_finite_window_definition': exact and not gap,
        'scope_limit': 'Historical B under the old definition; source-specification omission E under intended physics; the corrected generator theorem excludes this witness. No physical filter instability is asserted.',
        'conditional_ISS_with_unbounded_S_input_refuted': False,
        'nominal_filter_instability_established': False,
        'full_Normal_Live_source_admission_claimed': False,
        'full_source_transition_cover_materialized': False,
        'finite_startup_runtime_handoff_refuted': False,
        'P3_delta': 1e-18,
        'P4_PASS': False, 'P5_MAY_START': False,
        'deployment_arithmetic_qualification_claimed': False,
        'necessary_condition': 'Every indistinguishable admitted source class must have eventually bounded centered-S diameter; a uniform working tube requires a uniform diameter bound.',
    }


def validate(d):
    """Recompute the exact proof object, including coefficient maps and pins."""
    expected = build()
    return [key+' differs from recomputed certificate' for key in expected if d.get(key) != expected[key]]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d.update(validation_pass=not failures, validation_failures=failures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'classification': d['classification'], 'indefinite_target_refuted': d['bounded_all18_indefinite_target_refuted_under_finite_window_definition'], 'P4_PASS': False, 'failures': failures}, sort_keys=True))
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
