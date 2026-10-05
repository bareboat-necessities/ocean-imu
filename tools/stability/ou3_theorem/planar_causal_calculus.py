"""Analytical operation calculus for a causal planar cell, never a replay.

This module implements exact algebra used in the manuscript. It intentionally
has no time loop, observed-coefficient input, or default contraction constants.
The full shipping state and the discrete branch history remain required.
"""
from __future__ import annotations
from fractions import Fraction as F
from .matrix_certificates import add, identity, transpose
from .planar_linked_riccati_mean import product


def covariance_prediction_differential(f, p, df, dp, dq):
    return add(add(add(product(df, p, transpose(f)), product(f, dp, transpose(f))),
                   product(f, p, transpose(df))), dq)


def word_differential(factors, differentials):
    """D(L[n-1]...L[0]); preserve complete ordered suffix/prefix products."""
    if not factors or len(factors) != len(differentials):
        raise ValueError('nonempty paired factors required')
    n = len(factors[0])
    value = identity(n)
    delta = [[F(0) for _ in range(n)] for _ in range(n)]
    for a, da in zip(factors, differentials):
        delta = add(product(da, value), product(a, delta))
        value = product(a, value)
    return value, delta


def information_differential(y, sinv, dy, ds):
    """D(Y' S^-1 Y); Y=H*T for the same applied magnetic operation."""
    return add(add(product(transpose(dy), sinv, y), product(transpose(y), sinv, dy)),
               product(transpose(y), sinv, ds, sinv, y), -1)


def probe_word_differential(factors, differentials, root, droot):
    """Transport the physical two-column root, including its source-phase change."""
    value, delta = word_differential(factors, differentials)
    return product(value, root), add(product(delta, root), product(value, droot))


def information_perturbation_bound(*, center_stack_norm, stack_difference_norm,
                                   relative_innovation_radius):
    """Exact rational evaluation AFTER domain-dependent bounds are proved.

    Stack Ybar_j=Sbar_j^-1/2 Hbar_j Tbar_j and the matching differences.
    e bounds every relative innovation difference on the SAME placed window.
    Norm arguments must be rational certified upper bounds, not point norms.
    """
    b, d, e = map(F, (center_stack_norm, stack_difference_norm,
                     relative_innovation_radius))
    if min(b, d, e) < 0 or e >= 1:
        raise ValueError('nonnegative bounds and relative innovation radius <1 required')
    return (e*b*b + 2*b*d + d*d)/(1-e)


def certificate():
    return {
        'qualification': 'OU3_PLANAR_CAUSAL_CALCULUS_V1',
        'result_type': 'PROVED — analytical',
        'scope': 'operation identities and conditional cell inequalities; no invariant-cell certificate',
        'planar_state': {
            'mean_even_indices': [1, 4, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20],
            'covariance_block_dimensions': [12, 9],
            'mean_odd_zero_requires': 'planar initialization, zero lever arm/deheel, planar reference and regular parity-preserving branches',
            'retained_generators': ['raw Mahony q and integral', 'guard LP/detector/weight',
                                    'frequency/band/variance/period state', 'staged and filtered tuning',
                                    'S elapsed/period', 'AW clock/pending target',
                                    'reference acquisition/refinement/continuous HI sufficient statistics',
                                    'startup/watchdog/gravity/validity/LDLT/projection branches'],
        },
        'planar_reference_real_bound': {'norm_upper': '75', 'source_field_norm': '75',
                                        'startup_hard_iron_enabled': False,
                                        'continuous_information_gate_value': '0',
                                        'proof': 'weighted normalized pitch rotations fix e_y; (I-A^T A)e_y=0; canonical reference preserves weighted-mean norm',
                                        'reference_alignment_lower_bound': None,
                                        'float32_transfer_verified': False},
        'one_way_scope': 'Complementary tuner/S-period/AW-target only; reference/gate feedback is retained',
        'causal_center': 'same literal real-operation maps evaluated on center state and inherited same-history generator; no observed future H/G/K',
        'center_constructed_by_recurrence': True,
        'center_compact_forward_domain_certified': False,
        'gain_and_noise_calculus': 'planar-linked-riccati-mean-certificate.json',
        'mag_information_differential': "dI=sum(dY' S^-1 Y+Y' S^-1 dY-Y' S^-1 dS S^-1 Y)",
        'physical_probe_root_differential': 'dPhi=dT_word B_root+T_word dB_root; retain source phase and physical heading/BG roots',
        'mag_information_perturbation_bound': '(e*b^2+2*b*d+d^2)/(1-e)',
        'gauge_rule': 'retain C_Q alpha, covariance coupling and moving-projector differential',
        'physical_curvature_squared_bound': 'g^2*(beta^6/36+beta^4/4); nominal-chart and J bounds still required',
        'aw_rule': 'spectral divided differences away from zero; Lipschitz positive-part bound at active-face crossings',
        'projection_rule': 'piecewise derivatives on fixed faces; finite Lipschitz/remainder inclusions at boundaries',
        'source_uniform_storage_contraction': None,
        'certified_radius': None,
        'center_service_floor': None,
        'closure_status': 'OPEN',
        'missing_inputs': ['inherited startup/release cell', 'center reference alignment and nominal aw bounds',
                           'complete linked word storage inequality', 'all-prefix branch and clock containment',
                           'precision-normalized physical gauge/curvature', 'uniform center information floor'],
        'structures_preserved': ['21-state error and physical gauge', 'both actual covariance blocks',
                                 'literal prediction/Joseph/reset/AW/projection', 'same-history frontend and clocks',
                                 'actual accepted innovation-whitened magnetic probes'],
        'relaxations_introduced': ['fixed regular real-operation branches for derivatives',
                                  'explicit norm majorant only after complete linked word products'],
        'all_time_magnetic_service_verified': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
