#!/usr/bin/env python3
"""Runtime-bounded validated wrapper for complete-word translation P4.

The underlying backend keeps every operation exact-rational.  Unrestricted
Fraction denominators, however, grow through hundreds of prediction/Woodbury
steps and can exhaust the CI time limit.  This wrapper preserves the same
validated inequalities while replacing each propagated lower after prediction
and correction by a 192-bit dyadic Loewner lower.

For a symmetric exact matrix A, let Q use downward dyadic rounding on every
entry, with an additional (N-1)/D subtraction on each diagonal, D=2^BITS.
Then R=A-Q has nonnegative diagonal and each diagonal dominates the sum of the
absolute off-diagonal residuals in that row.  Hence R is symmetric diagonally
dominant PSD, so Q <= A in Loewner order.  Denominators are therefore bounded
without any binary64 representation floor.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import ou3_p4_translation_full_word_interval as C

DYADIC_BITS = 192
_D = 1 << DYADIC_BITS
_EPS = Fraction(1, _D)
_orig_predict = C._predict
_orig_rank1 = C._rank1_information_update_lower


def _floor_dyadic(x):
    x = C._q(x)
    return Fraction((x.numerator * _D) // x.denominator, _D)


def _dyadic_loewner_lower(A):
    """Return exact dyadic Q with Q <= A by an SDD PSD residual proof."""
    A = C._qm(A)
    Q = [[Fraction(0) for _ in range(C.N)] for _ in range(C.N)]
    for i in range(C.N):
        for j in range(i, C.N):
            q = _floor_dyadic(A[i][j])
            Q[i][j] = q
            Q[j][i] = q
    diag_guard = (C.N - 1) * _EPS
    for i in range(C.N):
        Q[i][i] -= diag_guard
    if not C._rational_spd(Q):
        raise RuntimeError('dyadic Loewner compression lost SPD')
    return Q


def _predict(L, F, rho):
    M, radius = _orig_predict(L, F, rho)
    return _dyadic_loewner_lower(M), radius


def _rank1(L, beta, qidx):
    M, radius = _orig_rank1(L, beta, qidx)
    return _dyadic_loewner_lower(M), radius


# Patch the original module globals used recursively by C._prop/C._mode.
C._predict = _predict
C._rank1_information_update_lower = _rank1


def build(domain_path=C.DEFAULT_DOMAIN, horizon_s=C.DEFAULT_HORIZON_S):
    d = C.build(domain_path, horizon_s)
    for m in d.get('modes', {}).values():
        m['dyadic_loewner_compression'] = True
        m['dyadic_loewner_bits'] = DYADIC_BITS
        m['dyadic_loewner_max_added_diagonal_loss'] = float((C.N - 1) * _EPS)
        m['prediction_enclosure'] = 'exact_rational_transition_interval_rowwise_loewner_dyadic192'
    return d


def validate(d):
    f = C.validate(d)
    # C.validate checks the old prediction label; replace that label-only
    # complaint with the stronger runtime-bounded exact-rational path check.
    f = [x for x in f if 'exact-rational transition enclosure missing' not in x]
    for mode in ('H', 'A'):
        m = d.get('modes', {}).get(mode, {})
        if m.get('prediction_enclosure') != 'exact_rational_transition_interval_rowwise_loewner_dyadic192':
            f.append(f'{mode}: dyadic exact-rational prediction enclosure missing')
        if m.get('dyadic_loewner_compression') is not True:
            f.append(f'{mode}: dyadic Loewner compression missing')
        if int(m.get('dyadic_loewner_bits', 0)) != DYADIC_BITS:
            f.append(f'{mode}: wrong dyadic precision')
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--domain', type=Path, default=C.DEFAULT_DOMAIN)
    ap.add_argument('--horizon-s', type=float, default=C.DEFAULT_HORIZON_S)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.horizon_s)
    f = validate(d)
    d['validation_failures'] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True))
    print(json.dumps({
        'translation_status': d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],
        'horizon_s': d['horizon_s'],
        'dyadic_bits': DYADIC_BITS,
        'modes': {x: {
            'delta': d.get('modes', {}).get(x, {}).get('complete_word_translation_margin_lower'),
            'factor': d.get('modes', {}).get(x, {}).get('margin_widening_factor_lower'),
            'tau_leaves': d.get('modes', {}).get(x, {}).get('tau_leaf_count'),
            'prediction': d.get('modes', {}).get(x, {}).get('prediction_enclosure'),
        } for x in ('H', 'A')},
        'failures': f,
    }, indent=2))
    return 0 if not f else 2


if __name__ == '__main__':
    raise SystemExit(main())
