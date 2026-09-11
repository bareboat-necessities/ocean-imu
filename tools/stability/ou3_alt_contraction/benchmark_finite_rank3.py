#!/usr/bin/env python3
"""Exact-arithmetic rank-3 benchmark; not a source enclosure or storage search."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import time
import tracemalloc

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as G


def solve3(S, rhs):
    rows = [list(a)+[b] for a, b in zip(S, rhs)]
    for j in range(3):
        pivot = next((i for i in range(j, 3) if rows[i][j]), None)
        if pivot is None: raise ValueError('singular regression innovation')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]; rows[j] = [v/scale for v in rows[j]]
        for i in range(3):
            if i != j:
                scale = rows[i][j]; rows[i] = [x-scale*y for x, y in zip(rows[i], rows[j])]
    return [row[-1] for row in rows]


def measure(function, repeat):
    start = time.perf_counter()
    for _ in range(repeat): value = function()
    seconds = (time.perf_counter()-start)/repeat
    tracemalloc.start()
    again = function()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    if again != value: raise AssertionError('non-deterministic exact arithmetic')
    return value, {'seconds_per_call': seconds, 'python_peak_bytes': peak}


def benchmark(repeat=3):
    if repeat < 1: raise ValueError('repeat must be positive')
    reports = {}
    for active in (False, True):
        u = [F((i % 4)-2, 20) for i in range(21)]
        P = G.plus(G.eye(21), [[x*y for y in u] for x in u])
        H, _ = G.residual_factors('accelerometer', [0, 0, 0],
                                  f_hat=[F(1, 5), F(-1, 4), -9], R_hat=G.eye(3))
        N, S = G.measurement_operands(P, G.eye(3), H, active_bias=active)
        K = [solve3(S, row) for row in N]
        def raw():
            return G.plus(G.plus(G.plus(P, G.mm(K, G.transpose(N)), -1),
                          G.mm(N, G.transpose(K)), -1), G.mm(G.mm(K, S), G.transpose(K)))
        direct, r = measure(raw, repeat)
        thin, t = measure(lambda: G.solved_joseph_covariance(P, N, S, K), repeat)
        if thin != direct: raise AssertionError('rank-3 rewrite changed the exact covariance')
        reports['A21' if active else 'H18'] = {
            'raw_Joseph': r, 'solve_constrained_rank3': t,
            'speedup': r['seconds_per_call']/t['seconds_per_call'],
            'exact_max_entry_gap': '0',
            'same_N_Sigma_K_retained': True, 'independent_interval_boxes_created': False,
            'source_uniform_enclosure_benchmark_performed': False,
            'subdivision_or_certificate_margin_improvement_claimed': False}
    return {'qualification': 'EXACT_ALGEBRA_ARITHMETIC_BENCHMARK_ONLY', 'modes': reports,
            'source_uniform_word_qualified': False, 'ALT_LIVE_PASS': False, 'ALT_STARTUP_PASS': False,
              'ALT_END_TO_END_PASS': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeat', type=int, default=3)
    args = parser.parse_args(); result = benchmark(args.repeat)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2))
