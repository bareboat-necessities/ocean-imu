#!/usr/bin/env python3
"""Full-scale P4 working set plus qualified fresh-Live reachability.

The coordinate radii remain the existing deterministic P4 working scales, not a
covariance-confidence set.  Fresh reachability is now supplied separately by
the uniformly qualified BRMM physical primitive/source graph: zero shipping
linear means, exact shared-S re-anchoring, and hard physical source bounds place
the fresh outer-wrapper error inside these same full-scale radii.
"""
from __future__ import annotations
import argparse,json,math
import ou3_p4_live_entry_graph as LIVE_ENTRY
import ou3_p4_qualified_fresh_entry as QUALIFIED_ENTRY
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'

def build(path:Path=CLOSURE):
    c=json.loads(Path(path).read_text())['hard_entry_search'];scale=float(c['minimum_certified_scale']);b={k:float(v) for k,v in c['base_coordinate_radii'].items()}
    if scale!=1.0 or list(map(float,c['candidate_scale_factors']))!=[1.0:]:raise RuntimeError('unreachable')
