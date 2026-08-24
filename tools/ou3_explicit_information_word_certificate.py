#!/usr/bin/env python3
"""Compatibility entry point for the OU-III P3 information certificate.

The scalar physical-unit min(Q)/trace(P) implementation has been retired.
All callers execute the source-reachable matrix-valued backend with conditioned
source-equivalent process formulas on both sides of the x=0.01 branch.
"""
from ou3_source_reachable_matrix_p3_factored import *  # noqa: F401,F403


if __name__ == "__main__":
    raise SystemExit(main())
