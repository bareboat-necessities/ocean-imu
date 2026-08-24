#!/usr/bin/env python3
"""Compatibility entry point for the OU-III P3 information certificate.

The scalar physical-unit min(Q)/trace(P) implementation has been retired.
All callers now execute the source-reachable matrix-valued backend in
``ou3_source_reachable_matrix_p3``.
"""
from ou3_source_reachable_matrix_p3 import *  # noqa: F401,F403


if __name__ == "__main__":
    raise SystemExit(main())
