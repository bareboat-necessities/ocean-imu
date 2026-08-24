#!/usr/bin/env python3
"""Compatibility entry point for the OU-III P3 information certificate.

The scalar physical-unit min(Q)/trace(P) path and the later scalar generalized
rho/max(Sigma) reduction are both retired.  All callers execute the
source-reachable direct matrix backend, with exact RL^{-1} congruence applied to
both Omega_word and Sigma_upper before validated interval LDLT.
"""
from ou3_source_reachable_matrix_p3_direct import *  # noqa: F401,F403


if __name__ == "__main__":
    raise SystemExit(main())
