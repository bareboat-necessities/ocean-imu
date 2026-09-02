#!/usr/bin/env python3
"""Temporary compatibility adapter for active P4 covariance primitives.

The historical P5 prefix implementation is retired.  Current P4 modules still
reference two helper names from its v2 module; both are forwarded directly to
the P4 covariance backend.  No P5 certificate/search logic lives here.
"""
from __future__ import annotations

from pathlib import Path

import ou3_p4_covariance_primitives as COV
import ou3_p5_full_h_prefix_cells as V1


def _tight_transition_and_Q(src: dict, domain: dict):
    return COV.tight_transition_and_Q(V1.N, src, domain)


def _corrected_initial_covariance(src: dict, domain_path: Path):
    return COV.initial_covariance(V1.N, src, domain_path)
