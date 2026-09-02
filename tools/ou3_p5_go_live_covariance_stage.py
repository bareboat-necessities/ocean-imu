#!/usr/bin/env python3
"""Compatibility shim for the active P4 goLive covariance seed.

The historical P5 pre-first-S search/certificate logic has been removed.
"""
from ou3_p4_golive_covariance import DEFAULT_DOMAIN, build, validate
