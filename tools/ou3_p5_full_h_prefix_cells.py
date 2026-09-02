#!/usr/bin/env python3
"""Compatibility aliases for active P4 covariance primitives.

The former full-H P5 prefix certificate/search implementation has been retired.
A few current P4 modules historically imported covariance helpers from this
module; only those aliases remain until their imports are renamed.
"""
from __future__ import annotations

from ou3_interval import matrix_identity
import ou3_p4_covariance_primitives as COV

MEKF = COV.MEKF
N = 18
TH = range(0, 3)
BG = range(3, 6)
V = range(6, 9)
P = range(9, 12)
SS = range(12, 15)
AW = range(15, 18)
BA = range(18, 18)

I = COV.I
down = COV.down
up = COV.up
_zero = COV.zero
_psd_tighten = COV.psd_tighten
_reset_covariance = COV.reset_covariance
_R_diag = COV.R_diag
_R_S = COV.R_S
_source_pb0 = COV.source_pb0
_rotation_step_box = COV.rotation_step_box
