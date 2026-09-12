"""Finite MagAutoTuner mean -> gauge-fixed world-reference boundary.

After the accepted-window accumulator has produced its raw weighted world mean
m, shipping removes the arbitrary horizontal yaw gauge by storing

    B_ref = (sqrt(m_x^2 + m_y^2), 0, m_z).

This module closes that algebraic boundary and its positive norm readiness gate.
It does NOT qualify the accepted-window mean or the binary32 sqrt itself; the
sqrt is a same-expression witness until native correspondence is attached.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


def R(x): return M.rational(x)

@dataclass(frozen=True)
class HorizontalSqrt:
    radicand:F
    value:F
    def __post_init__(self):
        a,v=R(self.radicand),R(self.value)
        if a<0 or v<0 or v*v!=a: raise ValueError('exact nonnegative horizontal sqrt witness required')
        object.__setattr__(self,'radicand',a); object.__setattr__(self,'value',v)

@dataclass(frozen=True)
class Result:
    mean:tuple
    world_reference:tuple
    ready:bool


def gauge_fix(mean,*,mag_norm_min,sqrt:HorizontalSqrt):
    m=tuple(M.vec(mean,3)); mn=R(mag_norm_min)
    if mn<0: raise ValueError('nonnegative magnetic norm threshold required')
    h2=m[0]*m[0]+m[1]*m[1]
    if not isinstance(sqrt,HorizontalSqrt) or sqrt.radicand!=h2:
        raise ValueError('horizontal sqrt witness detached from same magnetic mean')
    ref=(sqrt.value,F(0),m[2])
    ready=M.dot(ref,ref)>mn*mn
    return Result(m,ref,ready)


def readiness():
    return {
      'MagAutoTuner_horizontal_gauge_fix_materialized':True,
      'reference_y_component_exactly_zero':True,
      'ready_norm_threshold_materialized_without_free_reference':True,
      'accepted_window_mean_runtime_ancestry_attached':False,
      'horizontal_sqrt_binary32_attached':False,
      'finite_allFinite_branch_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
