"""Binary32 SpectralMSE ``apply_RS_tune_`` -> ``set_RS_noise`` commit graph.

For the deployed SpectralMSE law ``pseudo_update_information_rate_scale_()``
returns exactly ``1.0f``; cadence renormalization is only applied to Cubic.
Thus the shipping Live commit reduces to

    RSbase = clamp(tune_.RS_applied, min_R_S_, max_R_S_)
    s      = valid(rs_scale) ? min(rs_scale, 1.0f) : 1.0f
    rs_z   = (RSbase * 1.0f) * s
    sigma  = (rs_z*R_S_x_factor_, rs_z*R_S_y_factor_, rs_z)
    R_S    = diag(sigma.array().square())

The scalar binary32 graph is exact here.  The final Eigen expression execution
and compiler/vectorization correspondence remains a separate deployment fact;
this module does not claim bit identity for Eigen internals.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

FILTER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
MEKF=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
ONE=B.rn32(1)
DEFAULT_X_FACTOR=B.rn32(F(72,100)); DEFAULT_Y_FACTOR=B.rn32(F(72,100))
QUALIFICATION='OU3_ALT_RS_COMMIT_BINARY32_V1'


def clamp(x,lo,hi): return min(max(x,lo),hi)

def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


@dataclass(frozen=True)
class Commit:
    stored_RS:F
    min_RS:F
    max_RS:F
    rs_scale:F
    x_factor:F
    y_factor:F
    RSbase:F
    RSb:F
    rs_z:F
    sigma_x:F
    sigma_y:F
    sigma_z:F
    covariance_diag:tuple
    qualification:str=QUALIFICATION
    def __post_init__(self):
        names=('stored_RS','min_RS','max_RS','rs_scale','x_factor','y_factor','RSbase','RSb','rs_z','sigma_x','sigma_y','sigma_z')
        for n in names: object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS commit qualification')
        if not all(B.is_binary32(getattr(self,n)) for n in names): raise ValueError('RS commit stores binary32 scalars only')
        if self.max_RS<self.min_RS or self.min_RS<=0: raise ValueError('invalid RS clamp')
        if not 0<self.rs_scale<=1 or self.x_factor<0 or self.y_factor<0: raise ValueError('invalid RS commit scale/factors')
        if self.RSbase!=clamp(self.stored_RS,self.min_RS,self.max_RS): raise ValueError('RSbase detached from shipping clamp')
        if self.RSb!=B.mul(self.RSbase,ONE): raise ValueError('SpectralMSE information-rate scale must be literal one')
        if self.rs_z!=B.mul(self.RSb,self.rs_scale): raise ValueError('rs_z detached from rs_scale multiply')
        if self.sigma_x!=B.mul(self.rs_z,self.x_factor) or self.sigma_y!=B.mul(self.rs_z,self.y_factor) or self.sigma_z!=self.rs_z:
            raise ValueError('anisotropic sigma vector detached from source multiplications')
        cov=tuple(F(x) for x in self.covariance_diag); object.__setattr__(self,'covariance_diag',cov)
        expected=(B.mul(self.sigma_x,self.sigma_x),B.mul(self.sigma_y,self.sigma_y),B.mul(self.sigma_z,self.sigma_z))
        if cov!=expected: raise ValueError('R_S covariance diagonal detached from per-axis square graph')


def commit(stored_RS,*,min_RS,max_RS,rs_scale=ONE,x_factor=DEFAULT_X_FACTOR,y_factor=DEFAULT_Y_FACTOR):
    r=_q(stored_RS,'stored RS'); lo=_q(min_RS,'min RS'); hi=_q(max_RS,'max RS')
    s=_q(rs_scale,'rs_scale'); xf=_q(x_factor,'R_S_x_factor'); yf=_q(y_factor,'R_S_y_factor')
    if hi<lo or lo<=0: raise ValueError('invalid RS clamp')
    # Deployment theorem path supplies a valid positive scale; invalid caller
    # inputs take source fallback 1.0f and must be handled before this entry.
    if not 0<s<=1: raise ValueError('deployment RS scale must already be source-qualified in (0,1]')
    base=clamp(r,lo,hi)
    rsb=B.mul(base,ONE)  # explicit source operation; exact for finite normal base
    z=B.mul(rsb,s)
    sx=B.mul(z,xf); sy=B.mul(z,yf); sz=z
    cov=(B.mul(sx,sx),B.mul(sy,sy),B.mul(sz,sz))
    return Commit(r,lo,hi,s,xf,yf,base,rsb,z,sx,sy,sz,cov)


def _source_shape_matches():
    f=FILTER.read_text(); k=MEKF.read_text()
    return all(x in f for x in (
      'if (rs_law_ != RSAdaptationLaw::Cubic) return 1.0f;',
      'const float RSbase = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);',
      'const float RSb = RSbase * pseudo_update_information_rate_scale_();',
      'const float rs_z = RSb * s;',
      'rs_z * R_S_x_factor_', 'rs_z * R_S_y_factor_')) and all(x in k for x in (
      'void set_RS_noise(const Vector3& sigma_S)',
      'R_S = sigma_S.array().square().matrix().asDiagonal();'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_SpectralMSE_RS_commit_source_shape_matches':_source_shape_matches(),
      'SpectralMSE_information_rate_scale_is_literal_one':True,
      'stored_RS_clamp_binary32_materialized':True,
      'rs_scale_and_anisotropic_factor_binary32_multiplies_materialized':True,
      'per_axis_covariance_square_binary32_graph_materialized':True,
      'Eigen_set_RS_noise_execution_correspondence_closed':False,
      'compiler_vectorization_rounding_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
