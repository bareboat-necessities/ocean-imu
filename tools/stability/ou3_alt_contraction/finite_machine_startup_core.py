"""Source-owned default construction and ungauged binary32 CORE handoff.

The selected wrapper/application preserves the default initial gyro-bias,
linear and accel-bias covariance settings. Its bootstrap does not drive the
MEKF. The a_w construction block is irrelevant to this dependency slice:
goLive overwrites its entire rows/columns using the produced machine commit.
Noise settings need not equal the wrapper defaults for this argument.

No q/P result is an input. The ungauged seed is the actual carried private
observer quaternion. The two setter normalizations and down-axis covariance
normalization/fallback execute in a named scalar Eigen arithmetic profile.
The operation schedule is explicit; target compilation correspondence and the
machine magnetic pending-yaw producer remain separate obligations.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import re

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SQRT
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as FRESH
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as ENTRY
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE

PROFILE='OU3_ALT_DEFAULT_CONSTRUCTION_UNGAUGED_SCALAR_EIGEN_V1'
_TOKEN=object()
ROOT=Path(__file__).resolve().parents[3]


@lru_cache(maxsize=1)
def source_audit():
    paths=('src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h',
           'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h',
           'src/kalman_common/SeaStateFusionFilterCommon.h',
           'sensors/full_marine_ins/atomS3R_ins_kalman_ou3/atomS3R_ins_kalman_ou3.ino')
    texts=[(ROOT/p).read_text() for p in paths]
    wrapper,mekf,common,app=texts
    required=(
        (wrapper,('float Pq0       = 5e-4f;','float Pb0       = 1e-6f;',
            'cfg_.Pq0, cfg_.Pb0, cfg_.b0, cfg_.R_S_noise,',
            'impl_.updateFrontEnd(dt, gyro_body_ned, acc_body_ned);',
            'mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);',
            'mekf_->reset_aw_covariance_to_stationary();')),
        (mekf,('xext.setZero();','Pext.setZero();','T sigma_bacc0_ = T(0.004);',
            'const T sigma_v0 = T(1.0);','const T sigma_p0 = T(20.0);',
            'const T sigma_S0 = T(50.0);','set_quaternion_boat(q_bw);',
            'q.normalize();','qref.normalize();','u_down_body /= n;',
            'Pbase.template block<3,3>(3,3) = Matrix3::Identity() * Pb0;',
            'set_initial_linear_uncertainty(sigma_v0, sigma_p0, sigma_S0);')),
        (common,('enterCold();','applyTune();','mekf->reset_aw_covariance_to_stationary();')),
        (app,('Fusion::Config fcfg;','fusion_.begin(fcfg);')))
    for text,needles in required:
        if any(n not in text for n in needles):
            raise RuntimeError('default startup CORE source shape changed')
    if re.search(r'fcfg\.(?:Pq0|Pb0)\s*=',app) or any(name in app for name in
        ('set_initial_linear_uncertainty(', 'set_initial_acc_bias_std(', 'set_initial_acc_bias(')):
        raise RuntimeError('application overrides retained startup construction state')
    return {p:hashlib.sha256(t.encode()).hexdigest() for p,t in zip(paths,texts)}


@dataclass(frozen=True,init=False)
class Construction:
    """Sealed default retained-state constructor, not a covariance snapshot."""
    profile:str
    source_hashes:tuple
    Pq0:F
    Pb0:F
    ba_variance:F
    def __init__(self,*,_token=None):
        if _token is not _TOKEN:
            raise TypeError('default construction is produced only by literal_reset')
        object.__setattr__(self,'profile',PROFILE)
        object.__setattr__(self,'source_hashes',tuple(source_audit().items()))
        object.__setattr__(self,'Pq0',B.rn32(F(5,10000)))
        object.__setattr__(self,'Pb0',B.rn32(F(1,10**6)))
        ba=B.rn32(F(1,250)); object.__setattr__(self,'ba_variance',B.mul(ba,ba))


def literal_reset(): return Construction(_token=_TOKEN)


@dataclass(frozen=True)
class Operation:
    kind:str
    operands:tuple
    result:F


class Arithmetic:
    def __init__(self,mode):
        if mode not in ('separate','fma'): raise ValueError('literal compiler arithmetic mode required')
        self.mode=mode; self.operations=[]
    def op(self,kind,*xs):
        if kind=='sqrt': out=SQRT.sqrt32(xs[0])
        elif kind=='add': out=B.add(*xs)
        elif kind=='sub': out=B.sub(*xs)
        elif kind=='mul': out=B.mul(*xs)
        elif kind=='div': out=B.div(*xs)
        elif kind=='fma': out=B.fma(*xs)
        else: raise ValueError('unknown initialization operation')
        self.operations.append(Operation(kind,tuple(xs),out)); return out
    def norm(self,xs):
        # Eigen's fixed-size scalar reduction recursively splits at n//2.
        def sumsq(v):
            if len(v)==1: return self.op('mul',v[0],v[0])
            k=len(v)//2
            if len(v)==2 and self.mode=='fma':
                return self.op('fma',v[0],v[0],self.op('mul',v[1],v[1]))
            return self.op('add',sumsq(v[:k]),sumsq(v[k:]))
        return self.op('sqrt',sumsq(tuple(xs)))
    def qnorm(self,q):
        # Quaternion::coeffs storage order is x,y,z,w.
        return self.norm((q[1],q[2],q[3],q[0]))
    def normalized_q(self,q):
        n=self.qnorm(q)
        if not n> B.rn32(F(1,10**8)):
            raise ValueError('shipping quaternion setter rejects tiny seed')
        return tuple(self.op('div',x,n) for x in q)


def _attitude_covariance(a,q,tilt_sigma,yaw_sigma):
    op=a.op; w,x,y,z=q
    tx=op('mul',2,x); ty=op('mul',2,y); tz=op('mul',2,z)
    # The third column of Eigen Quaternion::toRotationMatrix().
    u=(op('add',op('mul',tz,x),op('mul',ty,w)),
       op('sub',op('mul',tz,y),op('mul',tx,w)),
       op('sub',1,op('add',op('mul',tx,x),op('mul',ty,y))))
    n=a.norm(u); tv=op('mul',tilt_sigma,tilt_sigma); yv=op('mul',yaw_sigma,yaw_sigma)
    if not n>B.rn32(F(1,10**8)):
        return tuple(tuple(yv if i==j else F(0) for j in range(3)) for i in range(3)),True
    u=tuple(op('div',v,n) for v in u)
    cov=[]
    for i in range(3):
        row=[]
        for j in range(3):
            yy=op('mul',u[i],u[j]); tt=op('sub',int(i==j),yy)
            ypart=op('mul',yv,yy)
            row.append(op('fma',tv,tt,ypart) if a.mode=='fma'
                       else op('add',op('mul',tv,tt),ypart))
        cov.append(row)
    return tuple(tuple(op('mul',F(1,2),op('add',cov[i][j],cov[j][i]))
                       for j in range(3)) for i in range(3)),False


@dataclass(frozen=True,init=False)
class Result:
    construction:Construction
    source:object
    mode:str
    state:CORE.State
    operations:tuple
    covariance_fallback:bool
    profile:str
    def __init__(self,construction,source,mode,state,operations,fallback,*,_token=None):
        if _token is not _TOKEN:
            raise TypeError('machine initialization result must be executed from startup')
        for k,v in (('construction',construction),('source',source),('mode',mode),('state',state),
                    ('operations',tuple(operations)),('covariance_fallback',fallback),('profile',PROFILE)):
            object.__setattr__(self,k,v)


def ungauged(startup,reference,active,*,mode):
    """Execute source-owned default initialization without a q/P input port."""
    from tools.stability.ou3_alt_contraction import finite_startup_joined_machine_history as START
    if not isinstance(startup,START.State) or not isinstance(startup.core_construction,Construction):
        raise TypeError('literal default construction ancestry required')
    if not isinstance(reference,CORE.Reference): raise TypeError('same physical Live reference required')
    source=getattr(startup,mode+'_source',None)
    if source is None or not source.vertical.initialized:
        raise ValueError('ungauged handoff requires actual initialized private observer')
    if startup.last_raw is None or any(getattr(startup.last_raw.physical,k)!=getattr(reference,k)
        for k in ('history_id','bias_root','bias_family')):
        raise ValueError('machine handoff detached from startup physical/bias lineage')
    a=Arithmetic(mode); q=tuple(source.vertical.q)
    if any(not B.is_binary32(x) for x in q): raise ValueError('startup proxy is not actual binary32')
    # zero heel: multiplication by the exact identity does not change finite
    # values; the setter still normalizes both before and after conjugation.
    first=a.normalized_q(q)
    qhat=a.normalized_q((first[0],-first[1],-first[2],-first[3]))
    tilt=B.rn32(F(35,1000)); yaw=B.rn32(F(15708,10000))
    att,fallback=_attitude_covariance(a,qhat,tilt,yaw)
    c=startup.core_construction
    cov=[[F(0)]*21 for _ in range(21)]
    for i in range(3):
        for j in range(3): cov[i][j]=att[i][j]
        cov[3+i][3+i]=c.Pb0
        cov[6+i][6+i]=F(1); cov[9+i][9+i]=F(400); cov[12+i][12+i]=F(2500)
        cov[18+i][18+i]=c.ba_variance
        for j in range(3): cov[15+i][15+j]=active.Sigma_aw[i][j]
    # All nonattitude means remain literal constructor zero through bootstrap;
    # initialization clears attitude error and preserves the closed H18 gate.
    entry=ENTRY.Result((F(0),)*21,tuple(tuple(row) for row in cov),qhat,active,
                       'Live',F(0),False,False,True)
    state=FRESH.build(entry,reference,scope=SCOPE.certified_scope())
    return Result(c,source,mode,state,a.operations,fallback,_token=_TOKEN)


def readiness():
    source_audit()
    return {'default_retained_construction_from_selected_wrapper_and_app':True,
            'bootstrap_zero_mean_and_non_aw_covariance_preserved':True,
            'ungauged_machine_proxy_to_full_CORE_producer_available':True,
            'caller_replacement_quaternion_or_covariance_port':False,
            'normalization_and_down_covariance_fallback_executed':True,
            'gauged_machine_pending_yaw_producer_available':False,
            'target_Eigen_compiler_operation_schedule_qualified':False,
            'universal_machine_handoff_quaternion_lower_norm_proved':False,
            'source_uniform_startup_CORE_qualified':False,
            'source_uniform_complete_600_step_word_qualified':False}
