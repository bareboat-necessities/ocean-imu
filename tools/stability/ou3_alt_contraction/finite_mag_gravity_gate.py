"""Startup magnetic gravity-alignment and admission gate for ALT.

Shipping advances the gravity certificate on IMU samples, not updateMag calls:

  a_world = q_proxy * a_body
  LPF(a_world; tau_world)
  warm = elapsed >= warmup
  align_sin = |LPF_xy|/|LPF| when warm, else 1
  aligned_branch = LPF.z < 0
  extreme_motion = |gyro|*rad2deg > threshold
  good += dt (cap 10 s) when all gates pass, else good -= 2 dt (floor 0)

An asynchronous startup updateMag call is admitted only after with-mag/delay and
settle gates.  The first call reaching that point latches ``eligible_t0`` even
when gravity is not yet trusted.  Accumulation may begin when either the gravity
hold is satisfied or the fallback duration from that same t0 has elapsed, and
only when a previous IMU packet exists.

The LPF exponential and norm/sqrt deployment arithmetic remain explicit
same-expression witnesses.  This is control/event identity, not source-bound or
startup-capture qualification.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


def R(x): return M.rational(x)
RAD2DEG=F(57295779513,1000000000)

@dataclass(frozen=True)
class SqrtWitness:
    radicand:F
    value:F
    def __post_init__(self):
        a,v=R(self.radicand),R(self.value)
        if a<0 or v<0 or v*v!=a: raise ValueError('exact nonnegative sqrt witness required')
        object.__setattr__(self,'radicand',a); object.__setattr__(self,'value',v)

@dataclass(frozen=True)
class LpfExpWitness:
    """Same-expression witness for exp(-dt/max(1e-3,tau))."""
    dt:F
    tau_used:F
    exp_value:F
    def __post_init__(self):
        d,t,e=R(self.dt),R(self.tau_used),R(self.exp_value)
        if d<=0 or t<F(1,1000) or e<0 or e>1:
            raise ValueError('valid LPF exponential witness required')
        object.__setattr__(self,'dt',d); object.__setattr__(self,'tau_used',t); object.__setattr__(self,'exp_value',e)
    @property
    def alpha(self): return F(1)-self.exp_value

@dataclass(frozen=True)
class Config:
    world_tau:F=F(12)
    world_warmup:F=F(5)
    max_align_sin:F=F(75,1000)
    hold_sec:F=F(2)
    fallback_sec:F=F(30)
    extreme_gyro_dps:F=F(30)
    mag_delay:F=F(7)
    proxy_settle:F=F(0)
    with_mag:bool=True
    def __post_init__(self):
        for n in ('world_tau','world_warmup','max_align_sin','hold_sec','fallback_sec','extreme_gyro_dps','mag_delay','proxy_settle'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative startup mag gate config required')
            object.__setattr__(self,n,v)
        if not isinstance(self.with_mag,bool): raise TypeError('literal with_mag branch required')

@dataclass(frozen=True)
class State:
    lpf:tuple=(F(0),F(0),F(0))
    lpf_initialized:bool=False
    world_elapsed:F=F(0)
    gravity_good:F=F(0)
    aligned_branch:bool=False
    eligible_t0:F|None=None
    def __post_init__(self):
        object.__setattr__(self,'lpf',tuple(M.vec(self.lpf,3)))
        for n in ('world_elapsed','gravity_good'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative gate clock required')
            object.__setattr__(self,n,v)
        if not isinstance(self.lpf_initialized,bool) or not isinstance(self.aligned_branch,bool):
            raise TypeError('literal LPF/branch flags required')
        if self.eligible_t0 is not None:
            t=R(self.eligible_t0)
            if t<0: raise ValueError('nonnegative eligibility origin required')
            object.__setattr__(self,'eligible_t0',t)

@dataclass(frozen=True)
class ImuResult:
    state:State
    world_accel:tuple
    align_sin:F
    gyro_dps:F
    gravity_good_now:bool

@dataclass(frozen=True)
class AdmissionResult:
    state:State
    admitted:bool
    gravity_trusted:bool
    fallback_ok:bool
    reason:str


def imu_step(state:State,cfg:Config,*,q_proxy_bw,acc_body,gyro_body,dt,
             lpf_exp:LpfExpWitness|None,
             lpf_norm:SqrtWitness|None=None,horizontal_norm:SqrtWitness|None=None,
             gyro_norm:SqrtWitness|None=None):
    if not isinstance(state,State) or not isinstance(cfg,Config): raise TypeError('gravity gate state/config required')
    q=tuple(M.vec(q_proxy_bw,4)); a=tuple(M.vec(acc_body,3)); g=tuple(M.vec(gyro_body,3)); d=R(dt)
    if d<=0: raise ValueError('positive IMU dt required')
    q2=M.dot(q,q)
    if q2==0: raise ValueError('zero startup proxy quaternion')
    # SENSOR.q_rotate is homogeneous in q and therefore implements the same
    # normalized quaternion rotation without a detached normalization operand.
    world=tuple(SENSOR.q_rotate(q,a))
    tau=max(F(1,1000),cfg.world_tau)
    if state.lpf_initialized:
        if not isinstance(lpf_exp,LpfExpWitness) or lpf_exp.dt!=d or lpf_exp.tau_used!=tau:
            raise ValueError('LPF exponential witness detached from same dt/tau')
        alpha=lpf_exp.alpha
        lp=tuple(state.lpf[i]+alpha*(world[i]-state.lpf[i]) for i in range(3))
    else:
        if lpf_exp is not None: raise ValueError('first LPF sample consumes no exponential witness')
        lp=world
    elapsed=state.world_elapsed+d
    warm=elapsed>=cfg.world_warmup
    aligned=(M.dot(lp,lp)>F(1,10**12) and lp[2]<0)
    if warm:
        n2=M.dot(lp,lp); h2=lp[0]*lp[0]+lp[1]*lp[1]
        if not isinstance(lpf_norm,SqrtWitness) or lpf_norm.radicand!=n2:
            raise ValueError('LPF norm witness detached from same world average')
        if lpf_norm.value<=F(1,10**6):
            align=F(1)
        else:
            if not isinstance(horizontal_norm,SqrtWitness) or horizontal_norm.radicand!=h2:
                raise ValueError('horizontal norm witness detached from same world average')
            align=min(F(1),max(F(0),horizontal_norm.value/lpf_norm.value))
    else:
        if lpf_norm is not None or horizontal_norm is not None:
            raise ValueError('pre-warm gravity gate consumes no alignment norm witnesses')
        align=F(1)
    g2=M.dot(g,g)
    if not isinstance(gyro_norm,SqrtWitness) or gyro_norm.radicand!=g2:
        raise ValueError('gyro norm witness detached from same IMU packet')
    gyro_dps=gyro_norm.value*RAD2DEG
    extreme=gyro_dps>cfg.extreme_gyro_dps
    good_now=(align<=cfg.max_align_sin and aligned and not extreme)
    good=min(F(10),state.gravity_good+d) if good_now else max(F(0),state.gravity_good-2*d)
    nxt=State(lp,True,elapsed,good,aligned,state.eligible_t0)
    return ImuResult(nxt,world,align,gyro_dps,good_now)


def startup_mag_admission(state:State,cfg:Config,*,wrapper_time,begun,have_last_imu,mag_ref_set=False):
    if not isinstance(state,State) or not isinstance(cfg,Config): raise TypeError('gravity gate state/config required')
    if not isinstance(begun,bool) or not isinstance(have_last_imu,bool) or not isinstance(mag_ref_set,bool):
        raise TypeError('literal wrapper control flags required')
    t=R(wrapper_time)
    if t<0: raise ValueError('nonnegative wrapper time required')
    if not begun: return AdmissionResult(state,False,False,False,'not_begun')
    if not cfg.with_mag: return AdmissionResult(state,False,False,False,'mag_disabled')
    if t<cfg.mag_delay: return AdmissionResult(state,False,False,False,'mag_delay')
    # Initial-reference acquisition only.  Once a reference exists, refinement
    # and continuous hard-iron use separate event machines.
    if mag_ref_set: return AdmissionResult(state,False,False,False,'reference_already_set')
    if t<cfg.proxy_settle: return AdmissionResult(state,False,False,False,'proxy_settle')
    t0=t if state.eligible_t0 is None else state.eligible_t0
    latched=State(state.lpf,state.lpf_initialized,state.world_elapsed,state.gravity_good,state.aligned_branch,t0)
    gravity_trusted=latched.gravity_good>=cfg.hold_sec
    fallback=(t-t0)>=cfg.fallback_sec
    if not gravity_trusted and not fallback:
        return AdmissionResult(latched,False,gravity_trusted,fallback,'gravity_or_fallback_wait')
    if not have_last_imu:
        return AdmissionResult(latched,False,gravity_trusted,fallback,'no_last_imu')
    return AdmissionResult(latched,True,gravity_trusted,fallback,'admitted')


def readiness():
    return {
      'world_acceleration_uses_same_startup_proxy_attitude':True,
      'world_frame_LPF_state_and_elapsed_persisted':True,
      'LPF_first_sample_and_exponential_update_materialized':True,
      'warmup_align_sin_and_right_side_up_branch_materialized':True,
      'extreme_gyro_veto_materialized':True,
      'gravity_good_hold_cap_and_2dt_recovery_materialized':True,
      'mag_delay_settle_and_eligible_t0_latch_materialized':True,
      'gravity_hold_or_fallback_admission_materialized':True,
      'have_last_imu_required_for_initial_accumulation':True,
      'LPF_exp_norm_binary32_attached':False,
      'raw_IMU_source_bounds_attached':False,
      'async_mag_call_schedule_source_attached':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
