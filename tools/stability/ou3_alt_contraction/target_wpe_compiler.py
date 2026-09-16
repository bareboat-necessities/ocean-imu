"""Select the pinned MCU WPE operation graph inside its persistent FMA track.

The target has a fixed mixed graph: separate HPF products, left-product FMA
for integrators/moments, MSUB for variances and omega, and FMA for log EMA.
The moment model already contains these same-operand choices. This module
constructs that actual choice on every reached branch, then calls the existing
machine relation to check its inclusion. No source coefficient is selected
independently and no histories are recombined.

The target object/assembly audit pins the compiled update method; its use by
the final firmware remains the separate whole-program compiler/link premise.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as L
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as X
from tools.stability.ou3_alt_contraction import finite_wpe_usable_binary32 as US
from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as U
from tools.stability.ou3_alt_contraction import target_wpe_libm as LIBM

QUALIFICATION='OU3_ALT_PINNED_WPE_COMPILER_FMA_TRACK_V1'
ROOT=Path(__file__).resolve().parents[3]
SOURCE_SHA='91c15942f24e08f1c39fa131f632d0974e6ab536a93a29a5af7650024ee41b7e'
LIMITS_SHA='91129a09c6fe328380410fee927936eeab19d695428328d7825aed06ca082abb'
PROBE='''#include "tuner/WavePeriodEstimator.h"
extern "C" void alt_wpe(WavePeriodEstimator* p,float dt,float a){p->update(dt,a);}
'''
CFG=M.Config(U.LAMBDA,F(4),F(20),F(180))
# Source-line compiler mapping. Each sequence retains the actual operand order.
# Register allocation is recorded in the audit as well; these operation lists
# identify the source arithmetic sites independently of register spelling.
EXPECTED={120:('add.s','sub.s','mul.s'),122:('add.s','sub.s','mul.s'),
  127:('mul.s','madd.s'),128:('mul.s','madd.s'),130:('add.s',),
  145:('sub.s','madd.s'),146:('mul.s','madd.s'),147:('mul.s','madd.s'),
  148:('mul.s','madd.s'),149:('mul.s','madd.s'),156:('msub.s',),
  158:('msub.s',),165:('msub.s',),222:('mul.s',),232:('sub.s',),
  258:('mul.s',),263:('sub.s','madd.s')}


@dataclass(frozen=True)
class State:
    moment:M.State=M.State()
    log:L.Track=L.Track()
    usable:bool=False
    def __post_init__(self):
        if not isinstance(self.moment,M.State) or not isinstance(self.log,L.Track):
            raise TypeError('persistent target WPE moment/log states required')
        if type(self.usable) is not bool:
            raise TypeError('literal target WPE usability state required')


@dataclass(frozen=True)
class Step:
    before:State
    state:State
    moment:M.StepResult
    log_step:object
    usable_step:US.Decision
    witnesses:X.ModeWitnesses
    dt:F
    vertical_accel:F


def initial():
    # The model's None diagnostic horizon means "not computed since reset";
    # shipping resets its diagnostic-only last horizon to min_horizon. This
    # slot is never read by any WPE transition. Every active state agrees.
    return State()


def step(before:State,*,dt,vertical_accel,libm):
    """Construct the target transition from SAME-argument library outputs.

    libm(kind, argument) supplies the actual exp/log/sqrt call result. Every
    argument is calculated here from this one history and checked against the
    qualified ERROR_PROFILE; repeated deterministic calls retain bit identity.
    Returns an existing ModeWitnesses object for the persistent .fma track.
    """
    if not isinstance(before,State): raise TypeError('target WPE State required')
    h,x=M._q(dt,'target WPE dt'),M._q(vertical_accel,'target WPE vertical input')
    if h!=U.DT or abs(x)>512:
        raise ValueError('target WPE step outside canonical clock/MEMS input profile')
    cache={}
    def call(kind,arg):
        arg=F(arg); value=M._q(libm(kind,arg),'target '+kind+' result')
        key=(kind,arg)
        if key in cache and cache[key]!=value:
            raise ValueError('deterministic target libm repeated argument changed bits')
        cache[key]=value
        if kind=='log':
            U.check_log(arg,value,scale_exponent=4,libm_profile=U.ERROR_PROFILE)
        else:
            {'exp':U.check_exp,'sqrt':U.check_sqrt}[kind](arg,value,libm_profile=U.ERROR_PROFILE)
        return value
    s=before.moment
    d=call('exp',-B.mul(CFG.lambda_,h))
    gain=B.div(B.sub(1,d),CFG.lambda_)
    hp1=B.mul(d,B.sub(B.add(s.hp1,x),s.accel_prev))
    hp2=B.mul(d,B.sub(B.add(s.hp2,hp1),s.hp1_prev))
    # Actual compiler rounds the gain product, then MADDs decay*old.
    v=B.fma(d,s.velocity,B.mul(gain,hp2))
    e=B.fma(d,s.elevation,B.mul(gain,v))
    keywords=dict(decay_exp=d,velocity_successor=v,elevation_successor=e)
    horizon_w=raw_binding=log_w=usable_w=None
    canonical=None
    if B.add(s.elapsed,h)>=B.div(3,CFG.lambda_):
        if before.log.log_period is None:
            canonical=F(6)
        else:
            canonical=call('exp',before.log.log_period)
            horizon_w=X.HorizonPeriodWitness(before.log.log_period,canonical)
        horizon=min(CFG.max_horizon_sec,max(CFG.min_horizon_sec,
                    B.mul(CFG.moment_horizon_periods,canonical)))
        md=call('exp',-B.div(h,horizon))
        a=B.sub(1,md); om=B.sub(1,a)
        av,ae=B.mul(a,v),B.mul(a,e)
        moments=M.MomentSuccessors(
            B.fma(om,s.weight,a),
            B.fma(om,s.velocity_mean,av),
            B.fma(om,s.velocity_sq,B.mul(av,v)),
            B.fma(om,s.elevation_mean,ae),
            B.fma(om,s.elevation_sq,B.mul(ae,e)))
        keywords.update(moment_decay_exp=md,moment_successors=moments)
        if moments.weight>M.WEIGHT_GATE:
            vm=B.div(moments.velocity_mean,moments.weight)
            em=B.div(moments.elevation_mean,moments.weight)
            vv=max(F(0),B.fma(-vm,vm,B.div(moments.velocity_sq,moments.weight)))
            ev=max(F(0),B.fma(-em,em,B.div(moments.elevation_sq,moments.weight)))
            keywords.update(velocity_var_successor=vv,elevation_var_successor=ev)
            if vv>M.VAR_GATE and ev>M.VAR_GATE:
                omega=B.fma(-CFG.lambda_,CFG.lambda_,B.div(vv,ev))
                keywords['omega_sq_successor']=omega
                if omega>M.OMEGA_GATE:
                    keywords['sqrt_omega']=call('sqrt',omega)
    # Existing relation checks every selected successor against its same-state
    # arithmetic alternatives, all early exits, and the raw-period ancestry.
    moment=M.step(s,CFG,dt=h,vertical_accel=x,
                  canonical_period=None if before.log.log_period is None else canonical,
                  libm_profile=U.ERROR_PROFILE,**keywords)
    log_after=before.log; log_step=None
    produced=moment.branch=='valid-period'
    if produced:
        raw=moment.raw_period
        log_raw=call('log',raw)
        raw_binding=X.RawLogBinding(raw,log_raw)
        if before.log.log_period is None:
            log_w=L.InitWitness(log_raw)
            log_step=L._init(before.log,log_w,True)
        else:
            sea=call('exp',before.log.log_period)
            horizon=min(L.HORIZON_MAX,max(L.HORIZON_MIN,B.mul(L.LOG_SMOOTH_PERIODS,sea)))
            decay=call('exp',-B.div(h,horizon))
            log_w=L.SmoothWitness(log_raw,sea,decay)
            log_step=L._smooth(before.log,log_w,True)
        log_after=log_step.after
        if not before.usable:
            usable_w=US.PeriodWitness(log_after.log_period,call('exp',log_after.log_period))
    usable=US.update(before.usable,produced_period=produced,
        log_period=log_after.log_period,elapsed=moment.state.elapsed,
        lambda_=CFG.lambda_,witness=usable_w)
    witness=X.ModeWitnesses(keywords,horizon_w,raw_binding,log_w,usable_w)
    return Step(before,State(moment.state,log_after,usable.after),moment,
                log_step,usable,witness,h,x)


def attach_to_product(before:X.State,selected:Step,*,separate:X.ModeWitnesses):
    """Lift one target step into the SAME persistent product .fma history."""
    if not isinstance(before,X.State) or not isinstance(selected,Step):
        raise TypeError('machine product and selected target transition required')
    projection=State(before.fma,before.logs.fma,before.fma_usable)
    if projection!=selected.before or before.cfg!=CFG:
        raise ValueError('target transition detached from persistent FMA pre-state')
    if before.libm_profile!=U.ERROR_PROFILE:
        raise ValueError('target transition requires qualified tolerant libm profile')
    after,*details=X.step(before,dt=selected.dt,vertical_accel=selected.vertical_accel,
                         separate=separate,fma=selected.witnesses)
    if State(after.fma,after.logs.fma,after.fma_usable)!=selected.state:
        raise AssertionError('target graph not included in persistent FMA history')
    return after,details


def _run(args):
    p=subprocess.run([str(x) for x in args],text=True,capture_output=True)
    if p.returncode: raise RuntimeError('target WPE audit command failed: '+p.stderr)
    return p.stdout


def arithmetic_sites(assembly):
    end=assembly.index('\t.size\t_ZN19WavePeriodEstimator6updateEff')
    assembly=assembly[:end]
    filematch=re.search(r'\.file\s+(\d+)\s+"[^"\n]*tuner/WavePeriodEstimator.h"',assembly)
    if not filematch: raise ValueError('WPE source debug map absent')
    wpe_id=int(filematch[1]); location=(0,0); sites=defaultdict(list)
    all_fused=[]
    for line in assembly.splitlines():
        match=re.match(r'\s*\.loc\s+(\d+)\s+(\d+)',line)
        if match: location=tuple(map(int,match.groups()))
        match=re.match(r'\s*(add\.s|sub\.s|mul\.s|madd\.s|msub\.s)\s+([^#]+)',line)
        if match:
            if match[1] in ('madd.s','msub.s'): all_fused.append(location)
            if location[0]==wpe_id:
                sites[location[1]].append((match[1],match[2].strip()))
    for source_line,expected in EXPECTED.items():
        if tuple(x[0] for x in sites[source_line])!=expected:
            raise ValueError('WPE target operation selection changed at source line '+str(source_line))
    expected_fused=[(wpe_id,line) for line,ops in EXPECTED.items()
                    for op in ops if op in ('madd.s','msub.s')]
    if sorted(all_fused)!=sorted(expected_fused):
        raise ValueError('unmodeled target WPE contraction site')
    return {str(k):v for k,v in sites.items()}


def audit(toolchain:Path,sdk:Path):
    source=ROOT/'src/tuner/WavePeriodEstimator.h'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=SOURCE_SHA:
        raise ValueError('shipping WPE source changed; re-audit target expression map')
    limits=ROOT/'src/tuner/SeaStateAdaptationLimits.h'
    if hashlib.sha256(limits.read_bytes()).hexdigest()!=LIMITS_SHA:
        raise ValueError('WPE dynamic horizon limits changed; re-audit profile')
    compiler=toolchain/'bin/xtensa-esp32s3-elf-g++'
    flags=shlex.split((sdk/'flags/cpp_flags').read_text())+['-Os','-funroll-loops','-fno-finite-math-only']
    with tempfile.TemporaryDirectory(prefix='ou3-wpe-compiler-') as directory:
        tmp=Path(directory); cpp=tmp/'probe.cpp'; assembly=tmp/'probe.s'; obj=tmp/'probe.o'
        cpp.write_text(PROBE)
        common=[compiler,*flags,'-ffile-prefix-map='+str(tmp)+'=/ou3-wpe-compiler',
                '-ffile-prefix-map='+str(ROOT)+'=/ocean-imu','-I'+str(ROOT/'src')]
        _run([*common,'-fverbose-asm','-S',cpp,'-o',assembly])
        _run([*common,'-c',cpp,'-o',obj])
        assembly_text=assembly.read_text()
        sites=arithmetic_sites(assembly_text)
        body=assembly_text.split('\t.size\t_ZN19WavePeriodEstimator6updateEff')[0]
        call_counts={name:len(re.findall(r'\bcall8\s+'+re.escape(name)+r'\s',body))
                     for name in ('expf','logf','sqrtf','__divsf3')}
        if call_counts!={'expf':6,'logf':1,'sqrtf':1,'__divsf3':11}:
            raise ValueError('WPE math-call namespace or branch sites changed')
        obj_sha=hashlib.sha256(obj.read_bytes()).hexdigest()
    profile=LIBM.profile_correspondence()
    if not profile['pinned_scalar_profile_attached']:
        raise ValueError('WPE MADD/MSUB instruction semantics not qualified')
    return {'qualification':QUALIFICATION,'shipping_WPE_source_sha256':SOURCE_SHA,
        'dynamic_horizon_limits_sha256':LIMITS_SHA,'math_call_sites':call_counts,
        'compiler_driver_sha256':hashlib.sha256(compiler.resolve().read_bytes()).hexdigest(),
        'sdk_cpp_flags_sha256':hashlib.sha256((sdk/'flags/cpp_flags').read_bytes()).hexdigest(),
        'compiler_flags':flags,
        'probe_object_sha256':obj_sha,'source_arithmetic_instruction_sites':sites,
        'selected_persistent_log_track':'fma',
        'integrator_and_moment_choice':'FMA(left factors, rounded right product)',
        'variance_and_omega_choice':'fused MSUB on same computed operands',
        'log_EMA_choice':'FMA(alpha, RN(log_raw-log_previous), log_previous)',
        'all_target_contraction_sites_mapped_to_existing_machine_relation':True,
        'target_step_constructor_selects_same_history_operands':True,
        'persistent_FMA_projection_preserved_by_transition_inclusion':True,
        'pinned_WPE_compiler_profile_selection_closed':True,
        'whole_firmware_compiler_and_link_correspondence_closed':False,
        'source_uniform_WPE_input_supply_qualified_here':False,
        'rounding_premise':profile['rounding_premise']}


def _check_pinned_compiler(report):
    path=Path(__file__).resolve().parents[1]/'ou3_alt_target_arithmetic_audit.json'
    target=json.loads(path.read_text())
    if (target.get('qualification')!='OU3_ALT_PINNED_TARGET_ARITHMETIC_AUDIT_V1'
            or any(report.get(key)!=target.get(key) for key in
                   ('compiler_driver_sha256','sdk_cpp_flags_sha256','compiler_flags'))):
        raise ValueError('WPE compiler or flags differ from pinned target build')


def readiness(*,report=None):
    path=Path(__file__).resolve().parents[1]/'ou3_alt_target_wpe_compiler.json'
    if report is None:
        report=json.loads(path.read_text())
    _check_pinned_compiler(report)
    source_sha=hashlib.sha256((ROOT/'src/tuner/WavePeriodEstimator.h').read_bytes()).hexdigest()
    limits_sha=hashlib.sha256((ROOT/'src/tuner/SeaStateAdaptationLimits.h').read_bytes()).hexdigest()
    if (report.get('qualification')!=QUALIFICATION
            or report.get('shipping_WPE_source_sha256')!=source_sha
            or report.get('dynamic_horizon_limits_sha256')!=limits_sha):
        raise ValueError('stale pinned WPE compiler qualification')
    profile=LIBM.profile_correspondence()
    return {'pinned_WPE_compiler_profile_selection_closed':(
        report.get('all_target_contraction_sites_mapped_to_existing_machine_relation') is True
        and profile['pinned_scalar_profile_attached']),
        'persistent_actual_target_track':'fma',
        'whole_firmware_compiler_and_link_correspondence_closed':False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--toolchain',type=Path,required=True)
    p.add_argument('--sdk',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(); result=audit(args.toolchain.resolve(),args.sdk.resolve())
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('selected_persistent_log_track',
        'pinned_WPE_compiler_profile_selection_closed','whole_firmware_compiler_and_link_correspondence_closed')}))


if __name__=='__main__': main()
