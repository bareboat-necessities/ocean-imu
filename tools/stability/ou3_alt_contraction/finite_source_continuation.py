"""Persistent COMPLETE-BRMM/BIAS/sensor ancestry for finite Live segments.

``finite_physical_prediction.PhysicalSegment`` proves an exact physical
recurrence but deliberately does *not* certify source admission. This module
adds the theorem-facing ancestry layer without converting the correlated
601-sample source relation into independent per-sample boxes.

A qualified continuation carries one O^601_BRMM source-history/generator root,
one Live S origin, one BIAS0/BIAS1/2 contract/token, consecutive source-cell and
physical-primitive identities, and the exact finite ``PhysicalSegment``. A
qualified raw IMU packet additionally carries persistent gyro- and
accelerometer-residual history tokens and is forced to use that segment's exact
predecessor physical endpoint and the ALT identity de-heel map. A qualified
physical endpoint is obtained only as the ``before`` or ``after`` endpoint of
one such admitted transition; asynchronous events can therefore attach to an
actual member of the correlated source continuation rather than a matching time
or history string alone.

The residual values remain explicit ISS/finite-horizon forcing. This module does
NOT reinterpret configured Racc/Rmag covariance standard deviations as hard
pathwise sensor-noise bounds. Quantitative sensor residual envelopes and target
conversion roundoff remain open. The remaining decisive bridge is estimator
ownership: every tuner/guard/measurement/prediction coefficient product in the
finite shipping event still has to be generated from this same qualified
continuation. Therefore this module cannot enable storage.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

import ou3_brmm_correlated_window_outer_enclosure as OUTER
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR

OUTER_RELATION = 'O^601_BRMM'
TRANSITIONS = 600
SAMPLES = 601
DT = F(1,200)
QUALIFICATION = 'OU3_ALT_FINITE_SOURCE_CONTINUATION_V1'


def R(x): return M.rational(x)


def _contract(name: str) -> BIAS.BiasFamilyContract:
    cs={c.name:c for c in BIAS.contracts()}
    if name not in cs: raise ValueError('unknown ALT bias family '+repr(name))
    return cs[name]


def _norm2(v):
    q=tuple(M.vec(v,3)); return sum(x*x for x in q)


def _le_float_bound_square(v, bound: float) -> bool:
    b=R(str(float(bound)))
    return _norm2(v) <= b*b


@dataclass(frozen=True)
class SourceRoot:
    history_id: str
    generator_id: str
    live_origin: F
    bias_family: str
    bias_parameter_token: str
    outer_relation: str = OUTER_RELATION
    qualification: str = QUALIFICATION
    def __post_init__(self):
        if any(not isinstance(x,str) or not x for x in (self.history_id,self.generator_id)):
            raise ValueError('persistent physical history/generator ids required')
        object.__setattr__(self,'live_origin',R(self.live_origin))
        c=_contract(self.bias_family)
        if self.bias_parameter_token != c.parameter_token:
            raise ValueError('bias parameter token detached from selected physical family')
        if self.outer_relation != OUTER_RELATION or self.qualification != QUALIFICATION:
            raise ValueError('wrong COMPLETE-BRMM finite-source qualification')


@dataclass(frozen=True)
class SensorDisturbanceRoot:
    source_root: SourceRoot
    gyro_residual_history_id: str
    accel_residual_history_id: str
    conversion_profile: str='ALT_ZERO_HEEL_RAW_BODY_V1'
    def __post_init__(self):
        if not isinstance(self.source_root,SourceRoot): raise TypeError('physical source root required')
        for x in (self.gyro_residual_history_id,self.accel_residual_history_id,self.conversion_profile):
            if not isinstance(x,str) or not x: raise ValueError('persistent sensor disturbance/profile id required')
        if self.conversion_profile!='ALT_ZERO_HEEL_RAW_BODY_V1':
            raise ValueError('unsupported theorem sensor conversion profile')


@dataclass(frozen=True)
class StepWitness:
    ordinal: int
    parent_source_cell_id: str
    source_cell_id: str
    primitive_in_id: str
    primitive_out_id: str
    def __post_init__(self):
        if not isinstance(self.ordinal,int) or isinstance(self.ordinal,bool) or not 1 <= self.ordinal <= TRANSITIONS:
            raise ValueError('source ordinal must be one of the 600 transitions')
        for x in (self.parent_source_cell_id,self.source_cell_id,self.primitive_in_id,self.primitive_out_id):
            if not isinstance(x,str) or not x: raise ValueError('persistent source/primitive ids required')


@dataclass(frozen=True)
class QualifiedPhysicalSegment:
    root: SourceRoot
    witness: StepWitness
    segment: PHYS.PhysicalSegment
    def __post_init__(self):
        if not isinstance(self.root,SourceRoot) or not isinstance(self.witness,StepWitness) or not isinstance(self.segment,PHYS.PhysicalSegment):
            raise TypeError('source root, step witness and exact PhysicalSegment required')
        s=self.segment; c=_contract(self.root.bias_family)
        if s.h != DT:
            raise ValueError('canonical O^601_BRMM continuation requires 5 ms transition')
        if s.before.live_origin != self.root.live_origin or s.after.live_origin != self.root.live_origin:
            raise ValueError('physical segment detached from one-time Live S origin')
        expected_before=self.root.live_origin + (self.witness.ordinal-1)*DT
        if s.before.time != expected_before or s.after.time != expected_before+DT:
            raise ValueError('segment clock detached from 601-sample source ordinal')
        phi=R(s.phi_true)
        lo=R(str(c.phi_true.lo)); hi=R(str(c.phi_true.hi))
        if not lo <= phi <= hi:
            raise ValueError('physical beta decay detached from selected BIAS family')
        cb=R(str(c.driver_component_bound))
        if any(abs(x)>cb for x in s.bias_driver):
            raise ValueError('bias driver exceeds family component contract')
        if not _le_float_bound_square(s.bias_driver,c.driver_norm_bound):
            raise ValueError('bias driver exceeds family norm contract')
        for beta in (s.before.beta,s.after.beta):
            if not _le_float_bound_square(beta,c.true_bias_norm_bound):
                raise ValueError('true accelerometer bias exceeds family hard norm contract')


@dataclass(frozen=True)
class QualifiedPhysicalEndpoint:
    """An endpoint that is proven to belong to one admitted source transition."""
    physical: QualifiedPhysicalSegment
    side: str
    def __post_init__(self):
        if not isinstance(self.physical,QualifiedPhysicalSegment):
            raise TypeError('qualified physical segment required')
        if self.side not in ('before','after'):
            raise ValueError("endpoint side must be 'before' or 'after'")
    @property
    def root(self): return self.physical.root
    @property
    def endpoint(self): return getattr(self.physical.segment,self.side)


@dataclass(frozen=True)
class QualifiedRawImuSample:
    physical: QualifiedPhysicalSegment
    sensor_root: SensorDisturbanceRoot
    raw: SENSOR.RawImuSample
    packet_id: str
    def __post_init__(self):
        if not isinstance(self.physical,QualifiedPhysicalSegment) or not isinstance(self.sensor_root,SensorDisturbanceRoot) or not isinstance(self.raw,SENSOR.RawImuSample):
            raise TypeError('qualified physical step, sensor root and raw IMU packet required')
        if self.sensor_root.source_root != self.physical.root:
            raise ValueError('sensor residual histories detached from COMPLETE-BRMM/BIAS root')
        if self.raw.physical != self.physical.segment.before:
            raise ValueError('raw IMU packet detached from qualified physical predecessor')
        if self.raw.deheel_body_to_internal != SENSOR.IDENTITY3:
            raise ValueError('ALT source-qualified packet requires zero-heel identity map')
        if not isinstance(self.packet_id,str) or not self.packet_id:
            raise ValueError('persistent raw packet identity required')


@dataclass(frozen=True)
class Continuation:
    root: SourceRoot
    steps: tuple[QualifiedPhysicalSegment,...]
    def __post_init__(self):
        if not isinstance(self.root,SourceRoot): raise TypeError('source root required')
        if len(self.steps)>TRANSITIONS: raise ValueError('more than 600 transitions')
        prev=None
        for k,q in enumerate(self.steps,1):
            if not isinstance(q,QualifiedPhysicalSegment) or q.root!=self.root:
                raise ValueError('source continuation restarted its root')
            if q.witness.ordinal!=k:
                raise ValueError('source continuation skipped/duplicated an ordinal')
            if k==1:
                if q.witness.parent_source_cell_id!='root':
                    raise ValueError('first source-cell parent must be root')
            else:
                if q.witness.parent_source_cell_id!=prev.witness.source_cell_id:
                    raise ValueError('source-cell parent/child ancestry broken')
                if q.witness.primitive_in_id!=prev.witness.primitive_out_id:
                    raise ValueError('physical primitive continuity broken')
                if q.segment.before!=prev.segment.after:
                    raise ValueError('finite physical endpoint continuity broken')
            prev=q
    @property
    def complete(self): return len(self.steps)==TRANSITIONS
    @property
    def next_ordinal(self): return len(self.steps)+1


def certified_root(*,history_id,generator_id,live_origin,bias_family):
    outer=OUTER.build(); failures=OUTER.validate(outer)
    if failures: raise RuntimeError('correlated COMPLETE-BRMM outer relation invalid: '+repr(failures))
    if not (outer['left_inclusion_closed'] and outer['same_history_required_for_entire_window']
            and outer['correlation_retained_across_samples'] and outer['sample_count']==SAMPLES):
        raise RuntimeError('correlated COMPLETE-BRMM theorem lost required whole-history relation')
    c=_contract(bias_family)
    return SourceRoot(history_id,generator_id,R(live_origin),c.name,c.parameter_token)


def append(cont: Continuation, *, witness: StepWitness, segment: PHYS.PhysicalSegment):
    if not isinstance(cont,Continuation): raise TypeError('qualified continuation required')
    q=QualifiedPhysicalSegment(cont.root,witness,segment)
    return Continuation(cont.root,cont.steps+(q,))


def qualify_raw_imu(physical:QualifiedPhysicalSegment,sensor_root:SensorDisturbanceRoot,
                    raw:SENSOR.RawImuSample,packet_id:str):
    return QualifiedRawImuSample(physical,sensor_root,raw,packet_id)


def endpoint(physical:QualifiedPhysicalSegment,side:str):
    return QualifiedPhysicalEndpoint(physical,side)


def begin(root: SourceRoot):
    return Continuation(root,())


def readiness():
    outer=OUTER.build(); of=OUTER.validate(outer); bias=BIAS.build(); bf=BIAS.validate(bias)
    if of or bf: raise RuntimeError('finite source ancestry prerequisite failed')
    return {
      'qualification':QUALIFICATION,
      'correlated_COMPLETE_BRMM_left_inclusion_consumed':outer['left_inclusion_closed'],
      'one_outer_history_required_for_all_601_samples':outer['same_history_required_for_entire_window'],
      'independent_per_sample_BRMM_boxes_forbidden':True,
      'one_bias_family_parameter_token_over_word_required':True,
      'BIAS0_BIAS1_BIAS2_contracts_available':bias['three_families_invoked_separately'],
      'finite_segment_exact_physical_recurrence_consumed':True,
      'bias_phi_driver_and_true_beta_hard_contracts_checked_per_segment':True,
      'source_cell_parent_child_and_primitive_continuity_checked':True,
      'complete_600_transition_continuation_shape_materialized':True,
      'qualified_async_endpoint_comes_from_admitted_transition':True,
      'raw_IMU_packet_bound_to_same_qualified_physical_predecessor':True,
      'persistent_gyro_and_accel_residual_history_tokens_required':True,
      'Racc_covariance_not_reinterpreted_as_hard_sensor_noise_bound':True,
      'quantitative_sensor_residual_ISS_envelope_attached':False,
      'sensor_conversion_binary32_attached':False,
      'finite_estimator_coefficients_bound_to_same_source_continuation':False,
      'finite_magnetic_source_bound_to_same_COMPLETE_BRMM_history':False,
      'all_configured_hybrid_branches_bound_to_finite_source_graph':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
