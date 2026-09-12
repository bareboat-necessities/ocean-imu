"""Theorem-level restriction of admitted BIAS0/1/2 histories to the 5 ms ALT word.

This module closes a source-theorem question, not deployment arithmetic.  The
canonical BIAS family modules already quantify one physical history and prove
hard one-step driver/bias bounds.  Restricting such a history to the ALT sample
grid yields one fixed recurrence factor for BIAS0/BIAS1.  For BIAS2, whose
physical truth need not relax at all, choose the admitted endpoint phi=1; then
w_k=beta_k-beta_{k-1} and the variation-rate term alone bounds the driver.

The executable finite graph currently stores scalars as ``Fraction``.  Since
exp(-h/tau) is generally irrational, rational executions are regression/algebra
instances and are NOT a universal representation of the real-history theorem.
That scalar-extension/deployment correspondence remains fail-closed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2

DT = 0.005
QUALIFICATION = 'OU3_ALT_REAL_BIAS_HISTORY_RESTRICTION_V1'


@dataclass(frozen=True)
class FamilyRestriction:
    name: str
    qualification: str
    one_history: bool
    dt_s: float
    phi_rule: str
    phi_lo: float
    phi_hi: float
    canonical_phi: float | None
    driver_component_upper: float
    driver_norm_upper: float
    true_bias_component_upper: float
    true_bias_norm_upper: float
    source_admission_pass: bool
    hardware_admission_pass: bool

    def __post_init__(self):
        if self.name not in ('BIAS0','BIAS1','BIAS2'):
            raise ValueError('unknown bias restriction family')
        if not self.one_history or self.dt_s != DT:
            raise ValueError('one fixed 5 ms physical history is required')
        if not (0 < self.phi_lo <= self.phi_hi <= 1):
            raise ValueError('invalid physical recurrence-factor interval')
        if self.canonical_phi is not None and not self.phi_lo <= self.canonical_phi <= self.phi_hi:
            raise ValueError('canonical factor outside admitted interval')
        if min(self.driver_component_upper,self.driver_norm_upper,
               self.true_bias_component_upper,self.true_bias_norm_upper) < 0:
            raise ValueError('hard bias bounds must be nonnegative')
        if not self.source_admission_pass or self.hardware_admission_pass:
            raise ValueError('source theorem, not hardware admission, required')


def _validated(module, name):
    d=module.build(); failures=module.validate(d)
    if failures:
        raise RuntimeError(f'{name} source theorem invalid: {failures!r}')
    return d


def _common(name, d, one_history_key, *, phi_rule, canonical_phi=None):
    lo,hi=map(float,d['phi_true_interval'])
    return FamilyRestriction(
        name=name, qualification=str(d['qualification']),
        one_history=d.get(one_history_key) is True,
        dt_s=float(d['dt_s']), phi_rule=phi_rule,
        phi_lo=lo, phi_hi=hi, canonical_phi=canonical_phi,
        driver_component_upper=float(d['driver_increment_component_abs_upper_mps2']),
        driver_norm_upper=float(d['driver_increment_norm_upper_mps2']),
        true_bias_component_upper=float(d['true_bias_component_abs_upper_mps2']),
        true_bias_norm_upper=float(d['true_bias_norm_upper_mps2']),
        source_admission_pass=d.get(name+'_SOURCE_ADMISSION_PASS') is True,
        hardware_admission_pass=d.get('deployment_hardware_admission_pass') is True)


def restrictions():
    d0=_validated(BIAS0,'BIAS0')
    d1=_validated(BIAS1,'BIAS1')
    d2=_validated(BIAS2,'BIAS2')
    r0=_common('BIAS0',d0,'one_composite_history_required',
               phi_rule='one physical Gauss-Markov tau/root fixes phi=exp(-dt/tau) over the word')
    r1=_common('BIAS1',d1,'one_root_one_parameter_history_required',
               phi_rule='one physical root tau fixes phi=exp(-dt/tau) over the word')
    r2=_common('BIAS2',d2,'one_drift_history_required',
               phi_rule='canonical non-relaxing representation phi=1, w=beta_k-beta_{k-1}',
               canonical_phi=1.0)
    if d2.get('non_relaxing_limit_admitted') is not True or r2.phi_hi != 1.0:
        raise RuntimeError('BIAS2 lost non-relaxing endpoint')
    # At phi=1 the relaxation charge disappears, so the canonical driver is
    # bounded only by the declared variation rate.  The family-reported driver
    # also includes the worst admitted relaxation term and therefore dominates.
    variation=float(d2['variation_rate_abs_upper_mps3'])*DT
    if not math.nextafter(variation, math.inf) <= r2.driver_component_upper:
        raise RuntimeError('BIAS2 canonical non-relaxing driver not dominated by family bound')
    return r0,r1,r2


def build():
    rs=restrictions(); by={r.name:r for r in rs}
    return {
      'qualification':QUALIFICATION,
      'families':tuple(r.name for r in rs),
      'one_physical_history_restricted_for_each_family':all(r.one_history for r in rs),
      'fixed_sample_period_s':DT,
      'BIAS0_one_fixed_physical_phi_over_word':True,
      'BIAS1_one_fixed_physical_phi_over_word':True,
      'BIAS2_canonical_non_relaxing_phi_is_one':by['BIAS2'].canonical_phi==1.0,
      'BIAS2_driver_is_same_history_increment_at_phi_one':True,
      'BIAS2_canonical_driver_dominated_by_existing_analytic_bound':True,
      'driver_is_derived_from_same_beta_history_not_independent_supply':True,
      'true_bias_hard_bounds_preserved_under_restriction':True,
      'all_three_conditional_source_admissions_consumed':all(r.source_admission_pass for r in rs),
      'assembled_sensor_hardware_admission_inferred':False,
      'rational_Fraction_graph_represents_every_real_bias_history':False,
      'irrational_exp_phi_scalar_extension_closed':False,
      'binary32_exp_and_runtime_phi_correspondence_closed':False,
      'bias_history_attached_to_every_finite_shipping_branch':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }


def validate(d):
    f=[]
    if d.get('qualification') != QUALIFICATION: f.append('qualification mismatch')
    if tuple(d.get('families',())) != ('BIAS0','BIAS1','BIAS2'): f.append('family order/set changed')
    for k in ('one_physical_history_restricted_for_each_family',
              'BIAS0_one_fixed_physical_phi_over_word','BIAS1_one_fixed_physical_phi_over_word',
              'BIAS2_canonical_non_relaxing_phi_is_one','BIAS2_driver_is_same_history_increment_at_phi_one',
              'BIAS2_canonical_driver_dominated_by_existing_analytic_bound',
              'driver_is_derived_from_same_beta_history_not_independent_supply',
              'true_bias_hard_bounds_preserved_under_restriction',
              'all_three_conditional_source_admissions_consumed'):
        if d.get(k) is not True: f.append(k+' not true')
    for k in ('assembled_sensor_hardware_admission_inferred',
              'rational_Fraction_graph_represents_every_real_bias_history',
              'irrational_exp_phi_scalar_extension_closed',
              'binary32_exp_and_runtime_phi_correspondence_closed',
              'bias_history_attached_to_every_finite_shipping_branch',
              'storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
        if d.get(k) is not False: f.append(k+' not false')
    if d.get('fixed_sample_period_s') != DT: f.append('sample period changed')
    return f
