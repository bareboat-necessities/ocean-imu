"""First literal Live IMU step from the exact startup-produced H18 state.

This is a substitution/composition layer, not a new filter approximation.  It
accepts only ``finite_startup_live_runtime_bridge.Result`` and requires the raw
IMU packet plus physical segment to begin at that exact fresh H18 physical
reference.  It then calls the already-proved ``finite_live_imu_prefix.step``
without modifying any event algebra.

Therefore the represented startup path now reaches one actual Live prediction /
S-service / accelerometer / tuner-WPE prefix with no synthetic Live root between
them.  Remaining open branches (tilt reset, async mag at this exact prefix,
finite precision, source-uniform continuation) stay fail-closed.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as BRIDGE
from tools.stability.ou3_alt_contraction import finite_live_imu_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS

@dataclass(frozen=True)
class Result:
    startup: BRIDGE.Result
    first_live: LIVE.Result


def step(startup:BRIDGE.Result,raw:SENSOR.RawImuSample,segment:PHYS.PhysicalSegment,**kwargs):
    if not isinstance(startup,BRIDGE.Result):
        raise TypeError('startup-to-Live runtime bridge result required')
    if not isinstance(raw,SENSOR.RawImuSample) or not isinstance(segment,PHYS.PhysicalSegment):
        raise TypeError('raw IMU sample and physical segment required')
    root=startup.state.mekf.reference
    if raw.physical != root or segment.before != root:
        raise ValueError('first Live packet/segment detached from exact fresh startup reference')
    out=LIVE.step(startup.state,raw,segment,**kwargs)
    return Result(startup,out)


def readiness():
    return {
      'exact_startup_bridge_consumed_as_first_Live_state':True,
      'first_Live_raw_packet_rooted_at_fresh_startup_reference':True,
      'first_Live_physical_segment_rooted_at_fresh_startup_reference':True,
      'first_prediction_S_accel_tuner_WPE_prefix_composed_from_startup':True,
      'synthetic_Live_root_removed_at_first_sample':True,
      'async_mag_at_first_prefix_composed':False,
      'firing_tilt_reset_at_first_prefix_composed':False,
      'deployment_finite_precision_closed':False,
      'source_uniform_indefinite_continuation_closed':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
