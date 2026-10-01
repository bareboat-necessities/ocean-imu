"""Analytical gauge bounds from the carried SLOW + FAST IMU contract."""
import json, math
from pathlib import Path

def slow_quiet_alias_envelope(window_s, gravity, accel_slow_amplitude, accel_slow_rate, gyro_slow_amplitude):
    vals=(window_s,gravity,accel_slow_amplitude,accel_slow_rate,gyro_slow_amplitude)
    if not all(math.isfinite(x) and x>0 for x in vals): raise ValueError("positive finite inputs required")
    chord=min(2*accel_slow_amplitude,accel_slow_rate*window_s)
    accel=2*math.asin(min(1.0,chord/(2*gravity)))
    gyro=min(math.pi,gyro_slow_amplitude*window_s)
    return {"accel_span_ceiling_rad":accel,"gyro_span_ceiling_rad":gyro,
            "joint_span_ceiling_rad":min(accel,gyro),"general_waveform":True,
            "fast_components_assumed_zero":True}

def qualified_fast_endpoint_envelope(window_s, gravity, accel_slow_amplitude, accel_slow_rate,
                                     accel_fast_amplitude, fast_cell_cap_0=None, fast_cell_cap_1=None):
    vals=(window_s,gravity,accel_slow_amplitude,accel_slow_rate,accel_fast_amplitude)
    if not all(math.isfinite(x) and x>0 for x in vals): raise ValueError("positive finite inputs required")
    slow=min(2*accel_slow_amplitude,accel_slow_rate*window_s)
    if fast_cell_cap_0 is None or fast_cell_cap_1 is None:
        return {"span_ceiling_rad":None,"temporal_qualification":"OPEN","full_joint_gauge_exclusion":False}
    caps=(fast_cell_cap_0,fast_cell_cap_1)
    if not all(math.isfinite(x) and 0<=x<=accel_fast_amplitude for x in caps):
        raise ValueError("fast endpoint caps outside amplitude ball")
    chord=slow+sum(caps)
    theta=math.pi if chord>=2*gravity else 2*math.asin(chord/(2*gravity))
    return {"span_ceiling_rad":theta,"reachable_accel_chord_upper_mps2":chord,
            "temporal_qualification":"CONDITIONAL","full_joint_gauge_exclusion":False}

def certificate():
    c=json.loads(Path(__file__).with_name("constants.json").read_text()); b=c["imu_bias"]; m=c["marine_motion"]
    examples={}
    for T in (17.,60.,100.):
        r=slow_quiet_alias_envelope(T,9.80665,b["B_a_s_mps2"],b["D_a_s_mps3"],b["B_g_s_rad_s"])
        examples[str(int(T))+"s"]={"joint_span_ceiling_rad":r["joint_span_ceiling_rad"],
                                  "joint_span_ceiling_deg":math.degrees(r["joint_span_ceiling_rad"])}
    return {"qualification":"OU3_TWO_TIMESCALE_GAUGE_BOUND_V1",
            "universal_slow_only_quiet_alias_bound_proved":True,
            "conditional_fast_endpoint_extension_proved":True,
            "marine_T_E_numeric":m["attitude_excitation"]["T_E_s"],
            "marine_theta_E_numeric":m["attitude_excitation"]["theta_E_rad"],
            "current_marine_slow_only_gauge_exclusion":False,
            "full_slow_fast_joint_gauge_exclusion":False,
            "stronger_marine_assumption_added":False,"examples":examples,"theorem_closed":False}

if __name__=="__main__": print(json.dumps(certificate(),indent=2,sort_keys=True))
