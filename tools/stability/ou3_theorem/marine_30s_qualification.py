"""Exact consequences of the user-qualified 30-s MARINE excitation profile."""
import math
T_E=30.0
THETA_E=math.radians(2.0)
T_P=30.0
P_E=0.03
V_MAX=5.5

def certificate():
    return {
      "qualification":"OU3_MARINE_30S_USER_QUALIFIED_V1",
      "T_E_s":T_E,"theta_E_rad":THETA_E,"theta_E_deg":2.0,
      "T_P_s":T_P,"P_E_m":P_E,
      "attitude_span_average_scale_rad_s":THETA_E/T_E,
      "displacement_velocity_chord_lower_mps":P_E/T_P,
      "signed_acceleration_chord_lower_m":
          max(0.0,P_E-T_P*V_MAX),
      "pointwise_acceleration_floor_inferred":False,
      "zero_translation_MOVING_excluded":True,
      "fast_IMU_temporal_qualification_closed":False,
      "full_joint_gauge_breaking_closed":False,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
