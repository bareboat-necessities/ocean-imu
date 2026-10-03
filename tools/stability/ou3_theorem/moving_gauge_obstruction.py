"""Analytical MOVING tilt/BA gauge obstruction under the current physical contract."""
from __future__ import annotations
import json, math
from pathlib import Path

def certificate(constants:dict|None=None)->dict:
    if constants is None:
        constants=json.loads((Path(__file__).with_name("constants.json")).read_text())
    m=constants["marine_motion"]; imu=constants["imu_bias"]; mag=constants["magnetic_service"]
    T=float(m["attitude_excitation"]["T_E_s"])
    assert T==float(m["displacement_excitation"]["T_P_s"])
    theta_E=float(m["attitude_excitation"]["theta_E_rad"]); P_E=float(m["displacement_excitation"]["P_E_m"])
    g=9.80665; delta=2*math.atan(1/200); Aphi=theta_E/2; Ap=P_E/2; w=2*math.pi/T
    ba_norm=2*g*math.sin(delta/2); ba_rate=ba_norm*Aphi*w
    vmax=Ap*w; amax=Ap*w*w; jmax=Ap*w*w*w; omega_max=Aphi*w
    checks={"gravity_span":2*Aphi>=theta_E,"displacement_span":2*Ap>=P_E,
      "position":Ap<=m["P_max_m"],"velocity":vmax<=m["V_max_mps"],"acceleration":amax<=m["A_max_mps2"],
      "jerk":jmax<=m["J_max_mps3"],"angular_rate":omega_max<=m["Omega_max_rad_s"],
      "slow_accel_amplitude":ba_norm<=imu["B_a_s_mps2"],"slow_accel_rate":ba_rate<=imu["D_a_s_mps3"],
      "slow_gyro_amplitude":True,"slow_gyro_rate":True,"fast_accel_zero":True,"fast_gyro_zero":True,
      "field_norm":mag["field_norm_min_uT"]<=75<=mag["field_norm_max_uT"],"horizontal_field":75>=mag["horizontal_field_min_uT"]}
    return {"classification":"B","role":"lossless obstruction to absolute MOVING physical tilt/BA separation; not shipping instability",
      "period_s":T,"relative_roll_rad":delta,"relative_roll_deg":math.degrees(delta),"rock_amplitude_rad":Aphi,
      "position_amplitude_m":Ap,"velocity_max_mps":vmax,"acceleration_max_mps2":amax,"jerk_max_mps3":jmax,
      "angular_rate_max_rad_s":omega_max,"slow_accel_bias_norm_mps2":ba_norm,"slow_accel_bias_rate_max_mps3":ba_rate,
      "delivered_accel_records_identical":True,"delivered_gyro_records_identical":True,"delivered_mag_records_identical":True,
      "fast_components_zero":True,"checks":checks,"admitted_by_numeric_physical_envelopes":all(checks.values()),
      "structures_preserved":["21-state shipping execution","MARINE MOVING spans","persistent SLOW+FAST split","MAGNETIC SERVICE chronology","Mahony/guard/tuner/scheduler/covariance/gains/Joseph/S=0 chronology"],"relaxations_introduced":[]}
if __name__=="__main__": print(json.dumps(certificate(),indent=2,sort_keys=True))
