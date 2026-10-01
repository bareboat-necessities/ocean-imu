"""Exact source-uniform span envelope for packet-indistinguishable attitude/BA gauge."""
import json, math
from pathlib import Path


def gauge_span_envelope(window_s):
    c=json.loads(Path(__file__).with_name("constants.json").read_text())
    g=9.80665
    b=c["imu_bias"]; m=c["marine_motion"]
    A=2.0*math.asin(b["B_a_mps2"]/(2.0*g))
    L=min(b["D_a_mps3"]/g,b["B_g_rad_s"],m["Omega_max_rad_s"])
    return min(2.0*A,L*window_s)


def certificate():
    c=json.loads(Path(__file__).with_name("constants.json").read_text())
    g=9.80665; b=c["imu_bias"]; m=c["marine_motion"]
    A=2.0*math.asin(b["B_a_mps2"]/(2.0*g))
    L=min(b["D_a_mps3"]/g,b["B_g_rad_s"],m["Omega_max_rad_s"])
    return {"A_g_rad":A,"peak_to_peak_cap_rad":2*A,"L_g_rad_s":L,
            "T_sat_s":2*A/L,"active_speed_bound":"D_a/g",
            "D_g_improves_arbitrary_interior_window_range":False,
            "formula":"min(4 asin(B_a/(2g)), T min(D_a/g,B_g,Omega_max))"}

if __name__=="__main__": print(json.dumps(certificate(),indent=2,sort_keys=True))
