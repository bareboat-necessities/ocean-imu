"""Candidate 30-s MARINE + 60-s FAST qualification consequences."""
from fractions import Fraction as F
from .imu_two_timescale_certificate import certificate as witness_certificate

def accumulation(T,B,H,C):
    n,rem=divmod(T,H)
    return min(B*T,n*C+min(B*rem,C))

def certificate():
    T=F(30); Ha=Hg=F(60); Ca=F(1,20); Cg=F(1,500)
    Baf=F(3,10); Bgf=F(1,50)
    wa=accumulation(T,Baf,Ha,Ca)
    wg=accumulation(T,Bgf,Hg,Cg)
    old=witness_certificate()["oscillatory_witness"]
    need_a=F(old["delivered_sample_necessary_accel_K_lower_mps"])
    need_g=F(old["delivered_sample_necessary_gyro_K_lower_rad"])
    return {
      "qualification":"OU3_CANDIDATE_JOINT_30S_60S_FAST_V1",
      "marine":{"T_E_s":"30","theta_E_deg":"2","T_P_s":"30","P_E_m":"0.03"},
      "fast":{"H_a_s":"60","C_a_mps":"0.05","H_g_s":"60","C_g_rad":"0.002"},
      "fast_30s_accumulation_cap_accel_mps":str(wa),
      "fast_30s_accumulation_cap_gyro_rad":str(wg),
      "old_oscillatory_witness_required_accel_cap_mps":str(need_a),
      "old_oscillatory_witness_required_gyro_cap_rad":str(need_g),
      "old_oscillatory_witness_excluded_by_accel":Ca<need_a,
      "old_oscillatory_witness_excluded_by_gyro":Cg<need_g,
      "old_oscillatory_witness_jointly_excluded":Ca<need_a or Cg<need_g,
      "full_same_history_gauge_breaking_closed":False,
      "paired_outer_bridge_closed":False,
      "hardware_qualification_pending":True,
      "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
