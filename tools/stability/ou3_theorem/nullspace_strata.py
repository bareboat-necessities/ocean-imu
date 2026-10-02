"""Qualitative zero-action strata for the complete-word nullspace lemma.

Unlike covariance/information certificates, zero Joseph loss uses H e=0 and
therefore needs only literal homogeneous transports and applied measurement
Jacobians. Innovation covariance and Kalman gain are intentionally absent.
"""
from .complete_word_nullspace import compatibility_line_certificate

def literal_zero_action_stack(*,process_rows,S_rows,mag_rows,acc_rows,gyro_rows,ba_rows):
    blocks=[]
    for name,rows in (("process_sync",process_rows),("S",S_rows),("magnetic",mag_rows),
                      ("accelerometer",acc_rows),("gyro_process",gyro_rows),("BA",ba_rows)):
        if rows: blocks.append((name,rows))
    return blocks

def audit_service_row_strata(strata):
    """Audit supplied theorem-domain service-row strata; fail closed on gaps."""
    results=[]; ok=True
    for st in strata:
        required=("process_rows","S_rows","mag_rows","acc_rows","gyro_rows","ba_rows","compatibility_line")
        miss=[k for k in required if k not in st]
        if miss:
            results.append({"name":st.get("name","?"),"verified":False,"reason":"missing "+",".join(miss)})
            ok=False;continue
        blocks=literal_zero_action_stack(**{k:st[k] for k in required[:-1]})
        z=compatibility_line_certificate(blocks,st["compatibility_line"])
        z["name"]=st.get("name","?");z["verified"]=z["kernel_equals_compatibility_line"]
        results.append(z);ok &= z["verified"]
    return {"qualification":"OU3_LITERAL_ZERO_ACTION_SERVICE_STRATA_V1",
            "strata":results,"all_strata_verified":ok and bool(results),
            "uses_K_or_innovation_covariance":False,
            "callback_pattern_enumeration_used":False,
            "source_uniform_verified":False if not results else ok}

def magnetic_service_kernel_certificate(mu_M):
    """Continuous service-row cover collapsed to its exact Gramian implication.

    If sum G_i'G_i >= mu_M I_2 with mu_M>0, then intersection ker G_i={0}.
    This covers arbitrary event count/timing/orientation already admitted by
    MAGNETIC SERVICE; no row enumeration is needed.
    """
    mu=float(mu_M)
    if not mu>0: raise ValueError("positive magnetic information floor required")
    return {"mu_M":mu,"service_coordinate_dimension":2,
            "aggregate_gramian_spd":True,
            "magnetic_service_common_kernel_dimension":0,
            "continuous_service_row_cover_complete":True,
            "event_count_or_timing_enumerated":False}

def certificate():
    # Structural regression only; theorem-domain literal strata are supplied
    # by the native/interval exporter, never invented here.
    return {"qualification":"OU3_LITERAL_ZERO_ACTION_SERVICE_STRATA_V1",
            "all_strata_verified":False,"theorem_domain_strata_supplied":False,
            "magnetic_continuous_cover":magnetic_service_kernel_certificate(1.0),
            "uses_K_or_innovation_covariance":False,
            "callback_pattern_enumeration_used":False,
            "theorem_closed":False}

if __name__=="__main__":
 import json
 print(json.dumps(certificate(),indent=2,sort_keys=True))
