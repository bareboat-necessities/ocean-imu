"""Exhaustive every-IMU-root physical +/- service replay on one literal stream.

FINITE DIAGNOSTIC ONLY.  This carries 200 simultaneously overlapping one-second
probe words through the exact exported F, Joseph K/H, reset G and magnetic
innovation S.  It does not replace outward joint-cell containment.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
from .planar_service_stream import records
from .planar_service_audit import correction_record,reset_matrix

def initial_phi(k:int)->np.ndarray:
    t=k*.005; pitch=.02*math.sin(math.pi*t/10)
    cp,sp=math.cos(-pitch),math.sin(-pitch)
    U=np.array([[cp,0,sp],[0,1,0],[-sp,0,cp]],float)
    sn=400/40001;cs=39999/40001
    p=np.zeros((21,4))
    for sign in range(2):
        d=U@np.array([0, -sn if sign else sn, cs])
        p[:3,2*sign]=d;p[3:6,2*sign+1]=.02*d
    return p

def exhaustive(path:Path):
    active={} ; vals=[]; sample_first=None;sample_last=None
    for kind,k,a in records(path):
        if kind==1:
            F=np.zeros((21,21)); # expand compressed F is full parity matrix
            from .planar_service_stream import expand
            F=expand(a[225:450])
            for q in active.values(): q["phi"]=F@q["phi"]
        elif kind in (2,3,4):
            P,H,R,S,K,PCt,res=correction_record(a); A=np.eye(21)-K@H
            if kind==3:
                L=np.linalg.cholesky((S+S.T)/2)
                for q in active.values():
                    Y=H@q["phi"]; W=np.linalg.solve(L,Y)
                    q["info"]+=W.T@W;q["mags"]+=1
            for q in active.values():q["phi"]=A@q["phi"]
        elif kind==5:
            G=reset_matrix(a[225:228])
            for q in active.values():q["phi"]=G@q["phi"]
        elif kind==7:
            sample_first=k if sample_first is None else sample_first;sample_last=k
            # Complete words ending at this sample, then seed this exact root.
            done=[root for root in active if k-root==200]
            for root in done:
                q=active.pop(root)
                ev=[]
                for j in (0,2):ev.append(float(np.linalg.eigvalsh((q["info"][j:j+2,j:j+2]+q["info"][j:j+2,j:j+2].T)/2).min()))
                vals.append((root,min(ev),q["mags"]))
            active[k]={"phi":initial_phi(k),"info":np.zeros((4,4)),"mags":0}
    if not vals:raise ValueError("no complete one-second roots")
    worst=min(vals,key=lambda x:x[1])
    return {"qualification":"OU3_PLANAR_EVERY_ROOT_SERVICE_DIAGNOSTIC_V1","result_type":"FINITE DIAGNOSTIC ONLY",
      "stream_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"sample_range":[sample_first,sample_last],
      "complete_one_second_roots":len(vals),"root_stride_samples":1,"worst_root_sample":worst[0],
      "minimum_physical_plus_minus_2x2_information":worst[1],"worst_root_accepted_magnetic_events":worst[2],
      "mu_M":1.0,"finite_margin_above_mu_M":worst[1]-1.0,
      "joint_cell_forward_invariant":False,"information_cell_perturbation_certified":False,
      "all_time_magnetic_service_verified":False,"theorem_closed":False}

def main():
 p=argparse.ArgumentParser();p.add_argument("--stream",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();o=exhaustive(a.stream);a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
